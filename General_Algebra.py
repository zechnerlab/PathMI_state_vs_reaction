#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 13:04:57 2024

@author: moor
"""

import sympy
import numpy as np
from scipy.integrate import ode

class Reaction:
    def __init__(self, name, reactantStoichiometry, productStoichiometry, rate, propensity):
        assert isinstance(name, str), "name must be a string"
        assert isinstance(reactantStoichiometry, list), "reactantStoichiometry must be a list"
        assert isinstance(productStoichiometry, list), "productStoichiometry must be a list"
        self.name = name
        self.reactantStoichiometry = np.array(reactantStoichiometry)  # Positive sign
        self.productStoichiometry = np.array(productStoichiometry)    # Positive sign
        self.stoichiometry = self.productStoichiometry - self.reactantStoichiometry
        self.rate = rate
        self.propensity = propensity
    
    def isChanging(self, index):
        return self.stoichiometry[index] != 0
    
    def __str__(self):
        return f"{self.name}: {self.reactantStoichiometry} -> {self.productStoichiometry} @ {self.rate} ({self.propensity})"
    
    def __repr__(self):
        return self.__str__()

class System:
    def __init__(self, numVariables, observedIndices, symbolic=False):
        assert isinstance(numVariables, int) and numVariables > 0, "numVariables must be a positive integer"
        assert isinstance(observedIndices, list) and all([ isinstance(i, int) and i >= 0 and i < numVariables for i in observedIndices ]), "observedIndices must be a list of integers between 0 and numVariables-1"
        assert isinstance(symbolic, bool), "symbolic must be a boolean"
        self.reactions = []
        self.lnaignored = []
        self.numVariables = numVariables
        self.observedIndices = observedIndices  # Indices of the observed variables
        self.hiddenIndices = [ i for i in range(numVariables) if i not in observedIndices ]
        self.numHiddenVariables = numVariables - len(observedIndices)
        self.symbolic = symbolic
    
    def getInitialK(self):
        return np.zeros((self.numHiddenVariables, self.numHiddenVariables))

    def isObserved(self, index):
        return index in self.observedIndices
    
    def isHidden(self, index):
        return not self.isObserved(index)
    
    def isSymbolic(self):
        return self.symbolic

    def addReaction(self, reaction, islnaignored=False):
        self.reactions.append(reaction)
        self.lnaignored.append(islnaignored)
        return self
    
    def addReactions(self, reactions, lnaignorednames=[]):
        for reaction in reactions:
            name = reaction.name
            if name in lnaignorednames:
                self.addReaction(reaction,True)
            else:
                self.addReaction(reaction,False)
        return self

    def isExclusivelyHidden(self, reaction):
        ret = True
        for i in self.observedIndices:
            if reaction.isChanging(i):
                ret = False
                break
        return ret
    
    def generateHiddenStoichiometryMatrix(self):    # S
        stoichiometryMatrix = np.zeros( (self.numHiddenVariables, len(self.reactions)) )
        for i in range(self.numHiddenVariables):
            for j in range(len(self.reactions)):
                hi = self.hiddenIndices[i]
                stoichiometryMatrix[i][j] = self.reactions[j].stoichiometry[hi]
        return stoichiometryMatrix
    
    def generateHiddenReactantStoichiometryMatrix(self):    # E
        reactantStoichiometryMatrix = np.zeros((self.numHiddenVariables, len(self.reactions)))
        for i in range(self.numHiddenVariables):
            for j in range(len(self.reactions)):
                hi = self.hiddenIndices[i]
                reactantStoichiometryMatrix[i][j] = self.reactions[j].reactantStoichiometry[hi]
        return reactantStoichiometryMatrix
    
    def generateConstantDiag(self):   # c
        c = [ reaction.rate for reaction in self.reactions ]
        return np.diag(c)
    
    def _sqrt(self, x):
        if self.symbolic:
            return sympy.sqrt(x)
        else:
            return np.sqrt(x)
    
    def generatePropensityDiag(self):   # sigma
        a = [ 0.0 if flag else self._sqrt(reaction.propensity) for (reaction, flag) in zip(self.reactions, self.lnaignored) ]
        return np.diag(a)
    
    def generateHiddenPropensityDiag(self):  # sigmaY
        a = [ 0 if self.isExclusivelyHidden(reaction) else self._sqrt(reaction.propensity) for reaction in self.reactions ]
        return np.diag(a)
    
    def generateInverseHiddenPropensityDiag(self):  # sigmaYinv
        a = [ 0 if self.isExclusivelyHidden(reaction) else 1.0/self._sqrt(reaction.propensity) for reaction in self.reactions ]
        return np.diag(a)
    
    def reacStoichMatrix(self):
        return np.transpose(sympy.Matrix([r.reactantStoichiometry for r in self.reactions]))
    
    def prodStoichMatrix(self):
        return np.transpose(sympy.Matrix([r.productStoichiometry for r in self.reactions]))
    
    def rhs(self, K):
        S = self.generateHiddenStoichiometryMatrix()
        E = self.generateHiddenReactantStoichiometryMatrix()
        c = self.generateConstantDiag()
        sigma = self.generatePropensityDiag()
        sigmaY = self.generateHiddenPropensityDiag()
        sigmaYinv = self.generateInverseHiddenPropensityDiag()

        A = S @ c @ E.T @ K
        B = S @ sigma
        G = K @ E @ sigmaYinv @ c + S @ sigmaY

        dK = A + A.T + (B @ B.T) - (G @ G.T)

        return dK
    
    def latexRhs(self, K, filename):
        assert self.isSymbolic # We need a symbolic system in order to be able to print its LaTeX form
        dK = self.rhs(K)
        with open(filename, 'w') as f:
            for i in range(self.numHiddenVariables):
                for j in range(i+1):
                    lhs = "\\frac{\\mathrm{d}}{\\mathrm{d}t} K_{%d,%d}" %(i,j)
                    rhs = sympy.latex(dK[i,j].simplify())
                    f.write(lhs + " = " + rhs + " \n")
    

def MatrixEquation(t, cov0, system, m):
    cov0 = cov0.reshape((m, m))
    dydt = system.rhs(cov0)           
    dydt=dydt.flatten()
    return dydt

    
def calculateCovariance(const, x_eq, reac_system, obsindx, dim, t_max, dt, lna_noise=False):
    '''
    Parameters
    ----------
    const : list of rate constants
    x_eq : list of steady state values of all chemical species
    system : function that generates the reaction system object
    dim : integer, dimension of covariance matrix; number of hidden species
    t_max : float, maximal time of integration
    dt : float, time step of integration

    Returns
    -------
    sol_y : numpy.ndarray, time x dim x dim array containing the resulting covariances of the integration

    '''
    
    if lna_noise:
        system = reac_system(const, x_eq, obsindx, lna_noise)
    else:
        system = reac_system(const, x_eq, obsindx)

    
    cov0_y = np.zeros((dim, dim)) 

    cov0_y = cov0_y.flatten()
    
    # set integrator, initial value, and parameters
    sol_y0 = ode(MatrixEquation).set_integrator('vode', method='bdf', order=4)
    # sol_y0.set_initial_value(cov0_y,0).set_f_params(S_y, R_y, S_y_prime, sigma, sigma_y, sigma_y_inv, const_matrix, dim)
    sol_y0.set_initial_value(cov0_y,0).set_f_params(system, dim)
    
    # integrator
    sol_y = [cov0_y]
    
    # sol_full = [np.array(M)]
    while sol_y0.successful() and sol_y0.t < t_max:
        sol_y.append(sol_y0.integrate(sol_y0.t+dt))
        
    sol_y = np.array(sol_y)
    
    sol_y = sol_y.reshape((len(sol_y), dim, dim))
    return sol_y

def integrator(integrateSpecies, iniconds, const, dt, t_max):
    y0 = iniconds
    
    sols0 = ode(integrateSpecies).set_integrator('vode', method='bdf', order=4)
    sols0.set_initial_value(y0,0).set_f_params(const,)
    
    # integrator
    sols = [y0]
        
    while sols0.successful() and sols0.t < t_max:
        sols.append(sols0.integrate(sols0.t+dt))    
    sols = np.array(sols)
    return sols





