# Multi-State Receptor System 
This repository contains Python code to calculate the trajectory mutual information rate between the ligand trajectory and the trajectory of the multi-state receptor as explained in the Paper "State- versus Reaction-Based Information Processing in Biochemical Networks" by Moor et al. 

## Content
- `multistatereceptor.py` contains the main code to calculate the trajectory mutual information rate for this case
- The file `General_Algebra.py` contains all the functions that are needed to use the main codes

For further information, contact Anne-Lena Moor (annemoor96@gmail.com). 

## Installation 
The provided code was written using Python v3.8.5 and uses the following libraries:
- [Numpy v1.19.2](https://www.numpy.org/)
- [Matplotlib v3.3.2](https://matplotlib.org/)


For installing Python and the required packages, one can use [Anaconda](https://www.anaconda.com/products/distribution#windows). Anaconda is a general package manager that contains all the for this code required packages. Instructions to install Python can be found [here](https://jupyter.readthedocs.io/en/latest/install.html). 

## Reproducing the figures of the paper 
The python script `multistatereceptor.py` contains the source code to calculate the path mutual information for the second case study in our paper. The scripts are commented in order to give a better understanding for the pipeline. Executing the pipeline in the source codes will lead to a reproduction of Figure 2 (b) in the main text. 

The code first calculates the mutual information rate between ligand and multi-state receptor depending on the number of receptor states and saves the resulting output file. For plotting, it loads the data contained in the output file, hence, it is required that the directory of code and file match. The rate of the single state receptor is calculated analytically and plotted as well.