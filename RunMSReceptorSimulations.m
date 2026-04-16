%% RunMSReceptorSimulations.m
% This is the main script that numerically simulates the state-based MI
% rate for a multistate receptor. To calculate this quantity, we introduce
% a fast downstream speces X, which can be produced from all active states.
% We then calculate the reaction-based MI rate between the ligand and X,
% which for very large production rates of X converges the state-based MI
% rate between the ligand and the receptor (i.e., the fast production
% reaction mimics a noise-free observer of the receptor on and off
% states).The script automatically generates the linear noise approximation
% and the relevant system matrices needed to solve the conditional
% covariance equations. The code is specific to the particular system and does not
% generalize to arbitrary systems. Note that the execution of the script
% can take several minutes as the number of ODEs that need to be solved scales
% quadratically with the number of receptor states. For N=40, for instance,
% the ODE system has >1500 dimensions.
% Author: Christoph Zechner (czechner (( AT )) sissa.it)

clear;
close all;


% Set this to 1 to compute the mutual information. Set this to 0 if you
% only want to plot from the already existing simulations (stored in
% "simulation_ms_receptor.mat".
runSimulations = 1;

if (runSimulations==1)

    % Time T should be large enough to achieve steady state. Note that very
    % large T can slow down the ODE solver.
    T = 5000;
    grid = linspace(0, T, 300);

    NVec = [1:10, 15:5:40]; % Number of active receptor states
    RTot = 10; % Total number of receptor copies
    c1 = 1; % Birth rate of ligand
    c2 = 0.01; % Death rate of ligand
    kOn = 1; % Activation rate

    % Birth and death rates of the downstream species X. Those are chosen fast to
    % converge to the state-based limit. Note that the degradation beta is only
    % needed to make X converge to a steady state. This parameter does
    % otherwise not enter the mutual information in the reaction-based
    % formalism, i.e., its particular value is irrelevant.
    gamma = 6000000;
    beta = 500000;

    % define parameters used for symbolic computation
    params = {'c1', 'c2', 'kOn', 'mu', 'gamma', 'beta'};


    % Run throught different N
    for u=1:length(NVec)

        % Set number of receptor states
        N = NVec(u);

        % Rate at which the receptor runs through the circle. This needs to be multiplied by N
        % to achieve a constant average circle time. This is why it is set
        % inside the loop.
        mu = 1*N;

        % Create stoihiometry of the receptor model
        [Pre, Post, c, X0, h] = GenerateReceptorModel(N, c1, c2, kOn, mu, gamma, beta);
        S = Post - Pre;

        %Generate conditional moments with only X observed. This will generate
        %a file "LNAMomentsODE_0.m", which defines the relevant ODE, which is solved later.
        GenerateLNAODEFilter(S, h, params, 0);

        %Generate conditional moments with both X and ligand observed. This
        %will generate a file "LNAMomentsODE_1.m", which defines the relevant
        %ODE, which is solved later.
        [X, h, A, BBT, F, BBT_Y] = GenerateLNAODEFilter(S, h, params, 1);

        for k=0:1
            condL = k;

            numStates = length(X0);

            % total number of states contains the downstream species. The inactive
            % state is eliminated algebraically. The number of "hidden" states is
            % therefore given by numStates - 2 - condL, where condL refers to the
            % ligand that is either conditioned on or not (which adds / removes one dimension).
            numHiddenStates = numStates - 2 - condL;

            % initialize system to steady state
            R0 = (c2*mu*RTot)/(N*c1*kOn + c2*mu);
            X0(1) = c1/c2;
            X0(2) = (1-R0)*gamma/beta;
            X0(3) = R0;
            X0(4:end) = (RTot-R0)/N;

            % covariance is initialized to zero
            Cov0 = zeros(numHiddenStates, numHiddenStates);
            M0 = [X0; reshape(Cov0, numHiddenStates^2, 1)];

            % set ODE solver options
            opts = odeset;
            opts.RelTol = 1e-9;
            opts.AbsTol = 1e-9;

            % Solve conditional moments
            if (condL==0)
                [~, M] = ode15s(@LNAMomentsODE_0, grid, M0, opts, c1, c2, kOn, mu, gamma, beta);
            elseif (condL==1)
                [~, M] = ode15s(@LNAMomentsODE_1, grid, M0, opts, c1, c2, kOn, mu, gamma, beta);
            end

            M = M';

            % Mutual information needs the total variance of the active state,
            % which can be obtained from the conditional covariance via projection
            W = zeros(1, numHiddenStates);
            W(2-condL:end) = 1; %select all active states and add them up

            CovCond = reshape(M(numStates+1:end, end), numHiddenStates, numHiddenStates);

            % Total conditional variance of the sum of the active states.
            CondVar_Bound(k+1) = W*CovCond*W';

        end


        % The total conditional variances need to be weighted by the relevant noise
        % term. This is obtained directly from the equation generator.
        sigmaS = BBT_Y;
        sigmaS = subs(sigmaS, X, M(1:numStates, end));
        sigmaS = double(subs(sigmaS, 'gamma', gamma));

        % Expression for the mutual information rate. It is calculated as the
        % difference between the conditional variances times
        % a weighting factor.
        mi_ms = gamma^2/2*(CondVar_Bound(1) - CondVar_Bound(2))/sigmaS;
        % Analytical expression for reaction-based MI rate
        mi_rb = -c2/2 + c2/2*sqrt(2*RTot * mu/N * kOn / (c2*mu/N + c1*kOn) + 1);
        % Analytical expression for state-based MI rate
        mi_sb = -c2/2 + c2/2*sqrt(RTot  * mu/N * kOn / (c2*mu/N + c1*kOn) + 1);

        fprintf('\nR0=%f, l=%f\n', RTot-sum(X0(4:end)), X0(1));
        fprintf('State-based: %f (nat/s)\n', mi_sb);
        fprintf('Reaction-based: %f (nat/s)\n', mi_rb);
        fprintf('Multi-state (N=%d): %f (nat/s)\n\n', N, mi_ms);

        mi_ms_vec(u) = mi_ms;

    end

    % Store results
    save("simulation_ms_receptor.mat");
end

% Create figure
load("simulation_ms_receptor.mat");
close all;
plot(NVec, mi_ms_vec, '.', 'MarkerSize',12); hold on;
plot([0, 41], mi_rb*ones(1, 2), '--', 'LineWidth', 1.5, 'Color', [0.3, 0.6, 0.3]);
plot([0, 41], mi_sb*ones(1, 2), '--', 'LineWidth', 1.5, 'Color', [0.6, 0.6, 0.6]);
xlim([0, 41]);
ylim([0.010, 0.019]);
xlabel('Number of Receptor States', 'Interpreter','tex');
ylabel('Mutual Information Rate', 'Interpreter', 'none');
legend('multi-state receptor', 'reaction-based', 'state-based');