% Thomas C. Day
% Script to find if fluctuations in resources alone can lead to
% trade-offs that either facilitate coexistence, or lead to the larger
% group out-competing the smaller group even when disadvantaged by mean.

% INPUTS ------------------------------------------------------------------
Msims   = 1;                          % number of simulations to run
FigViz  = 2;                            % T/F show the figures of the simulations
Nrounds = 2e2;                          % number of rounds of selection

% Compile traits for the two competing strategies:
alpha   = 2.5;                          % scaling exponent for encounter rate vs. size
phi     = 0.3;                          % packing fraction to turn trait -> size (biomass)
G_0     = 0.01;                         % encounter rate for the resident strategy [# encounters/doubling time]
delta   = 3e-2;                         % decay rate [#/doubling time]
Kappa   = 1e5;                          % carrying capacity of the population
Mu      = [1; 1];                       % maximum growth rate of the two strategies [# divisions/doubling time]
K_m     = [0.5; 0.5];                   % half-max number of encounters/capita to get max growth rate, Monod dynamics
Delta   = [-delta; -delta];             % death rates for each strategy, constant
DeathThresh  = [0.90; 0.90];            % the biomass threshold at which an individual dies
DivideThresh = [2.00; Inf];             % the biomass threshold at which an individual of a particular strategy divides
% -------------------------------------------------------------------------

%% Run each simulation ----------------------------------------------------
for mm = 1:Msims
    
    % Run a simulation with fluctuations:
    [n_frac, b_frac, ~, ~, NumInd] = RUN_SINGLE_SIMULATION_STATES(Nrounds, G_0, alpha, phi, 0, Mu, K_m, Delta, Kappa, DivideThresh, DeathThresh, FigViz);
        % n_track: the population list, where each entry corresponds to a
            % member of the population. 1=resident, 2=mutant.
        % b_track: biomass list, each entry corresponds to one member of
            % the population. The value is the number of cells for that
            % individual.
            
end