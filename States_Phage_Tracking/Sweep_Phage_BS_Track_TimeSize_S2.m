% Thomas C. Day
% Script to find the timescale and size scales at which first infection
% with phage occurs. Run the simulation until all groups of S2 have been
% infected.

% INPUTS ------------------------------------------------------------------
Msims   = 1e1;                          % number of simulations to run
FigViz  = 2;                            % T/F show the figures of the simulations
Nrounds = 218;                          % number of rounds of selection

% Compile traits:
alpha   = 2.7;                          % scaling exponent for encounter rate vs. size
phi     = 0.5;                          % packing fraction to turn trait -> size (biomass)
r0_F    = 0.5;                          % size of food particle
r0_P    = 0.1;                          % size of phage particle
delta   = 3e-2;                         % decay rate [#/doubling time]

% Phage parameters:
Track_Infection = 1;                    % T/F track the sfirst infection time of each multicellular individual
Infect_Prob = 1.0;                      % chance that the encountered phage infects the individual
Phage_Death_Prob = 1;                   % probability that an infection from a phage particle kills the individual
Phage_Burst_Size = [0,50,100,200];                 % number of viable phage per individual death per unit biomass

% Other parameters:
Kappa   = 1e4;                          % carrying capacity of the population
Mu      = [1; 1];                       % maximum growth rate of the two strategies [# divisions/doubling time]
K_m     = [0.5; 0.5];                   % half-max number of encounters/capita to get max growth rate, Monod dynamics
Delta   = [-delta; -delta];             % death rates for each strategy, constant
DeathThresh  = [0.90; 0.90];            % the biomass threshold at which an individual dies
DivideThresh = [2.00; Inf];             % the biomass threshold at which an individual of a particular strategy divides
% -------------------------------------------------------------------------

% Sweeping parameters:
% G0_F    = logspace(-2,0,5);             % encounter rate with food [# encounters/doubling time]
G0_F    = 1e-1;
G0_P    = logspace(-5,-1,5);             % encounter rate with infectious phage
TST_INF = struct('Time',[],'Size',[]);  % array enumerating properties of groups at first infection

%% Run each simulation ----------------------------------------------------
for bb = 1:length(Phage_Burst_Size)
    disp('------------------------');
    disp('Burst Size');
    disp('------------------------');
    BS = Phage_Burst_Size(bb);
    figure(1); clf;
    figure(2); clf;
    for pp = 1:length(G0_P)
        disp(pp);
        g0_f = G0_F;
        g0_p = G0_P(pp);
            
        for mm = 1:Msims
            
            % Run a simulation with fluctuations:
            [n_frac{bb,pp,mm}, b_frac{bb,pp,mm}, ~, ~, ~, NumInd, c_phage{bb,pp,mm}, tst_inf] = RUN_SINGLE_SIMULATION_STATES_PHAGE_TIMESIZETRACK(Nrounds, g0_f, g0_p, r0_F, r0_P,...
                alpha, phi, 0, Mu, K_m, Delta, Kappa, DivideThresh, DeathThresh, ...
                Track_Infection, Infect_Prob, Phage_Death_Prob, BS, FigViz);
                % n_track: the population list, where each entry corresponds to a
                    % member of the population. 1=resident, 2=mutant.
                % b_track: biomass list, each entry corresponds to one member of
                    % the population. The value is the number of cells for that
                    % individual.
            
            if FigViz
                % Show the phage concentration over time:
                figure(1);
                hold on; box on; set(gca,'linewidth',1);
                plot(c_phage{bb,pp,mm}./c_phage{bb,pp,mm}(1), '-','linewidth',2,'color','k');
                xlabel('Time [No. Divisions]');
                ylabel('Concentration of infectious phage');
                drawnow;
            
                % show the biomass fraction over time:
                figure(2);
                hold on; box on; set(gca,'linewidth',1);
                plot(n_frac{bb,pp,mm}, '-','linewidth',1,'color','k');
                plot(b_frac{bb,pp,mm},'-','linewidth',1,'color','r');
                title(num2str(g0_p));
                drawnow;
                
            end
                    
        end
        
        if Track_Infection
            ix = find(tst_inf(:,1) > 0); % find all instances where at least one infection occurred
            t_inf = tst_inf(ix,1);
            s_inf = tst_inf(ix,2);
            
            if FigViz == 3
                figure;
                tiledlayout('flow');
                nexttile;
                hold on; box on; set(gca,'linewidth',1);
                histogram(t_inf);
                xlabel('Time of First Infection');
                ylabel('Counts');
            
                nexttile;
                hold on; box on; set(gca,'linewidth',1);
                histogram(s_inf);
                xlabel('Biomass at First Infection');
                ylabel('Counts');
            end
        
        end

        % % Record the time and size distribution of infection:
        % TST_INF(ff,pp).Time = t_inf;
        % TST_INF(ff,pp).Size = s_inf;

    end
end

save('Phage_Sweep_BurstSize_PhageConc.mat');