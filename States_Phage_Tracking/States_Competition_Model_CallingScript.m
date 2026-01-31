% Thomas C. Day
% Script to find if fluctuations in resources alone can lead to
% trade-offs that either facilitate coexistence, or lead to the larger
% group out-competing the smaller group even when disadvantaged by mean.

% INPUTS ------------------------------------------------------------------
Msims   = 1e0;                          % number of simulations to run
FigViz  = 2;                            % T/F show the figures of the simulations
Nrounds = 3e2;                          % number of rounds of selection

% Compile traits for the two competing strategies:
alpha   = 2.7;                          % scaling exponent for encounter rate vs. size
phi     = 0.5;                          % packing fraction to turn trait -> size (biomass)
G0_F    = 1.0;                          % encounter rate for the resident strategy [# encounters/doubling time]
G0_P    = 0.01;                        % encounter rate with infectious phage
r0_F    = 0.5;                          % size of food particle
r0_P    = 0.1;                          % size of phage particle
delta   = 0;                            % decay rate [#/doubling time]

% Phage parameters:
Track_Infection = 1;                    % T/F track the simulation until the first infection, then break
Infect_Prob = 1.0;                      % chance that the encountered phage infects the individual
Phage_Death_Prob = 0;                   % probability that an infection from a phage particle kills the individual
Phage_Burst_Size = 100;                 % number of viable phage per individual death per unit biomass

% Other parameters:
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
    [n_frac, b_frac, n_track, b_track, p_track, NumInd, c_phage, tst_inf] = RUN_SINGLE_SIMULATION_STATES_PHAGE(Nrounds, G0_F, G0_P, r0_F, r0_P,...
        alpha, phi, 0, Mu, K_m, Delta, Kappa, DivideThresh, DeathThresh, ...
        Track_Infection, Infect_Prob, Phage_Death_Prob, Phage_Burst_Size, FigViz);
        % n_track: the population list, where each entry corresponds to a
            % member of the population. 1=resident, 2=mutant.
        % b_track: biomass list, each entry corresponds to one member of
            % the population. The value is the number of cells for that
            % individual.

    % Track the phage infections over time:
    fi1 = zeros(length(p_track),1); fi2 = fi1;
    for tt = 1:length(p_track)
        n1 = find(n_track{tt} == 1);
        n2 = find(n_track{tt} == 2);
        ixi1 = find(p_track{tt}(n1) > 0);
        ixi2 = find(p_track{tt}(n2) > 0);
        fi1(tt) = length(ixi1)/length(n1);
        fi2(tt) = length(ixi2)/length(n2);
    end

    
    if FigViz
        % Show the phage concentration over time:
        figure('units','centimeters','position',[18,3,5,5]);
        hold on; box on; set(gca,'linewidth',1);
        plot(c_phage./c_phage(1), '-','linewidth',2,'color','k');
        xlabel('Time [No. Divisions]');
        ylabel('Concentration of infectious phage');
    
        % Show the infected fraction over time:
        figure('units','centimeters','position',[3,13,10,10]);
        hold on; box on; set(gca,'linewidth',1);
        plot(fi1, '-','linewidth',2,'color',[.2,.4,.8]);
        plot(fi2, '-','linewidth',2,'color',[.8,.2,.2]);
        xlabel('Time [No. Divisions]');
        ylabel('Fraction of Infected Population');

        % Show the phage infection histogram at the end:
        figure('units','centimeters','position',[13,13,10,5]);
        n1 = find(n_track{end} == 1);
        n2 = find(n_track{end} == 2);
        M = -0.5:1:max(p_track{end});
        hold on; box on; set(gca,'linewidth',1);
        histogram(p_track{end}(n1),M,'facealpha',0.5,'normalization','pdf','edgecolor','none','facecolor',[.2,.4,.8]);
        histogram(p_track{end}(n2),M,'facealpha',0.5,'normalization','pdf','edgecolor','none','facecolor',[.8,.2,.2]);
        xlim([-0.5,max(p_track{end})]);
        xlabel('Number of infections');
        ylabel('PDF')
    end
            
end

%% S
if Track_Infection
    ix = find(tst_inf(:,1) > 0); % find all instances where at least one infection occurred
    t_inf = tst_inf(ix,1);
    s_inf = tst_inf(ix,2);
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