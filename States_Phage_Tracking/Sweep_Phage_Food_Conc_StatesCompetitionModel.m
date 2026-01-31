% Thomas C. Day
% Script to find if fluctuations in resources alone can lead to
% trade-offs that either facilitate coexistence, or lead to the larger
% group out-competing the smaller group even when disadvantaged by mean, as
% we sweep phage concentrations and burst sizes.

% This simulation's STOP condition is when state S2 goes extinct, or when
% the simulation has run for 72 hrs.

% INPUTS ------------------------------------------------------------------
Msims   = 10;                          % number of simulations to run
FigViz  = 0;                            % T/F show the figures of the simulations
Nrounds = 1e3;                          % number of rounds of selection
SaveData = 1;
save_folder = 'E:\2025\01_ProjectEncounter\202511XX_SimsRevisions\Phage_Food_Sweeps_Sims\';

% Compile traits:
alpha   = 2.7;                          % scaling exponent for encounter rate vs. size
phi     = 0.5;                          % packing fraction to turn trait -> size (biomass)
r0_F    = 0.5;                          % size of food particle [microns]
r0_P    = 0.1;                          % size of phage particle [microns]
delta   = 3e-2;                         % decay rate [#/doubling time]

% Phage parameters:
Track_Infection = 0;                    % T/F track the first infection time of each multicellular individual
Infect_Prob = 1.0;                      % chance that the encountered phage infects the individual
Phage_Death_Prob = 0.5;                 % probability that an infection from a phage particle kills the individual

% Other parameters:
Kappa   = 1e4;                          % carrying capacity of the population
Mu      = [1; 1];                       % maximum growth rate of the two strategies [# divisions/doubling time]
K_m     = [0.5; 0.5];                   % half-max number of encounters/capita to get max growth rate, Monod dynamics
Delta   = [-delta; -delta];             % death rates for each strategy, constant
DeathThresh  = [0.90; 0.90];            % the biomass threshold at which an individual dies
DivideThresh = [2.00; Inf];             % the biomass threshold at which an individual of a particular strategy divides

% -------------------------------------------------------------------------
% Sweeping parameters:
G0_F    = [.01,.03,.10,.20,.30,.50,1];             % encounter rate with food [# encounters/doubling time]
% G0_F    = 0.1;
G0_P    = [1e-5,1e-4,1e-3,3e-3,1e-2,3e-2,.1];             % encounter rate with infectious phage
% G0_P    = 3e-3;
% G0_P = 0;
BurstSize = [0, 50, 100, 200];     % phage burst size

% Pre-allocate measurements:
n_frac = cell(length(G0_F), length(G0_P), length(BurstSize), Msims);
b_frac = n_frac;
NumInd = n_frac;
c_phage = n_frac;

%% Run each simulation ----------------------------------------------------
for ff = 1:length(G0_F)
    disp('------------------------');
    disp(ff);
    for pp = 1:length(G0_P)
        disp(pp);
        g0_f = G0_F(ff);
        g0_p = G0_P(pp);

        for bb = 1:length(BurstSize)
            Phage_Burst_Size = BurstSize(bb);
            
            for mm = 1:Msims
                
                % Run a simulation with fluctuations:
                [n_frac{ff,pp,bb,mm}, b_frac{ff,pp,bb,mm}, n_track, b_track, NumInd{ff,pp,bb,mm}, c_phage{ff,pp,bb,mm}, tst_inf] = RUN_SINGLE_SIMULATION_STATES_PHAGE(Nrounds, g0_f, g0_p, r0_F, r0_P,...
                    alpha, phi, 0, Mu, K_m, Delta, Kappa, DivideThresh, DeathThresh, ...
                    Track_Infection, Infect_Prob, Phage_Death_Prob, Phage_Burst_Size, FigViz);
                    % n_track: the population list, where each entry corresponds to a
                        % member of the population. 1=resident, 2=mutant.
                    % b_track: biomass list, each entry corresponds to one member of
                        % the population. The value is the number of cells for that
                        % individual.
    
                   
    
                if FigViz == 2
                    figure('units','centimeters','position',[3,3,10,6]);
                    hold on; box on; set(gca,'linewidth',1);
                    plot(n_frac{ff,pp,bb,mm}, '-','linewidth',1,'color','k');
                    plot(b_frac{ff,pp,bb,mm}, '-','linewidth',1,'color',[.8,.2,.2]);
                    xlabel('Time (Rounds)');
                    ylabel('Fraction');
                    ylim([0,1]);
                end
                
                if FigViz == 3
                    % Show the phage concentration over time:
                    figure('units','centimeters','position',[18,3,5,5]);
                    hold on; box on; set(gca,'linewidth',1);
                    plot(c_phage{ff,pp,bb,mm}./c_phage{ff,pp,bb,mm}(1), '-','linewidth',2,'color','k');
                    xlabel('Time [No. Divisions]');
                    ylabel('Concentration of infectious phage');
                
                    % Show the infected fraction over time:
                    %{
                    figure('units','centimeters','position',[3,13,10,10]);
                    hold on; box on; set(gca,'linewidth',1);
                    plot(fi1, '-','linewidth',2,'color',[.2,.4,.8]);
                    plot(fi2, '-','linewidth',2,'color',[.8,.2,.2]);
                    xlabel('Time [No. Divisions]');
                    ylabel('Fraction of Infected Population');
            
                    %}
    
                    % Show the phage infection histogram at the end:
                    %{
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
                    %}
                end
                        
            end
            disp('Simulation Set DONE');
            
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
    
                % Record the time and size distribution of infection:
                TST_INF(ff,pp).Time = t_inf;
                TST_INF(ff,pp).Size = s_inf;
            
            end

            if FigViz == 4
                figure('units','centimeters','position',[3,3,10,6]);
                hold on; box on; set(gca,'linewidth',1);
                plot(n_frac{ff,pp,bb,1}, '-','linewidth',1,'color','k');
                plot(b_frac{ff,pp,bb,1}, '-','linewidth',1,'color',[.8,.2,.2]);
                xlabel('Time (Rounds)');
                ylabel('Fraction');
                title('First sim from this set');
                ylim([0,1]); drawnow;
            end

        end
    end
end
disp('SIMULATION COMPLETE');
disp('-------------------');

if SaveData
    save([save_folder, 'sims_PHAGE_FS_',date,'.mat']);
end


%% Load data, and calculate the time to extinction and the max biomass, also maybe the integrated biomass:
save_folder = 'E:\2025\01_ProjectEncounter\202511XX_SimsRevisions\Phage_Food_Sweeps_Sims\';
load([save_folder,'sims_PHAGE_FS_12-Dec-2025.mat']);

TimeToExtinction    = zeros(size(n_frac));
MaxBiomass          = zeros(size(n_frac));
IntegratedBiomass   = zeros(size(n_frac));
YNSimsSurviving     = zeros(size(n_frac));

% Make measurements:
for ff = 1:length(G0_F)
    for pp = 1:length(G0_P)
        for bb = 1:length(BurstSize)
            for mm = 1:Msims
                % Time to extinction:
                if n_frac{ff,pp,bb,mm}(end) == 0
                    TimeToExtinction(ff,pp,bb,mm) = length(n_frac{ff,pp,bb,mm});
                else
                    TimeToExtinction(ff,pp,bb,mm) = 218;
                end
            
                % Max Biomass Fraction:
                MaxBiomass(ff,pp,bb,mm) = max(b_frac{ff,pp,bb,mm});

                % Integrated Biomass
                IntegratedBiomass(ff,pp,bb,mm) = sum(b_frac{ff,pp,bb,mm});

                % Check the fraction of simulations that are still
                % surviving after T duration:
                T_choice = 72; % rounds, equivalent to 24 hrs
                if length(n_frac{ff,pp,bb,mm}) >= T_choice
                    if n_frac{ff,pp,bb,mm}(T_choice) == 0
                        YNSimsSurviving(ff,pp,bb,mm) = 0;
                    else
                        YNSimsSurviving(ff,pp,bb,mm) = 1;
                    end
                else
                    if n_frac{ff,pp,bb,mm}(end) == 0
                        YNSimsSurviving(ff,pp,bb,mm) = 0;
                    else
                        YNSimsSurviving(ff,pp,bb,mm) = 1;
                    end
                end
            end
        end
    end
end

% Pick a food concentration, vary PFU and burst size:
for ff = 3%1:length(G0_F)
    ii = 0;
    figure('units','centimeters','position',[3,3,17,17]);
    % axs = 
    for pp = 1:length(G0_P)
        row = floor(pp/7);
        for bb = 1:length(BurstSize)
            ii = ii + 1;
            subplot(7,4,ii); hold on; box on; set(gca,'linewidth',1);
            for mm = 1:Msims
                % Show the curves:
                plot(n_frac{ff,pp,bb,mm},'-','linewidth',1,'color',[.5,.5,.5]);
                plot(b_frac{ff,pp,bb,mm},'-','linewidth',1,'color',[.8,.2,.2]);
            end
            % xlabel('Time (rounds)');
            xlim([0,218]);
            ylim([0,1]);
            set(gca,'fontsize',7);
        end
    end
end

% Pick a burst size, show the varying food and PFU conc:
bix = 3;
figure('units','centimeters','position',[3,3,17,17]);
ii = 0;
for pp = 1:length(G0_P)
    for ff = 1:length(G0_F) 
        ii = ii + 1;           
        subplot(7,7,ii); hold on; box on; set(gca,'linewidth',1);
        for mm = 1:Msims
            % Show the curves:
            plot(n_frac{ff,pp,bix,mm},'-','linewidth',1,'color',[.5,.5,.5]);
            plot(b_frac{ff,pp,bix,mm},'-','linewidth',1,'color',[.8,.2,.2]);
        end
        % xlabel('Time (rounds)');
        xlim([0,218]);
        ylim([0,1]);
        set(gca,'fontsize',7);
    end
end

% Get sim averages:
Mean_Tau = mean(TimeToExtinction, 4); Err_Tau = std(TimeToExtinction, 0, 4);
Mean_Bmax = mean(MaxBiomass, 4); Err_Bmax = std(MaxBiomass, 0, 4);
Mean_Bint = mean(IntegratedBiomass, 4); Err_Bint = std(IntegratedBiomass, 0, 4);
FractionSimsSurviving = sum(YNSimsSurviving,4)./10;

% Show time to extinction for constant food and constant burst size:
figure('units','centimeters','position',[3,3,9,12]);

subplot(3,2,1);
hold on; box on; set(gca,'linewidth',1);
imagesc(1:4, G0_P, squeeze(Mean_Tau(3,:,:)));
xlabel('Burst Size');
ylabel('PFU');
set(gca,'layer','top');
xlim([.5,4.5]);
set(gca,'yscale','log');
ylim([5e-6,2e-1])
yticks([1e-5,1e-3, 1e-1]);
xticks([1,2,3,4]);
xticklabels({'0','50','100','200'});

subplot(3,2,2);
hold on; box on; set(gca,'linewidth',1);
imagesc(1:7, G0_P, squeeze(Mean_Tau(:,:,bix))');
xlabel('Food');
ylabel('PFU');
set(gca,'layer','top');
set(gca,'yscale','log');
ylim([5e-6,2e-1]);
yticks([1e-5,1e-3, 1e-1]);
xlim([.5,7.5]);
xticks([1,3,5,7]);
xticklabels({'0.01','0.1','0.3','1'})

subplot(3,2,3);
hold on; box on; set(gca,'linewidth',1);
imagesc(1:4, G0_P, squeeze(FractionSimsSurviving(3,:,:)));
xlabel('Burst Size');
ylabel('PFU');
set(gca,'layer','top');
xlim([.5,4.5]);
set(gca,'yscale','log');
ylim([5e-6,2e-1])
yticks([1e-5,1e-3, 1e-1]);
xticks([1,2,3,4]);
xticklabels({'0','50','100','200'});

subplot(3,2,4);
hold on; box on; set(gca,'linewidth',1);
imagesc(1:7, G0_P, FractionSimsSurviving(:,:,bix)');
xlabel('Food');
ylabel('PFU');
set(gca,'layer','top');
set(gca,'yscale','log');
ylim([5e-6,2e-1]);
yticks([1e-5,1e-3, 1e-1]);
xlim([.5,7.5]);
xticks([1,3,5,7]);
xticklabels({'0.01','0.1','0.3','1'})

subplot(3,2,5);
hold on; box on; set(gca,'linewidth',1);
imagesc(1:4, G0_P, squeeze(Mean_Bmax(3,:,:)));
xlabel('Burst Size');
ylabel('PFU');
set(gca,'layer','top');
xlim([.5,4.5]);
set(gca,'yscale','log');
ylim([5e-6,2e-1])
yticks([1e-5,1e-3, 1e-1]);
xticks([1,2,3,4]);
xticklabels({'0','50','100','200'});

subplot(3,2,6);
hold on; box on; set(gca,'linewidth',1);
imagesc(1:7, G0_P, Mean_Bmax(:,:,bix)');
xlabel('Food');
ylabel('PFU');
set(gca,'layer','top');
set(gca,'yscale','log');
ylim([5e-6,2e-1]);
yticks([1e-5,1e-3, 1e-1]);
xlim([.5,7.5]);
xticks([1,3,5,7]);
xticklabels({'0.01','0.1','0.3','1'})