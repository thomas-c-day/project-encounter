% Thomas C. Day
% Gillespie algorithm calling script for testing, but here we care only
% about the physical encounters, there is no growth associated with the
% simulation.

% Saving and sweeping parameters:
Msims       = [1,1,1,3,10,30,100];                  % number of simulations of each size
INIT_BIO    = [1e1,3e1,1e2,3e2,1e3,3e3,1e4];        % # cells in a group
BURST_SIZE  = [0,100,200];                   % burst size
PHAGE_CONC  = [1e4,3e4,1e5,3e5,1e6,3e6,1e7];        % phage concentration

% Simulation parameters:
figviz              = 0;            % Y/N show any figures
K                   = 1e4;          % population size/biomass in each strategy
tmax                = 120*60*60;    % [sec] set a maximum simulation time
initial_watervol    = 1e-1;         % [mL]
food_conc           = 1e0;          % [#/mL]
p_death             = 1.0;          % probability of dying/lysing to infectious phage
biomass_err         = 0;            % standard deviation for number of cells in the group
save_sims           = 1;
% save_folder         = 'C:\Users\daythoma\Documents\00_USC\02_Data\GillespieEncounters\Sweeping_Size_BurstSize_PhageConc\';
save_folder         = 'E:\Shared drives\Schwartzman Lab\Data\tom\Sweeping_Size_BurstSize_PhageConc_EMP\';
% save_folder         = 'F:\2025\01_ProjectEncounter\202511XX_SimsRevisions\Gillespie_SweepPhageSizeBurstSize\';

%% Sweeping different parameters:
for nn = 1:length(INIT_BIO)    
    msims_this_size = Msims(nn);
    for bb = 1:length(BURST_SIZE)
        for pp = 1:length(PHAGE_CONC)
            % Progress bar:
            disp([num2str(nn),' ',num2str(bb),' ',num2str(pp)]);

            save_folder_addendum = ['N=',num2str(nn,'%02.f'),'_B=',num2str(bb,'%02.f'),'_P=',num2str(pp,'%02.f'),'/'];
            mkdir([save_folder, save_folder_addendum]);
    
            initial_biomass = INIT_BIO(nn);
            burst_size = BURST_SIZE(bb);
            phage_conc = PHAGE_CONC(pp);
    
            % -------------------------------------------------------------------------
            % Initialize 50-50 SC and MC
            %{
            state_0 = 2*ones(K,1);
            state_0(1:floor(length(state_0)/2)) = 1;
            ix_2 = find(state_0==2);
            biomass_0 = 1*ones(size(state_0)) + 0*randn(size(state_0));
            biomass_0(ix_2) = biomass_err*randn(size(ix_2)) + initial_biomass;
            %}
    
            % Initialize 50-50 SC and MC by biomass:
            %{
            n_groups = floor(K./initial_biomass);
            state_0 = ones(K+n_groups,1);
            state_0(K+1:end) = 2;
            biomass_0 = state_0;
            biomass_0(K+1:end) = initial_biomass;
            ix_1 = find(state_0 == 1);
            ix_2 = find(state_0 == 2);
            biomass_1 = sum(biomass_0(ix_1));
            biomass_2 = sum(biomass_0(ix_2));      
            %}
    
            % Initialize with just the MC groups:
            n_groups = floor(K./initial_biomass);
            state_0 = 2*ones(n_groups,1);
            biomass_0 = initial_biomass * ones(n_groups,1);
            
            
            % -------------------------------------------------------------------------
            for mm = 1:msims_this_size
                % Run a simulation:
                LOG = GILLESPIE_ALGORITHM_FOOD_PHAGE_PHYSICS(state_0, biomass_0, food_conc, ...
                    phage_conc, p_death, burst_size, initial_watervol, tmax, K, figviz);
                disp(num2str(LOG(end).Time./60));
    
                % Save sim:
                if save_sims
                    save_name = ['SIM_N=',num2str(nn),'_B=',num2str(bb),'_P=',num2str(pp),'_M=',num2str(mm,'%03.f'),'.mat'];
                    save([save_folder,save_folder_addendum, save_name],'LOG','food_conc','phage_conc',...
                        'initial_biomass','K','tmax','initial_watervol','p_death','burst_size');
                end
    
                % Show sim results:
                if figviz == 2
                    figure;
                    for tt = 1:length(LOG)
                        clf;
                        hold on;
                        yyaxis left;
                        plot(LOG(tt).BANK_FOOD);
                        yyaxis right;
                        plot(LOG(tt).BANK_PHAGE);
                        plot(LOG(tt).BANK_DEATH,'kx','linewidth',1);
                        title([LOG(tt).Time/60]);
                        drawnow; pause(0.05);
                    end
                elseif figviz == 3
                    figure; hold on;
                    ix_1 = find(LOG(end).STATE == 1);
                    histogram(LOG(end).BANK_DEATH(ix_1));
                    line([LOG(end).BANK_DEATH(end), LOG(end).BANK_DEATH(end)],[.9,100],'linestyle','--','color','k','linewidth',1);
                    xlabel('Time to infection [sec]');
                    ylabel('Counts');
                end
            end
           
            % Ending the sweeping loops
        end
    end
end