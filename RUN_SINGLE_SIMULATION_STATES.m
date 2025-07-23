function [n_frac, b_frac, n_track, b_track, NumInd] = RUN_SINGLE_SIMULATION_STATES(Nrounds, G0, alpha, phi, Dct, Mu, K_m, Delta, Kappa, CellDivisionThreshold, CellDeathThreshold, FigViz)
    % Run one simulation of a competition between two different states
    % within one genotype over a number of selection rounds. The two
    % different states will represent a single-celled state and a
    % multi-celled state.

    % Initial conditions:
    n0 = [ones(1e3,1); 2*ones(1e3,1)];
        % State 1 is the single-celled strategy, State 2 is the
        % multi-celled strategy. Half the individuals are of State 2.

    % Notate the biomass of each individual:
    b0 = 0.10*randn(size(n0)) + 1; 
        % start off as each individual is single-celled, each cell is
        % slightly different in size

    % Variables to track:
    n_track{1,1} = n0;
    b_track{1,1} = b0;
    n_frac(1,1)  = length(find(n0==2))./length(n0);
    b_frac(1,1)  = sum(b0(find(n0==2)))/sum(b0);
    NumInd = length(b0);

    % Run competition experiment:
    for tt = 1:Nrounds
        
        % Progress bar:
        disp([num2str(tt), ' / ',num2str(Nrounds)]);

        % Checks to stop the simulation:
        if length(n_track{tt}) > Kappa
            % Carrying capacity is reached
            break;
        elseif length(n_track{tt}) <= 1
            % Extinction
            break;
        elseif n_frac(tt) > 0.99
            % Fixation
            break;
        elseif n_frac(tt) == 0
            % Other fixation
            break;
        else
            % Continue the simulation
            if FigViz == 1
                figure(2); clf; hold on;
                plot(n_track{tt});
                plot(b_track{tt});
                xlim([0,5000]);
                ylim([0,100]);
                title(tt);
                drawnow; pause(0.25);
            end
            
            % 0. Assign everybody their encounter rate according to size:
            ExpectedEncounters = ASSIGN_ENCOUNTER_RATES_STATES(n_track{tt}, b_track{tt}, G0, alpha, phi, Dct);
    
            % 1. All individuals grab the resources they can:
            ResourcesPerCapita = ALLOCATE_RESOURCES_SIMPLY_STATES(n_track{tt}, b_track{tt}, ExpectedEncounters);
            % ResourcesPerCapita = ALLOCATE_RESOURCES_RANDOMLY_STATES(n_track{tt}, b_track{tt}, ExpectedEncounters);
    
            % 2. Assign a growth rate according to their encounters:
            SpecificGrowthRates = SPECIFIC_GROWTH_RATES_STATES(n_track{tt}, ResourcesPerCapita, Mu, K_m);
    
            % 3. Assign a death rate:
            SpecificDeathRates = SPECIFIC_DEATH_RATES_STATES(n_track{tt}, Delta);
    
            % 4. Growth-Decay occurs:
            [n_track{tt+1}, b_track{tt+1}] = COMPETITION_ROUND_STATES(n_track{tt}, b_track{tt}, SpecificGrowthRates, SpecificDeathRates, CellDivisionThreshold, CellDeathThreshold);
    
            % Record measurables:
            NumInd = [NumInd; length(b_track{tt+1})];
            n_frac(tt+1,1) = length(find(n_track{tt+1}==2))/length(n_track{tt+1});
            b_frac(tt+1,1) = sum(b_track{tt+1}(find(n_track{tt+1}==2)))/sum(b_track{tt+1});
        end
        
    end

    % Show the results of this simulation:
    if FigViz == 1 || FigViz == 2
        figure('units','centimeters','position',[3,3,8,8]);
        hold on; box on; set(gca,'linewidth',1);
        plot(0:length(n_frac)-1, n_frac, 'k-','linewidth',2);
        plot(0:length(b_frac)-1, b_frac, 'r-','linewidth',2);
        ylim([0,1]);
        set(gca,'layer','top');
        xlabel('Rounds');
        ylabel('Fraction');
        set(gca,'fontsize',7);
        legend({'# Fraction','Biomass Fraction'});
        
        figure('units','centimeters','position',[11,3,8,8]); 
        hold on; box on; set(gca,'linewidth',1);
        plot(0:length(n_frac)-1, NumInd,'k-','linewidth',1);
        xlabel('Rounds');
        ylabel('Num Ind.');
    end


end