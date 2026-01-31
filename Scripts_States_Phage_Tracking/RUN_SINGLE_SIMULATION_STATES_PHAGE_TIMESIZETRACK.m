function [n_frac, b_frac, n_track, b_track, p_track, NumInd, c_track, InfectionArray] = RUN_SINGLE_SIMULATION_STATES_PHAGE_TIMESIZETRACK(Nrounds, G0_F, G0_P, R0_F, R0_P, alpha, phi, Dct, Mu, K_m, Delta, Kappa, CellDivisionThreshold, CellDeathThreshold, First_Infection_TrackYN, InfectProb, PhageDeathProb, PhageBurstSize, FigViz)
    % Run one simulation of a competition between two different states
    % within one genotype over a number of selection rounds. The two
    % different states will represent a single-celled state and a
    % multi-celled state. Individuals can encounter both resources and
    % phage particles.

    % Initial conditions:
    n0 = [ones(5e2,1); 2*ones(5e2,1)];
        % State 1 is the single-celled strategy, State 2 is the
        % multi-celled strategy. Half the individuals are of State 2.

    % Notate the biomass of each individual:
    b0 = 0.10*randn(size(n0)) + 1;
        % start off as each individual is single-celled, each cell is
        % slightly different in size

    % Variables to track:
    n_track{1,1} = n0; % track the individuals and their state
    b_track{1,1} = b0; % track the biomass of each individual
    p_track{1,1} = zeros(size(n0)); % track how many times each individual has been infected
    n_frac(1,1)  = length(find(n0==2))./length(n0); % number fraction of individuals in MC state
    b_frac(1,1)  = sum(b0(find(n0==2)))/sum(b0); % biomass fraction in MC state
    InfectionArray = zeros(1e3,2);
    NumInd = length(b0);

    % Run competition experiment:
    for tt = 1:Nrounds
        
        % Progress bar:
        %disp([num2str(tt), ' / ',num2str(Nrounds)]);

        % Checks to stop the simulation:
        if length(n_track{tt}) > Kappa
        % if sum(b_track{tt}) > Kappa % change this to be a carrying-capacity for biomass
            % Carrying capacity is reached
            break;
        elseif length(n_track{tt}) <= 1
            % Extinction
            break;
        elseif n_frac(tt) > 1%0.99
            % Fixation
            break;
        elseif n_frac(tt) == 0
            % Other fixation
            break;
        elseif ~any(p_track{tt}==0)
            % Everyone has been infected at least once
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
            ExpectedEncounters_Food  = ASSIGN_ENCOUNTER_RATES_STATES(n_track{tt}, b_track{tt}, G0_F, R0_F, alpha, phi, Dct);
            ExpectedEncounters_Phage = ASSIGN_ENCOUNTER_RATES_STATES(n_track{tt}, b_track{tt}, G0_P, R0_P, alpha, phi, Dct);
    
            % 1. All individuals grab the resources they can:
            % ResourcesPerCapita = ALLOCATE_RESOURCES_SIMPLY_STATES(n_track{tt}, b_track{tt}, ExpectedEncounters);
            Resources = ALLOCATE_RESOURCES_PHAGE_RANDOMLY_STATES(n_track{tt}, b_track{tt}, ExpectedEncounters_Food, ExpectedEncounters_Phage);
            % Resources

            % 2. Assign a growth rate according to their per-capita food encounters:
            ResourcesPerCapita = Resources(:,1)./b_track{tt};
            SpecificGrowthRates = SPECIFIC_GROWTH_RATES_STATES(n_track{tt}, ResourcesPerCapita, Mu, K_m);
                % Only depends on encounters with food.
    
            % 3. Assign a death rate:
            SpecificDeathRates = SPECIFIC_DEATH_RATES_STATES(n_track{tt}, Delta);
                % There is a baseline level of decay for each state. This
                % does not include encounters with phage.

            % 3.5 Get the number of phage infections occuring in this batch
            % and add it to the tracked number of times each group has
            % been infected:
            q = 1 - (1-InfectProb).^Resources(:,2);
            infection_TF = rand(size(q)) < q;
            p_track{tt+1} = zeros(size(infection_TF));
            %p_track{tt+1} = p_track{tt} + infection_TF; % accumulates the number of infections
    
            % 4. Growth-Decay occurs:
            [n_track{tt+1}, b_track{tt+1}, DeltaC_Phage] = COMPETITION_ROUND_STATES(n_track{tt}, b_track{tt}, ...
                SpecificGrowthRates, SpecificDeathRates, CellDivisionThreshold, CellDeathThreshold, ...
                PhageDeathProb, infection_TF, PhageBurstSize);
    
            % Increase the rate of encounters with phage according to the
            % changes in DeltaC:
            DeltaG = DeltaC_Phage*1e-7;
            G0_P = G0_P + DeltaG;
            c_track(tt) = G0_P/1e-7;

            % Record measurables:
            NumInd = [NumInd; length(b_track{tt+1})];
            n_frac(tt+1,1) = length(find(n_track{tt+1}==2))/length(n_track{tt+1});
            b_frac(tt+1,1) = sum(b_track{tt+1}(find(n_track{tt+1}==2)))/sum(b_track{tt+1});

            % Track how long until each individual of type 2 gets first infected:
            if First_Infection_TrackYN
                % Check if anyone has been infected:
                if any(infection_TF)
                    
                    % Who was infected? 
                    ix = find(infection_TF);
                    type_inf = n_track{tt}(ix);
                    ix = ix(type_inf == 2); % only care about type 2

                    % Determine if they have already been infected before
                    for ii = 1:length(ix)
                        if InfectionArray(ix(ii),1) > 0
                            % already been infected
                            continue;
                        else
                            % not previously infected
                            InfectionArray(ix(ii),1) = tt;
                            InfectionArray(ix(ii),2) = b_track{tt}(ix(ii));
                        end
                    end


                %{
                    if length(ix) > 1
                        rs = randsample(length(ix),1); % randomly sample one of the instances
                        type_inf = n_track{tt}(ix(rs));
                        size_inf = b_track{tt}(ix(rs));
                        time_inf = tt;
                    else
                        type_inf = n_track{tt}(ix);
                        size_inf = b_track{tt}(ix);
                        time_inf = tt;
                    end
                    TST_inf = [type_inf, size_inf, time_inf];
                    break;
                %}
                else
                    % Nobody has been infected, continue the simulation
                    continue;
                end
            end
            
        end
        
    end

end