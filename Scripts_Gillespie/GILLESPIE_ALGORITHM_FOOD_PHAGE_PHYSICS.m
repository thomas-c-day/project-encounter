function [LOG] = GILLESPIE_ALGORITHM_FOOD_PHAGE_PHYSICS(INITIAL_STATE, INITIAL_BIOMASS, FOOD_CONC,  PHAGE_CONC, P_DEATH, BURST_SIZE, INITIAL_WATERVOL, TMAX, K, FIGVIZ)
    % Thomas C. Day
    % Run a Gillespie-style algorithm for particle encounters between bacteria
    % and different types of particles. The particle classes are:
    % FOOD
    % PHAGE
    % The physical processes that lead to these encounters are primarily
    % diffusion, and our empirically measured size-dependent encounter kernel.
    
    % We track the following properties of the bacteria:
    % Number
    % Size
    % Amount of phage attached
    % Amount of food attached
    
    % We track the following properties of the particles:
    % Number concentration [#/mL]
    % Water volume [mL]

    % Initializing the arrays to track:
    STATE_BAC   = INITIAL_STATE;
    BIOMASS_BAC = INITIAL_BIOMASS;
    WATER_VOL   = INITIAL_WATERVOL;
    BANK_FOOD   = zeros(size(STATE_BAC));
    BANK_PHAGE  = zeros(size(STATE_BAC));
    BANK_DEATH  = Inf*ones(size(STATE_BAC));

    % Time intervals:
    t = 0;
	next_log_time = 0;
	log_time_interval = 20*60; % save frequency [secs]

    % Initialize LOG:
    LOG = struct('BIOMASS',BIOMASS_BAC,'STATE',STATE_BAC,'BANK_FOOD',BANK_FOOD,'BANK_PHAGE',BANK_PHAGE,'BANK_DEATH',BANK_DEATH,'Time',0);

    % While loop until simulation stop condition has been reached
    while t < TMAX

        % Progress bar:
        NUM_IND = length(STATE_BAC);
        % disp(['t = ' num2str(t/60, 2) ', no. microbes = ' num2str(NUM_IND)]);

        % Measurements ----------------------------------------------------
        if t > next_log_time
            LOG(end+1).BIOMASS  = BIOMASS_BAC;
            LOG(end).STATE      = STATE_BAC;
            LOG(end).BANK_FOOD  = BANK_FOOD;
            LOG(end).BANK_PHAGE = BANK_PHAGE;
            LOG(end).BANK_DEATH = BANK_DEATH;
            LOG(end).Time       = t;
            next_log_time       = next_log_time + log_time_interval;
        end

        % Stop condition --------------------------------------------------
        ix_2 = find(STATE_BAC == 2);
        IsInfected = BANK_DEATH(ix_2);
        if all(IsInfected ~= Inf)
            % the MS state is completely infected, stop the simulation
            % Store last timepoint:
            % disp('MS state got infected!');
            LOG(end+1).BIOMASS  = BIOMASS_BAC;
            LOG(end).STATE      = STATE_BAC;
            LOG(end).BANK_FOOD  = BANK_FOOD;
            LOG(end).BANK_PHAGE = BANK_PHAGE;
            LOG(end).BANK_DEATH = BANK_DEATH;
            LOG(end).Time       = t;
            break;
        end

        % GROWTH + REACTIONS ----------------------------------------------
        % Calculate rates of particles acquisition and assign a reaction
        % Calc encounter rates with food and phage
        radius_food = 5e-5;
        radius_phage = 1e-5;
        R_food  = computeFoodAcquisitionRate(BIOMASS_BAC, radius_phage, PHAGE_CONC);
        % R_phage = computePhageAcquisitionRate(BIOMASS_BAC, PHAGE_CONC);

        % Remove members who have already died:
        ix_dead = find(BANK_DEATH ~= Inf);
        R_food(ix_dead)  = 0;
        % R_phage(ix_dead) = 0;

        % Cumulative rate:
        R_total = sum(R_food);% + sum(R_phage); % sum up all the individual rates
        R_net = R_food;

        % Draw time to the next event:
        dt = exprnd(1/R_total);
        t  = t + dt; % increment the time

        % Growth/decay all individuals: -------------------------------------------
        % [BIOMASS_BAC, BANK_FOOD] = MONOD_GROWTH(STATE_BAC, BIOMASS_BAC, BANK_FOOD, MU_MAX, K_M, DELTA, YIELD, dt);

        % Choose which type of reaction to perform: -----------------------
        x = [sum(R_food)./R_total, 1];
        a = rand;
        b = a < x;
        % reaction_ix = find(b,1,'first');
        reaction_ix = 1;

        % Do reaction: ----------------------------------------------------
        if reaction_ix == 0
            % Food acquisition
            [~, ~, BANK_FOOD] = TRACK_UPDATE_ARRAYS(BANK_FOOD, R_food, FOOD_CONC, WATER_VOL);

        else
            % Phage acquisition
            % R_net = R_phage + R_food;
            [~, IX_IND, BANK_PHAGE] = TRACK_UPDATE_ARRAYS(BANK_PHAGE, R_net, PHAGE_CONC, WATER_VOL);
            
            % Individual that acquired phage dies with some probability,
            % and increases the concentration of phage:
            if rand < P_DEATH % probability of dying to infectious phage
                % disp(['KILLING # ',num2str(IX_IND)]); pause(0.1);
                
                num_phage_add = BURST_SIZE * BIOMASS_BAC(IX_IND);
                num_phage_net = PHAGE_CONC * WATER_VOL;
                num_phage_net = num_phage_net + num_phage_add;
                PHAGE_CONC = num_phage_net/WATER_VOL;

                % Kill
                BANK_DEATH(IX_IND) = t;
                % STATE_BAC(IX_IND) = [];
                % BIOMASS_BAC(IX_IND) = 0;
                % BANK_FOOD(IX_IND) = [];
                % BANK_PHAGE(IX_IND) = [];
            end

        end

        % Do cell division ------------------------------------------------
        %{
        ix_1 = find(STATE_BAC == 1);
        for ii = 1:length(ix_1)
            if BIOMASS_BAC(ix_1(ii)) > CELL_DIV_THRESH
                new_size = BIOMASS_BAC(ix_1(ii))/2;
                new_food = 0; %BANK_FOOD(ix_1(ii))/2;

                % Updates:
                BIOMASS_BAC(ix_1(ii)) = new_size;
                BIOMASS_BAC = [BIOMASS_BAC; new_size];
                STATE_BAC = [STATE_BAC; 1];
                BANK_FOOD(ix_1(ii)) = new_food;
                BANK_FOOD = [BANK_FOOD; new_food];
                BANK_PHAGE = [BANK_PHAGE; 0];
            end
        end
        %}

        % Do cell death ---------------------------------------------------
        %{
        ix_kill = find(BIOMASS_BAC < CELL_DEATH_THRESH);
        BIOMASS_BAC(ix_kill) = [];
        STATE_BAC(ix_kill) = [];
        BANK_FOOD(ix_kill) = [];
        BANK_PHAGE(ix_kill) = [];
        %}

        % Show options:
        if FIGVIZ==1
            figure(1);
            % hold on;
            % plot(t, FOOD_CONC,'ko','linewidth',1);
            hold on;
            plot(t, PHAGE_CONC,'rx','linewidth',1);
            % yyaxis right;
            % hold on;
            % % plot(BANK_FOOD);
            % plot(BANK_PHAGE);
            % plot(t, WATER_VOL, 'bs','linewidth',1);
            % title(['Time: ',num2str(t/60),' [min]']);
            drawnow;
        end
        

    end

end