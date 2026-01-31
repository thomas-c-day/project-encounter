function [n_next, b_next, Delta_C] = COMPETITION_ROUND_STATES(n, b, G, D, D_T, X_T, Phage_Death_Prob, N_Phage, Phage_Burst_Size)
    % n: population list
    % b: biomass list
    % G: growth rate list
    % D: death rate list
    % D_T: threshold for division into multiple individuals
    % X_T: threshold for death below a certian size
    % N_Phage: number of phage particles attached to each individual
    % PDP: Probability of death by phage (percent chance for each encounter
    % to kill the individual).

    % Growth-decay for each individual:
    db    = b .* (G + D);
    b_now = b + db;
    n_now = n;

    % Once they are below a certain size, they die:
    X_Thresh = X_T(n);
    ix_x1 = find(b_now < X_Thresh);

    % Include interactions with phage particles: There is a chance that
    % encounters with phage kills the whole group immediately.
    PDP = 1 - (1-Phage_Death_Prob).^N_Phage;
    ix_x2 = find(rand(size(PDP)) < PDP);
    % disp([' # Phage deaths = ',num2str(length(ix_x2))]);

    % Increase the concentration of phage in the system via the burst size:
    Delta_C = sum(b_now(ix_x2) * Phage_Burst_Size);

    % Combine the death lists:
    ix_x = unique([ix_x1; ix_x2]);
    % ix_x = ix_x1;

    % Do the death:
    b_now(ix_x) = [];
    n_now(ix_x) = [];

    % Single cell states will split into new cells if they get big enough:
    D_Thresh = D_T(n_now);
    ix_d = find(b_now > D_Thresh);

    % Do the cell division:
    n_next = n_now;
    b_next = b_now;
    for nn = 1:length(ix_d)
        n_next = [n_next; n_now(ix_d(nn))];         % add new individuals to the list
        b_next(ix_d(nn)) = b_next(ix_d(nn))/2;      % split the dividing individual in half, i.e. into a single cell
        b_next = [b_next; b_next(ix_d(nn))];        % ... and append the other half
    end

end