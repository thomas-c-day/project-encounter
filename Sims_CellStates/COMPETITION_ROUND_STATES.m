function [n_next, b_next] = COMPETITION_ROUND_STATES(n, b, G, D, D_T, X_T)
    % n: population list
    % b: biomass list
    % G: growth rate list
    % D: death rate list
    % K: carrying capacity

    % Growth-decay for each individual:
    db    = b .* (G + D);
    b_now = b + db;
    n_now = n;

    % Once they are below a certain size, they die:
    X_Thresh = X_T(n);
    ix_x = find(b_now < X_Thresh);
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

    % Multi-celled states will revert to single-celled states under some
    % conditions:

    % Single-celled states will switch to multi-celled states under some
    % conditions:

end
