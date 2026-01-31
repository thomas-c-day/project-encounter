function [Amount_Of_Resources] = ALLOCATE_RESOURCES_PHAGE_RANDOMLY_STATES(n, b, G_F, G_P)
    % Thomas C. Day
    % Randomly assign the amount fo food (F) and phage (P) that each
    % individual encounters during a time step. Drawn from a Poisson
    % distribution for each.
    % n: population list
    % b: amount of biomass (number of cells) for each ind. in the pop.
    % G_F: Mean encounter rate of food
    % G_P: Mean encounter rate with phage
    F = zeros(length(n), 2); % this will be the food and phage allocated to each individ.
    for nn = 1:length(n)
        F(nn,1) = poissrnd(G_F(nn), 1); % food that each unit achieves is drawn from a poisson distribution
        F(nn,2) = poissrnd(G_P(nn), 1); % phage that each unit encounters
    end
    % b_divide = [b, ones(size(b))];
    Amount_Of_Resources = F;
        % Food is contained in the first column, phage is contained in the
        % second column.

end