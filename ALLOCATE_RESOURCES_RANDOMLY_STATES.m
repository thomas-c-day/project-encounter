function [Amount_Of_Resources] = ALLOCATE_RESOURCES_RANDOMLY(n, b, G)
    % Amount of resource that each unit acquires is drawn from a Poisson
    % distribution with mean given by mean encounter rate.
    % n: population list
    % b: amount of biomass (number of cells) for each ind. in the pop.
    % G: Mean Encounter Rate, depends on strategy
    F = zeros(size(n)); % this will be the food allocated to each individ.
    for nn = 1:length(n)
        F(nn) = poissrnd(G(nn), 1); % food that each unit achieves is drawn from a poisson distribution
    end
    Amount_Of_Resources = F./b; % scale to measure resources per capita
end