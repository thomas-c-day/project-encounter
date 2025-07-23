function [Specific_Growth_Rates] = SPECIFIC_GROWTH_RATES(n, F, mu, K)
    % According to resources, assign a specific growth rate to each member
    % of the population. Monod dynamics, dependent on resources per capita.

    Specific_Growth_Rates = mu(n) .* F ./(F + K(n));

end