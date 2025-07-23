function [Amount_Of_Resources] = ALLOCATE_RESOURCES_SIMPLY_STATES(n, b, G)
    % Amount of resource is deterministic, it is just exactly the rate at
    % which the strategy encounters the resource.
    F = zeros(size(n));
    for nn = 1:length(n)
        F(nn) = G(nn);
    end
    Amount_Of_Resources = F./b;
end