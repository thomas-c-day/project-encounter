function [radius_bac] = BIOMASS_TO_RADIUS(biomass_bac, packing_fraction)
    % Input either an individual or a vector of biomasses (ncells).
    % Output the list of individual radii in cm.
    
    vol_cell = 4/3 * pi * (5e-5)^3; % cubic cm
    vol_bac = vol_cell * biomass_bac ./ packing_fraction; % in cubic cm
    vol_bac(biomass_bac == 1) = vol_cell;
    radius_bac = ((3*vol_bac)/(4*pi)).^(1/3); % in cm
end