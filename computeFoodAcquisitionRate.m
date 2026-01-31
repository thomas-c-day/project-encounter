function [Encounter_Rate] = computeFoodAcquisitionRate(biomass_bac, radius_food, food_conc)

    % Parameters (these will eventually be called in to the function but
    % for now they are hard-coded):
    packing_fraction = 0.4;
    lambda = 2.8;
    A = 2e-2; % typical values from bead experiments over 30 mins.
    A = A / (30*60); % convert to rate per second;
    A = A ./ 1e7; % convert by number of beads per mL concentration
    A = A * 1e12; % cubic microns per mL
    % radius_food = 5e-5; % cm

    % Convert biomass (ncells) to length scale:
    radius_bac = BIOMASS_TO_RADIUS(biomass_bac, packing_fraction);

    % Calculate encounter kernel for each individual:
    % radius_bac = 5e-4*ones(size(radius_bac)); % cm
    Encounter_Kernel = A * (radius_bac + radius_food).^lambda; % in mL/sec
    Encounter_Kernel(radius_bac == 0) = 0;
    Encounter_Rate = food_conc .* Encounter_Kernel; % in #/sec

end