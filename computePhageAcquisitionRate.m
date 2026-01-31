function [Encounter_Rate] = computePhageAcquisitionRate(biomass_bac, phage_conc)

    % Parameters:
    D_phage = 2e-8; % cm^2/sec, using Stokes-Einstein
    % D_bac = 4e-9;
    packing_fraction = 0.5; % packing fraction of cells in groups
    radius_phage = 1e-5; % cm

    % Convert biomass (ncells) to length scale:
    radius_bac = BIOMASS_TO_RADIUS(biomass_bac, packing_fraction); % in cm

    % Calculate the diffusion constant for bacteria at a given size:
    r = radius_bac * 1e-2; % in m
    kT = 4e-21; % Joules
    nu = 1e-3; % [Pa*s]
    D_bac = kT./(6*pi*nu*r); % in m^2/sec
    D_bac = 1e4 * D_bac; % convert to cm^2/sec
    % D_phage = 1e4 * kT./(6*pi*nu*(radius_phage*1e-2));

    % Calculate encounter kernel via diffusive process [cm^3/sec]:
    Encounter_Kernel = 4 * pi * (D_phage + D_bac) .* (radius_phage + radius_bac);
    Encounter_Kernel(radius_bac == 0) = 0;

    % Convert to rate [#/sec]:
    Encounter_Rate = phage_conc .* Encounter_Kernel;

end