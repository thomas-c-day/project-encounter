% Thomas C. Day
% This simulation will model encounters between (non-growing) bacteria and
% phage. The phage are characterized by their prevalence (concentration)
% and burst size. The bacteria can either be dispersed as single cells or
% clumped into an aggregate. They encounter phage accordingly.
%
% When bacteria encounter phage, they immediately die and release more
% phage, according to the phage burst size, thus increasing the net
% concentration of phage.
%
% Groups can be modeled as either entirely dying, or having some kind of
% immune response that allows some of the cells to remain living, still
% clumped into groups.
%
% Single cells will be modeled as either completely inert, or as having
% motility that is uncorrelated to the phage, i.e. there is no chemotaxis.
%
% Which strategy will allow bacterial cells to persist longer?
%
% -------------------------------------------------------------------------
% Inputs:
T_MAX           = 1e3*24*60*60; % seconds
UM3_TO_M3       = 1e-18;
NCells          = 1e3;
C_PFU           = 1e5; % PFUs/mL
BurstSize       = 200; % typical burst size value of phage
PackingFraction = 0.5;
VolumeCell      = UM3_TO_M3 * 4/3*pi*0.5^3; % cubic microns
RadiusCell      = 0.5e-6; % [meters]
VolumeGroup     = NCells * VolumeCell/PackingFraction; % cubic meters
RadiusGroup     = (3*VolumeGroup/(4*pi)).^(1/3); % meters
RadiusPhage     = 1e-7; % meters
Dphage          = 1e-11; % [m^2/s]
Dbacteria       = 8e-10; % [m^2/s]

% -------------------------------------------------------------------------
% Run the experiment for single cells


% Set up parameters:
ind_size_list = RadiusCell*ones(NCells, 1); % in meters
t = 0;
counter = 1;
SimTracker = zeros(length(ind_size_list)+1, 3);
EncounterModelParams.C_PFU    = C_PFU;
EncounterModelParams.BASERATE = 1;
EncounterModelParams.LAMBDA   = 3;
EncounterModelParams.D_PHAGE  = Dphage;
EncounterModelParams.D_BACTERIA = Dbacteria;
EncounterModelParams.R_PHAGE  = RadiusPhage;

while ~isempty(ind_size_list)
    % Progress bar:
    disp(['t = ' num2str(t, 2) ', no. inds = ' num2str(length(ind_size_list))]);

    % Calculate all encounter rates:
    [K,P] = ENCOUNTER_MODEL(ind_size_list, EncounterModelParams);
    
    % Draw time to next event:
    R = sum(K);
    p = K./R;
    dt = exprnd(1/R);
    t = t + dt;
    counter = counter + 1;

    % Draw which individual is getting infected:
    idx = randsample(length(ind_size_list), 1, true, p);

    % Do the infection:
    total_burst = BurstSize;
    ind_size_list(idx) = []; % individual is killed
    EncounterModelParams.C_PFU = EncounterModelParams.C_PFU + total_burst; % add phage to the pool
    
    % Track the simulation:
    SimTracker(counter,:) = [t, length(ind_size_list), EncounterModelParams.C_PFU];
end
SimTracker(1,:) = [0, NCells, C_PFU];

% Show sim results:
figure;
tiledlayout('flow');
ax1 = nexttile;
hold on; box on; set(gca,'linewidth',1);
plot(SimTracker(:,1), SimTracker(:,2), '-','markersize',12,'linewidth',1,'color',[.2,.4,.8]);
xlabel('Time [s]');
ylabel('Num. Ind');
ax2 = nexttile;
hold on; box on; set(gca,'linewidth',1);
plot(SimTracker(:,1), SimTracker(:,3), '-','markersize',12,'linewidth',1,'color',[.8,.2,.2]);
xlabel('Time [s]');
ylabel('PFUs/mL');

% Run simulation for groups of cells:
% Set up parameters:
ind_size_list = RadiusGroup; % in meters
t = 0;
counter = 1;
SimTrackerG = zeros(length(ind_size_list)+1, 3);
EncounterModelParams.C_PFU    = C_PFU;
EncounterModelParams.BASERATE = 1;
EncounterModelParams.LAMBDA   = 3;
EncounterModelParams.D_PHAGE  = Dphage;
EncounterModelParams.D_BACTERIA = 0;
EncounterModelParams.R_PHAGE  = RadiusPhage;

while ~isempty(ind_size_list)
    % Progress bar:
    disp(['t = ' num2str(t, 2) ', no. inds = ' num2str(length(ind_size_list))]);

    % Calculate all encounter rates:
    [K,P] = ENCOUNTER_MODEL(ind_size_list, EncounterModelParams);
    
    % Draw time to next event:
    R = sum(K);
    p = K./R;
    dt = exprnd(1/R);
    t = t + dt;
    counter = counter + 1;

    % Draw which individual is getting infected:
    idx = randsample(length(ind_size_list), 1, true, p);

    % Do the infection:
    total_burst = BurstSize;
    ind_size_list(idx) = []; % individual is killed
    EncounterModelParams.C_PFU = EncounterModelParams.C_PFU + total_burst; % add phage to the pool
    
    % Track the simulation:
    SimTrackerG(counter,:) = [t, length(ind_size_list), EncounterModelParams.C_PFU];
end
SimTrackerG(1,:) = [0, NCells, C_PFU];

axes(ax1);
hold on;
line([SimTrackerG(2,1), SimTrackerG(2,1)], [0,NCells], 'linestyle','--','color',[.2,.4,.8],'linewidth',1);

% -------------------------------------------------------------------------
% Functions

function [K,P] = ENCOUNTER_MODEL(IND_SIZE_LIST, MODEL_PARAMS)
    % Generate an encounter kernel:
    K = MODEL_PARAMS.C_PFU * (MODEL_PARAMS.BASERATE * (IND_SIZE_LIST + MODEL_PARAMS.R_PHAGE).^MODEL_PARAMS.LAMBDA + ...
        4*pi*(MODEL_PARAMS.D_PHAGE + MODEL_PARAMS.D_BACTERIA) * (IND_SIZE_LIST + MODEL_PARAMS.R_PHAGE) );
    P = K > 0;
end

% -------------------------------------------------------------------------