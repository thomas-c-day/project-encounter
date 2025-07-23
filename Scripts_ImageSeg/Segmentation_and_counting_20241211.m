% Thomas C. Day, 2024
% Segmenting transmission images, then gfp images, to count the number of
% particles on each aggregate.

% Routine built specifically for the EDF images.
% INPUTS ------------------------------------------------------------------
% Inputs for TD processing:
SEGTHRESH_AG    = 0.90;                             % confidence threshold
AREA_MIN        = 10;                               % um^2
DILATION_SIZE   = 3;                                % pixels
XYRES           = 0.44;                             % um per pixel
AREA_MIN        = floor(AREA_MIN / (XYRES^2));      % convert area to number of px

% Inputs for GFP processing:
% SEGTHRESH_GFP   = 1e3;                             % GFP threshold, early images
SEGTHRESH_GFP   = 5e3;                               % later threshold
EXP_BEAD_DIAM   = 1;                                 % Expected bead diameter in microns
EXP_BEAD_DIAM   = floor(EXP_BEAD_DIAM ./ XYRES);     % Diameter, in pixels, expected of a bead

% Sub-sampling:
SubSample       = 1;                    % T/F for subsampling the set, useful for speed and visualization [0=No, 1=yes]
d               = 3e3:7e3;              % Range of the subsample in pixels
d_x             = d;
d_y             = d;

% Showing/saving:
FIGVIZ     = 4;                    % Show the figures [0=No, ...
                                                     % 1=Show results of DIC segmentation, ...
                                                     % 2=Show results of GFP segmentation, ...
                                                     % 3=Show Fluorescence per unit area plot,...
                                                     % 4=Show overlaid images and segmentations after all processing]
SAVEDATA   = 0;                    % T/F save data [0=No, 1=Yes]

% Folders where images are saved:
TDfolder   = 'F:\00_USC\03_Data_ProjectEncounter\20241211_Micobeads\TD_c\';
SPfolder   = 'F:\00_USC\03_Data_ProjectEncounter\20241211_Micobeads\TD_c\Probs\';
GFPfolder  = 'F:\00_USC\03_Data_ProjectEncounter\20241211_Micobeads\GFP\';
% savefolder = 'D:\00_USC\03_Data_ProjectEncounter\20241211_Micobeads\20241230\';
% TDfolder   = 'F:\00_USC\03_Data_ProjectEncounter\20250402_Microbeads\TD\';
% SPfolder   = 'F:\00_USC\03_Data_ProjectEncounter\20250402_Microbeads\Probs\';
% GFPfolder  = 'F:\00_USC\03_Data_ProjectEncounter\20250402_Microbeads\GFP\';
% savefolder = 'F:\00_USC\03_Data_ProjectEncounter\20250402_Microbeads\';

% Datafiles:
TDfilelist  = dir([TDfolder,'*.tif']);              % list of the brightfield images
SPfilelist  = dir([SPfolder,'*Probabilities.tif']); % list of the probability masks output from Ilastik
GFPfilelist = dir([GFPfolder,'*.tif']);             % list of the gfp images
% -------------------------------------------------------------------------

% MAIN --------------------------------------------------------------------
% Loop through images:
for ff = 1:length(GFPfilelist)
    
    % Progress bar:
    disp([num2str(ff),' / ',num2str(length(GFPfilelist))]);
    disp('Initializing...');
    disp('Reading files...');
    disp(SPfilelist(ff).name);
    disp(TDfilelist(ff).name);

    % Read images:
    IMG_TD   = imread([TDfolder, TDfilelist(ff).name]);
    IMG_SP   = tiffreadVolume([SPfolder,SPfilelist(ff).name]);
    IMG_GFP  = imread([GFPfolder,GFPfilelist(ff).name]);

    % Subsampling for testing, useful to boost speed:
    if SubSample
        IMG_TD  = IMG_TD(d_y,d_x);
        IMG_SP  = IMG_SP(d_y,d_x,:);
        IMG_GFP = IMG_GFP(d_y,d_x);
    end

    % TD Processing: ------------------------------------------------------
    [BW_TD, BOUNDS] = DIC_PROCESS(IMG_TD, IMG_SP, SEGTHRESH_AG, AREA_MIN, DILATION_SIZE, FIGVIZ);

    % GFP processing: -----------------------------------------------------
    [Xpk, Ipk, wiGFP] = GFP_PROCESS(IMG_GFP, SEGTHRESH_GFP, EXP_BEAD_DIAM, FIGVIZ);

    % Count particles inside CCs: -----------------------------------------
    [CC, AggregateArea, NumParticles, IntegratedFL] = COUNT_PARTICLES_INSIDE_CCs(BW_TD, Ipk, wiGFP);
    FL_per_Area = IntegratedFL./AggregateArea;

    if FIGVIZ == 3
        % Special conditions depending on bead size:
        ExpBeadArea = pi*(EXP_BEAD_DIAM/2 + 6)^2;
        
        figure; hold on;
        plot(AggregateArea, IntegratedFL, 'k.');
        line([ExpBeadArea, ExpBeadArea], [0,max(IntegratedFL)],'linestyle','--','color','k','linewidth',1);
        xlabel('Area');
        ylabel('Integrated FL');
        set(gca,'xscale','log');
    
        figure; hold on;
        plot(FL_per_Area,'k.');
        xlabel('Individual');
        ylabel('FL/Area');
    end
    
    % Show results of segmentation: ---------------------------------------
    GFP_LO = min(IMG_GFP,[],'all');
    GFP_HI = max(IMG_GFP,[],'all');
    if FIGVIZ == 4
        figure;
        imshowpair(imadjust(IMG_TD,[0.1,0.4]), imadjust(IMG_GFP,[0,0.2]),'blend');
        hold on; box on; set(gca,'linewidth',1);
        for bb = 1:length(BOUNDS)
            bf = sgolayfilt(BOUNDS{bb},3,11);
            plot(bf(:,2), bf(:,1),'r-','linewidth',1);
        end
        plot(Xpk(1,:), Xpk(2,:), 'b.','markersize',8);
        RP = regionprops(CC,'Circularity','Centroid');
        for ii = 1:length(FL_per_Area)
            text(RP(ii).Centroid(1), RP(ii).Centroid(2), num2str(NumParticles(ii)),'fontsize',10);
            % text(RP(ii).Centroid(1), RP(ii).Centroid(2), num2str(FL_per_Area(ii)),'fontsize',10);
        end
        drawnow;
    end

    % Save to file: -------------------------------------------------------
    if SAVEDATA
        OUTPUT.AggregateArea = AggregateArea;
        OUTPUT.NumParticles  = NumParticles;
        OUTPUT.IntegratedFL  = IntegratedFL;
        OUTPUT.FL_per_Area   = FL_per_Area;
        OUTPUT.CC            = CC;
        OUTPUT.Bounds        = BOUNDS;
        OUTPUT.Ipk           = Ipk;
        OUTPUT.XYres         = XYRES;
        OUTPUT.InputVals     = [SEGTHRESH_AG, SEGTHRESH_GFP, AREA_MIN, DILATION_SIZE, EXP_BEAD_DIAM(ff)];
        save([savefolder, SPfilelist(ff).name(1:end-4),'.mat'],'-struct','OUTPUT');
    end

    % Close all images:
    % close all;

end