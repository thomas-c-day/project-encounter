% Thomas C. Day
% Validation and analysis of segmentations from confocal data.

% Load ellipsoid fits, make measurements of average cell size:
seg_folder = 'E:\My Drive\01_Data\Project_Encounter\PackingMeasurements\02_Segs\';
save_folder = 'E:\My Drive\01_Data\Project_Encounter\PackingMeasurements\03_Plots_Results\';
filelist = dir([seg_folder,'*.mat']);
FIGVIZ = 1;

% Record resolutions per voxel for each image:
XYres = repmat(0.135, 10, 1);
XYres(6) = 0.144;

% Put on a size filter:
MinSize = 0.05; % cubic microns

% Preallocate measurements:
MeanCellVol = zeros(length(filelist),1);
ErrCellVol  = MeanCellVol;
MeanCellAsp = zeros(length(filelist),1);
ErrCellAsp  = MeanCellAsp;
GroupV      = MeanCellVol;
GroupVEst   = MeanCellVol;
GroupN      = MeanCellVol;
GroupR      = MeanCellVol;
GroupPhi    = MeanCellVol;
GroupAS     = cell(length(filelist),1);

% Loop through each aggregate:
for ff = 1%:length(filelist)
    disp(['Aggregate # ',num2str(ff),' / ',num2str(length(filelist))]);

    load([seg_folder,filelist(ff).name]);

    % Record Cell Volumes and Aspect Ratios fit from Ellipsoids:
    CellVol = zeros(length(Part4.Ellipsoids),1);
    CellAsp = CellVol;
    for ii = 1:length(Part4.Ellipsoids)
        CellVol(ii) = XYres(ff)^3 * 4/3 * pi * Part4.Ellipsoids(ii).S(1,1) * Part4.Ellipsoids(ii).S(2,2) * Part4.Ellipsoids(ii).S(3,3);
        CellAsp(ii) = Part4.Ellipsoids(ii).S(1,1) / Part4.Ellipsoids(ii).S(3,3);
    end

    % Record cell volumes from a convex hull of the surface:
    V = zeros(length(Part4.Surfaces),1);
    for kk = 1:length(Part4.Surfaces)
        [~,V(kk,1)] = convhull(XYres(ff) * Part4.Surfaces(kk).vertices);
    end
    ix = find(V < MinSize);

    % Weird fits filter:
    ix2 = find(CellVol < 0);
    ix3 = find(CellAsp < 0);
    ix4 = find(CellAsp > 10);
    ix_net = unique([ix; ix2; ix3; ix4]);

    % Generate colors based on weird ellipsoid fits:
    Colors = repmat(.5, length(Part4.Ellipsoids), 3);
    for ii = 1:length(ix_net)
        Colors(ix_net(ii),:) = [1,0,0];
    end
    
    % Filter out the weird fits:
    FilteredSurfaces = Part4.Surfaces;
    FilteredEllipsoids = Part4.Ellipsoids;
    FilteredSurfaces(ix_net) = [];
    FilteredEllipsoids(ix_net) = [];
    CellVol(ix_net) = [];
    CellAsp(ix_net) = [];    

    % Make final measurements:
    CellCenters = zeros(length(FilteredEllipsoids),3);
    CellLengths = zeros(length(FilteredEllipsoids),3);
    for kk = 1:length(FilteredEllipsoids)
        CellCenters(kk,:) = XYres(ff) * [FilteredEllipsoids(kk).T(1,end), FilteredEllipsoids(kk).T(2,end), FilteredEllipsoids(kk).T(3,end)];
        CellLengths(kk,:) = XYres(ff) * [FilteredEllipsoids(kk).S(1,1), FilteredEllipsoids(kk).S(2,2), FilteredEllipsoids(kk).S(3,3)];
    end
    GroupCOM        = mean(CellCenters);
    CellCenters     = CellCenters - GroupCOM;
    rCOM            = (CellCenters)./vecnorm(CellCenters, 2, 2);
    CellCenters_Sh  = CellCenters + mean(CellLengths(:,1)).*rCOM;
    CellCenters_BS  = CellCenters_Sh + GroupCOM;

    % Estimate a group volume from z-projection cross-section:
    CellCenters_Zproj = CellCenters_Sh;
    CellCenters_Zproj(:,3) = 0;
    [~,CrossSectionalArea] = convhull(CellCenters_Zproj(:,1), CellCenters_Zproj(:,2));
    estimated_radius = sqrt(CrossSectionalArea/pi);

    % Record net measurements:
    MeanCellVol(ff) = mean(CellVol); ErrCellVol(ff) = std(CellVol);
    MeanCellAsp(ff) = mean(CellAsp); ErrCellAsp(ff) = std(CellAsp);
    GroupN(ff)      = length(FilteredEllipsoids);
    GroupAS{ff}     = alphaShape(CellCenters_BS, 5);
    GroupV(ff)      = volume(GroupAS{ff});
    GroupVEst(ff)   = 4/3*pi*estimated_radius.^3;
    GroupR(ff)      = ( (3 * GroupV(ff)) / (4*pi) ).^(1/3);
    GroupPhi(ff)    = GroupN(ff) * MeanCellVol(ff) / GroupV(ff);

    % Optional, show the good and weird fits together:
    if FIGVIZ
        figure; 
        ax = axes('position',[0,.1,1,.9]);
        hold on; box on; set(gca,'linewidth',1);
        for kk = 1:length(Part4.Surfaces)
            patch('faces',Part4.Surfaces(kk).faces,'vertices', XYres(ff) * Part4.Surfaces(kk).vertices, 'facecolor',Colors(kk,:), 'edgecolor','none','facealpha',0.7);
            % text(Part4.Ellipsoids(kk).T(1,end), Part4.Ellipsoids(kk).T(2,end), Part4.Ellipsoids(kk).T(3,end), num2str(kk));
        end
        plot(GroupAS{ff}, 'facecolor',[.5,.5,.5],'edgecolor','k','facealpha',0.2);
        view(3); axis equal;
        camlight; material dull;
        % print([save_folder,filelist(ff).name(1:end-4)],'-dpng','-r500');
        % close(gcf);
    end
end

%% Make plots:

% plotfolder = 'E:\My Drive\01_Data\Project_Encounter\PackingMeasurements\03_Plots_Results\';
% load([plotfolder,'DATA_Figure_04_Inset.mat']);

% Cell volume:
figure; hold on; box on; set(gca,'linewidth',1);
errorbar(1:length(GroupN), MeanCellVol, ErrCellVol, 'k.-','linewidth',1,'markersize',15);
xlabel('Aggregate #');
ylabel('Cell Volume [um^3]');
% print([plotfolder,'CellVolumes'],'-dpng','-r500');

% Cell Aspect Ratio:
figure; hold on; box on; set(gca,'linewidth',1);
errorbar(1:length(GroupN), MeanCellAsp, ErrCellAsp, 'k.-','linewidth',1,'markersize',15);
xlabel('Aggregate #');
ylabel('Cell Asp');
% print([plotfolder,'CellAsps'],'-dpng','-r500');

% Fit a regression to Group volume and number of cells:
c_lin = polyfit(GroupV, GroupN,1);
disp(['Regression equation is y = ' num2str(c_lin(1)) '*x + ' num2str(c_lin(2))]);
y_fit = polyval(c_lin, GroupV);

% Fit a regression to estimated group volume and number of cells:
c_est = polyfit(GroupVEst, GroupN, 1);
disp(['Regression equation is y = ' num2str(c_est(1)) '*x + ' num2str(c_est(2))]);
y_est = polyval(c_est, GroupVEst);
residuals = GroupN - y_est;
residuals_std = std(residuals);
n = length(GroupN);
x_mean = mean(GroupVEst);
x_var = var(GroupVEst);
m_std_error = residuals_std/sqrt(x_var*(n-1))

% Fit a regression to the log-transformed to regress the power of the
% relationship:
LogV = log10(GroupV);
LogN = log10(GroupN);
c_log = polyfit(LogV, LogN, 1);
yfitlog = polyval(c_log, LogV);
residuals = LogN - yfitlog;
residuals_std = std(residuals);
n = length(LogV);
x_mean = mean(LogV);
x_var = var(LogV);
m_std_error = residuals_std/sqrt(x_var*(n-1));
disp(['Power Law Regression is y = ' num2str(c_log(1)) '*x + ' num2str(c_log(2))])
disp(['Slope fit: m = ',num2str(c_log(1)),' +- ',num2str(m_std_error), ' , 1 sigma']);

% Group Ncells vs Volume:
figure; 
hold on; box on; set(gca,'linewidth',1);
plot(GroupV, GroupN, 'k.','markersize',15);
plot(GroupV, y_fit, 'k-','linewidth',1);
xlabel('Group Volume [um^3]');
ylabel('Ncells');
% print([plotfolder,'Ncells_vs_Vols'],'-dpng','-r500');

% Packing Fraction:
figure; hold on; box on; set(gca,'linewidth',1);
plot(GroupR, GroupPhi, 'k.','markersize',15);
for nn = 1:length(GroupR)
    text(GroupR(nn), GroupPhi(nn), num2str(nn));
end
xlabel('Radius [um]');
ylabel('\phi');
% print([plotfolder,'PackingFractions'],'-dpng','-r500');

% Actual volume from confocal vs. estimated volume from cross-section:
figure;
hold on; box on; set(gca,'linewidth',1);
plot(GroupV, GroupVEst, 'k.','markersize',12);
plot(GroupV, GroupV, 'k-','linewidth',1);
xlabel('Confocal volume');
ylabel('Estimated volume');

% Plot estimate volume vs ncells:
figure; 
hold on; box on; set(gca,'linewidth',1);
plot(GroupVEst, GroupN, 'k.','markersize',15);
plot(GroupVEst, y_est, 'k-','linewidth',1);
xlabel('Group Volume [um^3]');
ylabel('Ncells');
% print([plotfolder,'Ncells_vs_Vols'],'-dpng','-r500');

% save([plotfolder,'DATA_Figure_04_Inset.mat'],'c_lin','c_log','c_est','m_std_error','GroupAS','GroupV','GroupVEst','GroupN','GroupR','GroupPhi','MeanCellVol','MeanCellAsp','ErrCellVol','ErrCellAsp');
% close all;