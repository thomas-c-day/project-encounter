% Thomas C. Day
% Figure generation script for Project Encounter.
% February version after re-doing some experiments.

PRINT_FIGURES = 1;
drive_folder = 'G';
printfolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\'];
res = '-r400';
big_width = 17.8;
lil_width = 8.7;
med_width = 11.4;

%% Figure 1:
% Bacteria make aggregates of different sizes. Size matters for encounters.

clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

% Generating an image of multicellular aggregates:
img_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\Imgs\'];
IMG_HOW = imread([img_folder,'Culture_Narrow.png']);
IMG_MONTAGE = imread([img_folder, 'Montage_12B01_1e4_10umscalebar.png']);
IMG_MONTAGE = IMG_MONTAGE(:,367:end);
IMG_SPAN = imread([img_folder, 'SpanSymbol.png']);

% Getting a size histogram:
size_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\08_SizeDistribution\'];
S24 = load([size_folder,'Example_1e7_100rpm_24hrs.mat']);
S48 = load([size_folder,'Example_1e7_100rpm_48hrs.mat']);
nbins = calcnbins(S24.EqDiameter,'all');
[Counts24, Edges24] = histcounts(S24.EqDiameter, nbins.scott,'normalization','pdf');
Bins24 = mean([Edges24(1:end-1); Edges24(2:end)]);
nbins = calcnbins(S48.EqDiameter,'all');
[Counts48, Edges48] = histcounts(S48.EqDiameter, nbins.scott,'normalization','pdf');
Bins48 = mean([Edges48(1:end-1); Edges48(2:end)]);

% Append data for boxplots:
EQD = [S24.EqDiameter'; S48.EqDiameter'];
GRP = [zeros(length(S24.EqDiameter),1); ones(length(S48.EqDiameter),1)];

% 3D Encounter kernel (theory): -------------------------------------------
kT    = 1.38e-23 * 300;          % thermal energy in [J]
nu    = 1e-3;                    % water dynamic viscosity in [Pa*s]
mu    = 1e-6;                    % water kinematic viscosity (nu/rho)
DRhob = 25;                      % difference in density between particles and water [kg/m^3]
DRhor = 25;                      % difference in density between patches and water [kg/m^3]
g     = 10;                      % acceleration due to gravity [m/s^2]
eps   = 1e-5;                    % ocean energy dissipation rate, [W/kg]
rb    = logspace(-10,-3,1e3);    % bacteria agg radius in [m]
rr    = logspace(-10,-3,1e3);    % resource patch radius in [m]
Conv2UM = (1e6)^3;               % Conversion factor from cubic meters to cubic microns
Conv2CM = (1e2)^3;               % Conversion factor from cubic meters to cubic centimeters (mL)

% Make patch size matrix:
[Rb, Rr] = meshgrid(rb, rr);

% Diffusive Kernel:
Db = kT./(6*pi*nu*Rb);
Dr = kT./(6*pi*nu*Rr);
G_D = 4*pi*(Db + Dr) .* (Rb + Rr);

% Buoyancy Kernel:
Ub = 2 * DRhob * g * Rb.^2 / (9*nu);
Ur = 2 * DRhor * g * Rr.^2 / (9*nu);
G_B = pi * (Rb + Rr).^2 .* abs(Ub-Ur);

% Turbulence Kernel, chosen for a particular epsilon value, and for well below the Kolmogorov limit:
G_T1 = 1.3 * (Rb + Rr).^3 .* sqrt(eps./mu);
G_T2 = 1.37*pi.*eps ^ (1/3) * (Rb+Rr).^(7/3);
Lk = 0.5*(mu^3./eps).^(1/4);
ix1 = Rb+Rr <= Lk;
ix2 = Rb+Rr > Lk;
G_T = G_T1.*ix1 + G_T2.*ix2;

% Make a swimming kernel, too: only one item swims
% G_S = 4*pi*SwimStrength .* (Rb + Rr);

% Net kernel:
% G_N = Conv2UM * (G_D + G_B + G_T);
G_N = Conv2CM * (G_D + G_B + G_T);

% Estimate per-capita:
phi    = 0.4;              
rc     = 5e-7;
Ncells = phi*(Rb.^3)./(rc.^3);
G_Npc  = G_N./Ncells;

% Back-of-the-envelope of how many collisions between individual cells:
c_cells = 1e6; % cells per mL
BOE = 1/2 * (G_N(572,572)) * c_cells * c_cells;

NetKernel = cat(3, G_D, G_B, G_T);
[~,PhaseMap] = max(NetKernel,[],3);

% Find a regime where the diffusion kernel dominates:
IX = G_D > (G_T + G_B);
B_DIFF = bwboundaries(IX);
B_DIFF = B_DIFF{1};
X_DIFF = 1e6 * rr(B_DIFF(:,1));
Y_DIFF = 1e6 * rb(B_DIFF(:,2));

% % Find diffusion-dominated regimes for each of 0.01, 1, & 100 um diameter
% % resources:
% ix001 = find(X_DIFF < 0.005);
% ix010 = find(X_DIFF < 0.5);
% ix100 = find(X_DIFF < 50);

% 1D Slices: --------------------------------------------------------------
% Inputs:
kT   = 1.38e-23 * 300;              % thermal energy in [J]
nu   = 1e-3;                        % water dynamic viscosity in [Pa*s]
mu   = 1e-6;                        % water kinematic viscosity (nu/rho)
DRho = DRhob;                          % difference in density between cells and water [kg/m^3]
g    = 10;                          % acceleration due to gravity [m/s^2]
eps  = logspace(-10,0,11);          % ocean energy dissipation rate, [W/kg]
r0   = [0.5e-8, 0.5e-6, 0.5e-4];    % particle radius in meters
rv   = logspace(-7,-2,1000);        % varying microbe size in [m]

% Choose some values to plot on the main plot:
EpsPlotVal = 6;                 % this is the chosen value to plot for epsilon
eps(EpsPlotVal)

% Loop through different particle sizes:
G_D1d = zeros(length(r0), length(rv));
G_B1d = G_D1d;
G_T1d = G_D1d;
for dd = 1:length(r0)

    % Diffusive Kernel:
    D0 = kT/(6*pi*nu*r0(dd));
    Dv = kT./(6*pi*nu*rv);
    G_D1d(dd,:) = 4*pi*(D0 + Dv) .* (r0(dd) + rv);
    
    % Buoyancy Kernel:
    U0 = 2 * DRho * g * r0(dd)^2 / (9*nu);
    Uv = 2 * DRho * g * rv.^2 / (9*nu);
    G_B1d(dd,:) = pi * (r0(dd) + rv).^2 .* abs(U0-Uv);
    
    % Turbulence Kernel, chosen for a particular epsilon value:
    G_T1 = 1.3 * (r0(dd) + rv).^3 .* sqrt(eps(EpsPlotVal)./mu);
    G_T2 = 1.37*pi.*eps(EpsPlotVal)^(1/3) * (r0(dd)+rv).^(7/3);
    Lk = 0.5*(mu^3./eps).^(1/4);
    ix1 = find(r0(dd)+rv < Lk(EpsPlotVal));
    ix2 = find(r0(dd)+rv >= Lk(EpsPlotVal));
    G_T1d(dd,:) = [G_T1(ix1), G_T2(ix2)];
    
    % Net kernel:
    G_N1d = G_D + G_B + [G_T1(ix1), G_T2(ix2)];
end

% Make interesting zone to shade:
InterestingZoneLims = 10;
LkLo = Lk./InterestingZoneLims;
LkHi = 10*Lk;
LineLims = [1e-19, 1e-5];
Xpatchvals = 1e6*[Lk(EpsPlotVal)/InterestingZoneLims, InterestingZoneLims*Lk(EpsPlotVal), InterestingZoneLims*Lk(EpsPlotVal), Lk(EpsPlotVal)/InterestingZoneLims];
Ypatchvals = Conv2CM*[LineLims(1), LineLims(1), LineLims(2), LineLims(2)];




% FIGURE ------------------------------------------------------------------
figure('units','centimeters','position',[3,3,big_width,12]);

% PLOT A ------------------------------------------------------------------
% Heatmap scaling:
ax3_label = axes('position',[0,.48,.08,.50]);
imshow(IMG_SPAN);

% Heatmap
ax3 = axes('position',[.13,.48,.34,.50]);
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','log');
set(gca,'yscale','log');
imagesc(1e6*rb, 1e6*rr, log10(G_N));

% Plot interesting contours
% y = x
plot(1e6*rr, 1e6*rr, '--','linewidth',1,'color','k');
% Diffusion-dominant
plot(X_DIFF, Y_DIFF, '-.','linewidth',1,'color','k');

% Lines depicting cross-sections shown in B:
line([1e6*min(rb), 1e6*max(rb)], [0.5,0.5], 'linestyle',':','color','k','linewidth',0.5);
line([1e6*min(rb), 1e6*max(rb)], [50,50], 'linestyle',':','color','k','linewidth',0.5);
line([1e6*min(rb), 1e6*max(rb)], [0.005,0.005], 'linestyle',':','color','k','linewidth',0.5);

% Labels:
xlabel('Microbe Size [\mum]');
ylabel('Resource Size [\mum]');
set(gca,'layer','top');
set(gca,'fontsize',7);
xlim(1e6*[1e-7,max(rb)]);
ylim(1e6*[min(rr),max(rr)]);
xticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]);
yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]);

% Colors:
MinRange = min(log10(G_N),[],'all');
MaxRange = max(log10(G_N),[],'all');
Ncolors  = 2*(ceil(MaxRange)-floor(MinRange));
CMap = turbo(Ncolors);
clim([floor(MinRange), ceil(MaxRange)]);
colormap(ax3,CMap);
cb = colorbar('AxisLocationMode','manual','position',[.48,.48,.015,.50]);
set(cb,'linewidth',1);
cb.Label.String = 'log10(Encounter Kernel [cm^3/s])';

% Annotations:
A1 = annotation('textbox',[.20,.78,.10,.10],'string','I','fontsize',10,'fontweight','bold','edgecolor','none');
A2 = annotation('textbox',[.20,.55,.10,.10],'string','II','fontsize',10,'fontweight','bold','edgecolor','none');
A3 = annotation('textbox',[.40,.75,.10,.10],'string','III','fontsize',10,'fontweight','bold','edgecolor','none');
A4 = annotation('textarrow',[.35,.39],[.86,.85],'string','y=x',...
    'fontsize',8,'color','k','headstyle','none','linestyle','none','textrotation',30);
A5 = annotation('textarrow',[.31,.31],[.55,.60],'string','Diffusion',...
    'fontsize',8,'color','k','headstyle','none','linestyle','none','textrotation',-55);

% PLOT B ------------------------------------------------------------------
% 1D Encounter Kernels:
% 100 um diameter resource
axBi = axes('position',[.63,.82,.36,.16]); 
hold on; box on; set(gca,'linewidth',1);
patch(Xpatchvals, Ypatchvals, [.8,.8,.8],'facealpha',1,'edgecolor','none');
% p0=line(1e6*[Lk(EpsPlotVal),Lk(EpsPlotVal)], Conv2CM*[LineLims],'linestyle','--','color','k','linewidth',1);

% Plot interesting contours from A:
line(.5*[1e2,1e2], Conv2CM*[LineLims],'linestyle','--','color','k','linewidth',1); % y = x
% line([Y_DIFF(ix100), Y_DIFF(ix100)], Conv2CM*[LineLims],'linestyle','-.','color','k','linewidth',1); % Diffusion-dominant

p1=plot(1e6*rv, Conv2CM*G_D1d(3,:), '-','linewidth',2,'color',[.2,.4,.8]);
p2=plot(1e6*rv, Conv2CM*G_B1d(3,:), '-','linewidth',2,'color',[.9,.7,.2]);
p3=plot(1e6*rv, Conv2CM*G_T1d(3,:), '-','linewidth',2,'color',[.8,.2,.2]);
set(gca,'layer','top');
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([1e-1, 2e2]);
ylim(Conv2CM*[1e-18, 1e-10]);
set(gca,'fontsize',7);
set(gca,'ticklength',[.02,.02]);
set(gca,'xticklabel',[]);
yticks([1e-11, 1e-8, 1e-5]);
set(gca,'YMinorTick','Off');

% 1 micron diameter resource
axBii = axes('position',[.63,.66,.36,.16]);
hold on; box on; set(gca,'linewidth',1);
patch(Xpatchvals, Ypatchvals, [.8,.8,.8],'facealpha',1,'edgecolor','none');
% p0=line(1e6*[Lk(EpsPlotVal),Lk(EpsPlotVal)], Conv2CM*[LineLims],'linestyle','--','color','k','linewidth',1);
line(.5*[1e0,1e0], Conv2CM*[LineLims],'linestyle','--','color','k','linewidth',1);
% line([Y_DIFF(ix010), Y_DIFF(ix010)], Conv2CM*[LineLims],'linestyle','-.','color','k','linewidth',1); % Diffusion-dominant

p1=plot(1e6*rv, Conv2CM*G_D1d(2,:), '-','linewidth',2,'color',[.2,.4,.8]);
p2=plot(1e6*rv, Conv2CM*G_B1d(2,:), '-','linewidth',2,'color',[.9,.7,.2]);
p3=plot(1e6*rv, Conv2CM*G_T1d(2,:), '-','linewidth',2,'color',[.8,.2,.2]);
set(gca,'layer','top');
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([1e-1, 2e2]);
ylim(Conv2CM*[1e-18, 1e-10]);
ylabel('Encounter kernel [cm^3/s]')
set(gca,'fontsize',7);
set(gca,'ticklength',[.02,.02]);
set(gca,'xticklabel',[]);
yticks([1e-11, 1e-8, 1e-5]);
set(gca,'YMinorTick','Off');

% 0.01 um diameter resource
axBiii = axes('position',[.63,.50,.36,.16]);
hold on; box on; set(gca,'linewidth',1);
patch(Xpatchvals, Ypatchvals, [.8,.8,.8],'facealpha',1,'edgecolor','none');
% p0=line(1e6*[Lk(EpsPlotVal),Lk(EpsPlotVal)], Conv2CM*[LineLims],'linestyle','--','color','k','linewidth',1);
line(.5*[1e-2,1e-2], Conv2CM*[LineLims],'linestyle','--','color','k','linewidth',1);
% line([Y_DIFF(ix001), Y_DIFF(ix001)], Conv2CM*[LineLims],'linestyle','-.','color','k','linewidth',1); % Diffusion-dominant

p1=plot(1e6*rv, Conv2CM*G_D1d(1,:), '-','linewidth',2,'color',[.2,.4,.8]);
p2=plot(1e6*rv, Conv2CM*G_B1d(1,:), '-','linewidth',2,'color',[.9,.7,.2]);
p3=plot(1e6*rv, Conv2CM*G_T1d(1,:), '-','linewidth',2,'color',[.8,.2,.2]);
set(gca,'layer','top');
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([1e-1, 2e2]);
ylim(Conv2CM*[1e-18, 1e-10]);
xlabel('Microbe Radius [\mum]');
set(gca,'fontsize',7);
set(gca,'ticklength',[.02,.02]);
yticks([1e-11, 1e-8, 1e-5]);
set(gca,'YMinorTick','Off');

% Legend:
legend([p1,p2,p3],{'Diff.','Buoy.','Turb.'},'location','northwest');

% Annotations:
annotation('textbox',[.63,.93,.36,.05],'string','100\mum diam.','fontsize',7,'edgecolor','none','horizontalalignment','center');
annotation('textbox',[.63,.77,.36,.05],'string','1\mum diam.','fontsize',7,'edgecolor','none','horizontalalignment','center');
annotation('textbox',[.63,.61,.36,.05],'string','0.01\mum diam.','fontsize',7,'edgecolor','none','horizontalalignment','center');


%{
patch(Xpatchvals, Ypatchvals, [.8,.8,.8],'facealpha',1,'edgecolor','none');
p0=line(1e6*[Lk(EpsPlotVal),Lk(EpsPlotVal)], Conv2CM*[LineLims],'linestyle','--','color','k','linewidth',1);
p1=plot(1e6*rv, Conv2CM*G_D, '-','linewidth',3,'color',[.2,.4,.8]);
p2=plot(1e6*rv, Conv2CM*G_B, '-','linewidth',3,'color',[.9,.7,.2]);
p3=plot(1e6*rv(ix1), Conv2CM*G_T1(ix1), '-','linewidth',3,'color',[.8,.2,.2]);
p4=plot(1e6*rv(ix2), Conv2CM*G_T2(ix2), '-','linewidth',3,'color',[.5,.2,.2]);
% plot(1e6*rv, Conv2CM*G_N1d, '-','linewidth',2,'color','k');
set(gca,'layer','top');
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([1e-1, 1e4]);
ylim(Conv2CM*[1e-18, 1e-5]);
xlabel('Microbe Size [\mum]');
ylabel('Encounter kernel [cm^3/s]')
set(gca,'fontsize',7);
set(gca,'ticklength',[.03,.03]);
%}

% % Inset: Interesting fluid regime
% ax4i = axes('position',[.66,.75,.14,.20]);
% hold on; box on; set(gca,'linewidth',1);
% patch([eps,fliplr(eps)],1e6*[LkLo, fliplr(LkHi)],[.8,.8,.8],'edgecolor','none','facealpha',1);
% plot(eps, 1e6*Lk, 'k--','linewidth',1);
% line([min(eps),max(eps)], [100, 100], 'linestyle','-','linewidth',1,'color','k');
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% set(gca,'layer','top');
% set(gca,'ticklength',[.03,.03]);
% xlabel('\epsilon [W/kg]');
% ylabel('\eta [\mum]');
% xlim([min(eps), max(eps)]);
% set(gca,'fontsize',7);

% PLOT C (left) -----------------------------------------------------------
axCL = axes('position',[0,0.02,.45,.38]);
% hold on; box on; set(gca,'linewidth',1); set(gca,'layer','top');
imshow(IMG_HOW);
T = annotation('textbox');
T.Position = [.30,.10,.10,.04];
T.String = 'Sodium alginate';
T.FontSize = 8;
T.EdgeColor = 'none';
T.BackgroundColor = 'none';

% PLOT C (right) ----------------------------------------------------------
axCR = axes('position',[.41,0.02,.25,.38]);
imshow(IMG_MONTAGE);

% PLOT D ------------------------------------------------------------------
% Histogram of sizes:
axD = axes('position',[.72,.07,.27,.32]);
hold on; box on; set(gca,'linewidth',1);
plot(Bins24, Counts24, '-','linewidth',2,'color',[.3,.3,.3]);
plot(Bins48, Counts48, '-','linewidth',2,'color',[.6,.6,.6]);
xlabel('Group Diameter [\mum]');
ylabel('PDF');
set(gca,'yscale','linear');
set(gca,'layer','top');
set(gca,'fontsize',7);
xlim([0,200]);
xticks([0,25,50,75,100,125,150,175]);
yticks([0,.005,.015,.025]);
legend({'24 hrs','48 hrs'},'location','northeast'); legend('boxoff');

% -------------------------------------------------------------------------
% Panel Labels:
annotation('textbox',[.10,.90,.10,.10],'string','a','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.53,.90,.10,.10],'string','b','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.02,.33,.10,.10],'string','c','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
% annotation('textbox',[.37,.33,.10,.10],'string','d','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.65,.33,.10,.10],'string','d','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')

% -------------------------------------------------------------------------
% Print figure:
if PRINT_FIGURES
    print([printfolder,'Figure_01'],'-dpng',res);
    exportgraphics(gcf,'Figure_01.pdf');
    close(gcf);
end


%% Figure 2:
% Mean encounter rates scale with size like a power law.

% Inputs:
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;
FilterFLPerArea = 5e3;
MinSize         = 4;
MaxSize         = 50;
FilterSize      = [MinSize, MaxSize];
MicrobeadDiam   = 1;
img_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\Imgs\'];
fig_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\'];
img_c = imread([img_folder, '12B01_3d_Aggregate_5umscalebar.png']);
% img_b = imread([img_folder, 'Microbeads_and_Boundaries_Example_50umscalebar.png']);
img_b = imread([img_folder,'Glamour_Shot_Zoomout.png']);
img_a = imread([img_folder,'Composite-1.png']);

% Crop the images a bit:
img_b_MicronsPerPixel = 1004./size(img_b,1);
img_b = permute(img_b,[2,1,3]);
img_b = fliplr(img_b(10:7e3,:,:));
% (2e3:7e3,2e3:7e3,:);
% img_a = img_a(1:1088, 11:1098,:);

% Loading and analysis ----------------------------------------------------
% Load a single example:
X = LOAD_DATA_BIN_CALIBRATE([fig_folder, '02_Timelapse\A\', 'A_090mins_001_EDF.nd2_TD_c_Probabilities.mat'], FilterFLPerArea, FilterSize, MicrobeadDiam);

% Load a range of concentrations:
Cbins = [12,7,2];
Clist = dir([fig_folder,'01_ChangingConcentration\New_Data\','*.mat']);
for cc = 1:length(Cbins)
    Y(cc) = LOAD_DATA_BIN_CALIBRATE([fig_folder,'01_ChangingConcentration\New_Data\',Clist(Cbins(cc)).name], FilterFLPerArea, FilterSize, MicrobeadDiam);
end

% Load a range of shaking speeds
RPMs = [50,100,150];
SSlist = dir([fig_folder,'04_ChangingShakingSpeed\Newest_Data\Mins_060\','*.mat']);
Slope = zeros(length(RPMs),1);
SlopeE = Slope;
for ss = 1:length(SSlist)
    Z(ss) = LOAD_DATA_BIN_CALIBRATE([fig_folder,'04_ChangingShakingSpeed\Newest_Data\Mins_060\',SSlist(ss).name], 5e3, FilterSize, MicrobeadDiam);
    Slope(ss) = Z(ss).B_Particles(2);
    SlopeE(ss) = Z(ss).m_std_error;
    PreFactor(ss) = Z(ss).B_Particles(1);
end

% Load timetracks:
t_list = dir([fig_folder,'02_Timelapse\A\','*.mat']);
Time = [0,30,60,90,120,245];

% Track a few bins:
ix_R05 = 8;
ix_R10 = 18;
ix_R15 = 24;
ix_R20 = 28;
ix_track = [ix_R05; ix_R10; ix_R15; ix_R20];

% Within tracked bins, get the average number of particles attached.
P = zeros(length(ix_track), length(t_list));
for tt = 1:length(t_list)
    W = LOAD_DATA_BIN_CALIBRATE_TIMETRACKS([fig_folder,'02_Timelapse\A\',t_list(tt).name], FilterFLPerArea, MinSize, 1);
    disp(length(W.AggregateRadiusFilt));
    disp(t_list(tt));
    P(:,tt) = W.Avg_Rad(ix_track, 3);
    S(:,tt) = W.Avg_Rad(ix_track, 4)./sqrt(W.N_Rad(ix_track));
end
Radii_avg = W.Avg_Rad(ix_track,1);
Radii_err = W.Avg_Rad(ix_track,2);

% Set up the theoretical result:
modelfun = @(parameters,t) ( 1/parameters(2) * parameters(1) * 1e7 * (1 - exp(-parameters(2)*t)) );

% % Previous regression results:
% FitTable = [.54, .063, 1.1e-8; ...
%             .41, .020, 2.5e-8; ...
%             .47, .020, 7.5e-8; ...
%             .51, .015, 1.0e-7];

% Do a regression to the model:
gamma_0 = 1e-8*ones(length(ix_track));
beta_0  = 0.01;
fit_params = zeros(2, length(ix_track));
for ii = 1:length(ix_track)
    parameters0 = [gamma_0(ii), beta_0];
    mdl = fitnlm(Time, P(ii,:), modelfun, parameters0);
    fit_params(:,ii) = mdl.Coefficients.Estimate;
end
FitTable = fit_params';

% FIGURE ------------------------------------------------------------------
figure('units','centimeters','position',[3,3,big_width,12]);

% PLOT A: -----------------------------------------------------------------
% Microbeads and aggregates in suspension, and confocal image:
ax_a = axes('position',[0,.50,.33,.50]);
imshow(img_b); hold on;
x_loc = 100; y_loc = 6.7e3; ScaleBarLength = 100;
quiver(x_loc, y_loc, ScaleBarLength/img_b_MicronsPerPixel, 0, 'ShowArrowHead', 'off','color','k','linewidth',3);

% Box annotation:
T = annotation('rectangle');
T.Position = [.2325,.608,.02,.03];
T.EdgeColor = 'k';
T.LineWidth = 1;
annotation('line',[.2525,.34],[.608,.55],'linewidth',1,'color','k');
annotation('line',[.2525,.34],[.638,.77],'linewidth',1,'color','k');
T = annotation('rectangle');
T.Position = [.34,.55,.15,.22];
T.EdgeColor = 'k';
T.LineWidth = 1;

% One aggregate confocal image:
ax_bt = axes('position',[.34,.77,.15,.23]);
imshow(img_c);

% Zoomed-in image of an aggregate:
ax_bb = axes('position',[.34,.55,.15,.22]);
imshow(imrotate(img_a,-90));

% PLOT B: -----------------------------------------------------------------
% Power law regression for one dataset, log-log scale:
ax_c = axes('position',[.55,.56,.32,.43]);
hold on; box on; set(gca,'linewidth',1);
plot(X.AggregateRadiusFilt, X.NumParticlesFilt, 'o','markersize',3,'markerfacecolor',[.65,.65,.65],'markeredgecolor','none');
errorbar(X.Avg_Particles(:,1), X.Bins_Particles, 0.6, 0.6, X.Avg_Particles(:,2), X.Avg_Particles(:,2), 'd','markersize',6,'markerfacecolor','k','markeredgecolor','none','linewidth',1,'color','k');
p3=plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',1,'color',[.8,.2,.2]);
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([3,40]);
ylim([5e-1, 1e2]);
xticks([1,5,10,30,100]);
yticks([1,3,10,30,100]);
set(gca,'fontsize',7);
set(gca,'TickLength',[.02,.02]);
set(gca,'layer','top');
xlabel('Group Radius [\mum]');
ylabel('# Microbeads');

% Text additions:
legend([p3],{'p = A(r_a + r_b)^{\lambda}'},'location','northwest','fontsize',8);
legend('boxoff');
text(8,40,['\lambda = ',num2str(X.B_Particles(2),2),' \pm ',num2str(3*X.m_std_error,2)],'fontsize',6);

% Side panel, histogram of particles attached to aggregates:
% ax_ci = axes('position',[.87,.62,.11,.12]);
ax_ci = axes('position',[.87,.56,.12,.43]);
hold on; box on; set(gca,'linewidth',1);
barh(X.Bins_Particles, X.N_Particles, 1, 'facecolor',[.65,.65,.65],'edgecolor','none','facealpha',1,'linewidth',1);
xlabel('Counts');
set(gca,'fontsize',7);
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([5e-1,2e3]);
ylim([5e-1,100]);
xticks([1,10,100,1e3]);
yticks([1,3,10,30,100]);
yticklabels({'','','','',''});
set(gca,'ticklength',[.03,.03]);
set(gca,'layer','top');

% PLOT C: -----------------------------------------------------------------
ax_d = axes('position',[.05, .06, .30, .42]);
hold on; box on; set(gca,'linewidth',1);
Colors = abyss(length(Cbins));
for cc = 1:length(Cbins)
    errorbar(Y(cc).Avg_Particles(:,1), Y(cc).Bins_Particles, 0.6, 0.6, Y(cc).Avg_Particles(:,2), Y(cc).Avg_Particles(:,2), 'd','markersize',6,'markerfacecolor',Colors(cc,:),'markeredgecolor','k','linewidth',1,'color','k');
    ph(cc) = plot(Y(cc).X_Particles(:,2), Y(cc).R_Particles, '-','linewidth',1,'color',Colors(cc,:));
end
xlabel('Group radius [um]');
ylabel('No. Beads');
xlim([3, 40]);
xticks([5,10,30]);
ylim([6e-1, 100]);
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'layer','top');
set(gca,'ticklength',[.02,.02]);
set(gca,'fontsize',7);
legend(ph,{'C=1e6 [beads/mL]','C=5e6','C=1e7'},'location','none','fontsize',7,'position',[.10,.32,.10,.10]);
legend('boxoff');

% PLOT D: -----------------------------------------------------------------
ax_e = axes('position',[.40,.06,.33,.42]);
hold on; box on; set(gca,'linewidth',1);
Colors = bone(length(ix_track)+1);
for ii = 1:length(ix_track)
    y_reg  = modelfun(FitTable(ii,:), Time);
    t_plot = linspace(0,max(Time),1e3);
    y_plot = modelfun(FitTable(ii,:), t_plot);
    p(ii) = errorbar(Time, P(ii,:), S(ii,:), 's-','markersize',7,'linewidth',1,'markerfacecolor',Colors(ii,:), 'markeredgecolor','none','color',Colors(ii,:),'capsize',6);
    % p(ii)  = plot(Time, P(ii,:), 's-','markersize',7,'linewidth',2,'color',Colors(ii,:),'markerfacecolor',Colors(ii,:), 'markeredgecolor','none');
    plot(t_plot, y_plot, 'k-','linewidth',1);
end
set(gca,'layer','top');
xlabel('Time [min]');
ylabel('No. Beads');
LegendLabels = {};
for ii = 1:length(ix_track)
    LegendLabels{ii} = [num2str(Radii_avg(ii),3), ' \pm ', num2str(Radii_err(ii),2), '\mum'];
end
LegendLabels
legend([p(1), p(2), p(3), p(4)], LegendLabels,'location','none','fontsize',6,'position',[.42,.35,.15,.10]);
legend('boxoff');
set(gca,'fontsize',7);

% PLOT E: -----------------------------------------------------------------
Colors = jet(3);
ax_f = axes('position',[.78,.06,.20,.42]);
hold on; box on; set(gca,'linewidth',1);
for ss = 1:length(Z)
    bar(ss, Slope(ss), 0.75, 'facecolor',Colors(ss,:), 'edgecolor','k','linewidth',1,'facealpha',1);
end
errorbar(1:length(SSlist), Slope, 2*SlopeE, '.','markersize',1,'color','k','linewidth',1,'capsize',6);
xlim([.5, 3.5]);
ylim([2, 3]);
xlabel('Shaking Speed [rpm]');
ylabel('Power Law Value, \lambda');
xticks([1,2,3]);
xticklabels({'50','100','150'});
yticks([2,2.5,3])
set(gca,'layer','top');
set(gca,'fontsize',7);

% Inset: tracking epsilon with shaking speed:
%{
ax_fi = axes('position',[.78,.32,.20,.16]);
hold on; box on; set(gca,'linewidth',1);
% plot(n0, eps0, 'k-','linewidth',1);
% for ii = 1:length(RPMs)
%     plot(RPMs(ii), eps(ii), 'o','markersize',8,'markerfacecolor',Colors(ii,:),'markeredgecolor','k','linewidth',1);
% end
set(gca,'yscale','log');
set(gca,'fontsize',7);
ylim([1e-7,1e-2]);
yticks([1e-6,1e-4,1e-2]);
set(gca,'ticklength',[.03,.03]);
xlabel('Shaking speed [rpm]');
ylabel('\epsilon [W/kg]')
%}

% -------------------------------------------------------------------------
% Panel Labels:
annotation('textbox',[.01,.92,.03,.10],'string','a','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.52,.92,.10,.10],'string','b','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.05,.40,.05,.10],'string','c','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.365,.40,.10,.10],'string','d','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.75,.40,.10,.10],'string','e','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')

% Print:
% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Figure_02'],'-dpng',res);
    exportgraphics(gcf,'Figure_02.pdf');
    close(gcf);
end

%% Figure 3:
% Fluctuations in encounters are predicted by Poisson considerations.

% Inputs:
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;
FilterFLPerArea = 5e3;
FilterSize      = [4,Inf];
MicroBeadDiam   = 1;
fig_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\02_Timelapse\A\'];
X = LOAD_DATA_BIN_CALIBRATE([fig_folder, 'A_120mins_001_EDF.nd2_TD_c_Probabilities.mat'], ...
    FilterFLPerArea, FilterSize, MicroBeadDiam);

% Generate Poisson predictions:
[f, PCPoiss, chi2] = DO_POISSON_STUFF(X, 3, 0);
figure('units','centimeters','position',[3,3,lil_width,12]);

% PLOT A: -----------------------------------------------------------------
% Stacked Plot:
AddConstant = 0.075;
YTicks = [];
YLabels = cell(1,1);
ax_a = axes('position',[.09,.40,.53,.59]);
hold on; box on; set(gca,'linewidth',1);
x = [];
y = [];

% Generate list of ticks, colors, and labels:
YTicks = [];
for bb = 1:length(PCPoiss)
    if ~isempty(PCPoiss{bb})
        YTicks = [YTicks, AddConstant*bb];
        YLabels{length(YTicks)} = [num2str(X.Avg_Rad(bb,1),2), ' \mum'];
    end
end

% Generate colors:
Colors = turbo(length(YLabels));

% Plot the curves:
bb_list_temp = [11,21,31,41];

for ii = 1:2:length(PCPoiss)
    bb = ii;
    if ~isempty(PCPoiss{bb})
        plot(f{bb}(:,1), f{bb}(:,2) + AddConstant*bb, 'o-','linewidth',2,'markersize',4,'color',Colors(bb,:),'markerfacecolor',Colors(bb,:),'markeredgecolor','none');
        plot(PCPoiss{bb}(1,:), PCPoiss{bb}(2,:) + AddConstant*bb, '-','linewidth',1,'markersize',12,'color','k');

        % Find the peak of the Poisson distribution for each size bin:
        [val, ix] = max(PCPoiss{bb}(2,:));
        x = [x, PCPoiss{bb}(1,ix)];
        y = [y, val + AddConstant*bb];
    end
end

% Colormap:
colormap(ax_a, Colors);
clim([1,length(YTicks)]);
cb = colorbar('manual');
cb.Position = [.63,.407,.02,.56];
cb.Ticks = 1:5:length(YTicks);
cb.TickLabels = YLabels(1:5:end);
cb.LineWidth = 1;

% Axis scaling and visualization:
set(gca,'xscale','linear');
xlim([0,50]);
ylim([0,AddConstant*(length(Colors)+1)]);
ylabel('PDF');
xlabel('# Microbeads');
set(gca,'fontsize',7);
yticks(YTicks);
yticklabels({});
set(gca,'layer','top');
set(gca,'yscale','linear');
set(gca,'ticklength',[.01,.05]);

% Outside: some panels
ax_a_i = axes('position',[.75,.40,.24,.12]);
hold on; box on; set(gca,'linewidth',1);
bbi=bb_list_temp(1);
plot(f{bbi}(:,1), f{bbi}(:,2), '.-','linewidth',2,'markersize',12,'color',Colors(bbi,:));
plot(PCPoiss{bbi}(1,:), PCPoiss{bbi}(2,:), '-','linewidth',1,'markersize',12,'color','k');
xlim([0,6]);
ylim([0,1.3*max(f{bbi}(:,2))]);
set(gca,'layer','top');
set(gca,'fontsize',7);
xticks([0,2,4,6]);
yticks([]);
% xticklabels({});
set(gca,'ticklength',[.03,.03]);

ax_a_ii = axes('position',[.75,.563,.24,.12]);
hold on; box on; set(gca,'linewidth',1);
bbii=bb_list_temp(2);
plot(f{bbii}(:,1), f{bbii}(:,2), '.-','linewidth',2,'markersize',12,'color',Colors(bbii,:));
plot(PCPoiss{bbii}(1,:), PCPoiss{bbii}(2,:), '-','linewidth',1,'markersize',12,'color','k');
xlim([0,12]);
ylim([0,1.3*max(f{bbii}(:,2))]);
set(gca,'layer','top');
set(gca,'fontsize',7);
xticks([0:3:10]);
% xticklabels({});
yticks([]);
set(gca,'ticklength',[.03,.03]);

ax_a_iii = axes('position',[.75,.717,.24,.12]);
hold on; box on; set(gca,'linewidth',1);
bbiii=bb_list_temp(3);
plot(f{bbiii}(:,1), f{bbiii}(:,2), '.-','linewidth',2,'markersize',12,'color',Colors(bbiii,:));
plot(PCPoiss{bbiii}(1,:), PCPoiss{bbiii}(2,:), '-','linewidth',1,'markersize',12,'color','k');
xlim([0,30]);
ylim([0,1.4*max(f{bbiii}(:,2))]);
set(gca,'layer','top');
set(gca,'fontsize',7);
xticks([0:8:30]);
% xticklabels({});
yticks([]);
set(gca,'ticklength',[.03,.03]);

ax_a_iv = axes('position',[.75,.87,.24,.12]);
hold on; box on; set(gca,'linewidth',1);
bbiv=bb_list_temp(4);
plot(f{bbiv}(:,1), f{bbiv}(:,2), '.-','linewidth',2,'markersize',12,'color',Colors(bbiv,:));
plot(PCPoiss{bbiv}(1,:), PCPoiss{bbiv}(2,:), '-','linewidth',1,'markersize',12,'color','k');
xlim([0,40]);
ylim([0,1.4*max(f{bbiv}(:,2))]);
set(gca,'layer','top');
set(gca,'fontsize',7);
xticks([0:10:35]);
xticklabels({'0','10','20','30','40'});
yticks([]);
set(gca,'ticklength',[.03,.03]);

% Add annotation as axis label:
T = annotation('textbox');
T.String = '# Microbeads';
T.FontSize = 7;
T.Position = [.74,.34,.26,.03];
T.BackgroundColor = [1,1,1];
T.EdgeColor = 'none';
T.HorizontalAlignment = 'center';
T.VerticalAlignment = 'middle';

% Add annotations as size class labels:
Ti = annotation('textbox');
Ti.FontSize = 7;
Ti.Position = [.74,.47,.26,.034];
Ti.BackgroundColor = 'none';
Ti.EdgeColor = 'none';
Ti.HorizontalAlignment = 'center';
Ti.VerticalAlignment = 'middle';
Ti.String = YLabels{bbi};

Tii = annotation('textbox');
Tii.FontSize = 7;
Tii.Position = [.74,.64,.26,.034];
Tii.BackgroundColor = 'none';
Tii.EdgeColor = 'none';
Tii.HorizontalAlignment = 'center';
Tii.VerticalAlignment = 'middle';
Tii.String = YLabels{bbii};

Tiii = annotation('textbox');
Tiii.FontSize = 7;
Tiii.Position = [.74,.80,.26,.034];
Tiii.BackgroundColor = 'none';
Tiii.EdgeColor = 'none';
Tiii.HorizontalAlignment = 'center';
Tiii.VerticalAlignment = 'middle';
Tiii.String = YLabels{bbiii};

Tiv = annotation('textbox');
Tiv.FontSize = 7;
Tiv.Position = [.74,.95,.26,.034];
Tiv.BackgroundColor = 'none';
Tiv.EdgeColor = 'none';
Tiv.HorizontalAlignment = 'center';
Tiv.VerticalAlignment = 'middle';
Tiv.String = YLabels{bbiv};

% PLOT B: -----------------------------------------------------------------
ix_show = find(X.N_Rad > 3);
ax_b = axes('position',[.09,.07,.42,.26]);
hold on; box on; set(gca,'linewidth',1);
p1=plot(X.Avg_Rad(ix_show,1), X.Avg_Rad(ix_show,3), 'd','markersize',5,'markeredgecolor','none','markerfacecolor','k');
p2=plot(X.Avg_Rad(ix_show,1), X.Avg_Rad(ix_show,4).^2, 's','markersize',5,'markeredgecolor','none','markerfacecolor',[.2,.4,.8]);
xlabel('Radius [\mum]');
ylabel('Measure');
set(gca,'xscale','linear');
set(gca,'yscale','log');
hleg=legend([p1,p2],{'\langle P\rangle','\sigma^2'},'location','southeast'); legend('boxoff');
set(gca,'fontsize',7);
set(gca,'ticklength',[.02,.02]);
xlim([3,30]);
ylim([3e-1,3e2]);
xticks([5:5:30]);
yticks([1,3,10,30,100]);

% PLOT C: -----------------------------------------------------------------
var_by_mean = X.Avg_Rad(ix_show,4).^2./X.Avg_Rad(ix_show,3);
% test_statistic = (X.Avg_Rad(ix_show,4).^2 - X.Avg_Rad(ix_show,3))./X.Avg_Rad(ix_show,3);
COV = X.Avg_Rad(ix_show,4)./X.Avg_Rad(ix_show,3);
test_statistic = var_by_mean;
ix_outside = find(test_statistic < 0.5 | test_statistic > 2);
ix_inside = find(test_statistic >= 0.5 & test_statistic <= 2);
ax = axes('position',[.59,.07,.39,.26]);
hold on; box on; set(gca,'linewidth',1);
% plot(X.Avg_Rad(ix_show,3), sqrt(X.Avg_Rad(ix_show,3)),'k--','linewidth',1);
% plot(X.Avg_Rad(ix_show,3), X.Avg_Rad(ix_show,4),'o','markersize',3,'markerfacecolor','k','markeredgecolor','none');
plot([0,50],[.5,.5], 'k--','linewidth',1);
plot([0,50],[2,2], 'k--','linewidth',1);
plot(X.Avg_Rad(ix_show(ix_inside), 1), test_statistic(ix_inside), 'o','markersize',3,'markeredgecolor','k','linewidth',1);
plot(X.Avg_Rad(ix_show(ix_outside), 1), test_statistic(ix_outside), 'x','markersize',4,'markeredgecolor','r','linewidth',1);
xlim([3,30]);
ylim([0,10]);
xlabel('Radius [\mum]');
ylabel('\sigma^2/\langle P\rangle');
set(gca,'fontsize',7);
xticks([5:5:30]);
yticks(1:2:10);

% -------------------------------------------------------------------------
% Panel Labels:
annotation('textbox',[.03,.92,.03,.10],'string','a','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.11,.25,.03,.10],'string','b','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.61,.25,.03,.10],'string','c','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Figure_03'],'-dpng',res);
    exportgraphics(gcf,'Figure_03.pdf');
    close(gcf);
end

%% Figure 4:
% Mean encounters per capita decrease, but consistency in encounters per
% capita increases.

% Load data:
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;
fig_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\02_Timelapse\A\'];
reg_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\06_ConfocalSegmentations\'];
FilterFLPerArea = 5e3;
FilterSize      = [4,Inf];
MicroBeadDiam   = 1;
X = LOAD_DATA_BIN_CALIBRATE([fig_folder, 'A_090mins_001_EDF.nd2_TD_c_Probabilities.mat'], FilterFLPerArea, FilterSize, MicroBeadDiam); 
Y = load([reg_folder, 'DATA_Figure_04_inset.mat']);
Y.c_est;

% Use V vs. N regression to estimate the number of cells per aggregate in
% microbead data:
Y.c_est(2) = 0;
AggVolume = 4/3*pi*X.AggregateRadiusFilt.^3;
AggNcells = polyval(Y.c_est, AggVolume);

% Revision addition: Estimate the number of cells just in the peripheral
% shell of the group:
DeltaR = [1, 3, 6, 10]; % micron
RadiusCore = X.AggregateRadiusFilt - DeltaR;
ShellVolume = 4/3*pi * (X.AggregateRadiusFilt.^3 - RadiusCore.^3);
ShellNcells = polyval(Y.c_est, ShellVolume);
DR_Ill = imread(['G:\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\Imgs\','DeltaR_Ill_lumpy.png']);

% Accounting for error in the slope:
m_std_error = 0.021574;
m_min = Y.c_est(1) - m_std_error;
m_max = Y.c_est(1) + m_std_error;
AggNcells_min = m_min * AggVolume;
AggNcells_max = m_max * AggVolume;
N_error = AggNcells_max - AggNcells;

% Measure particles per capita:
PPerCapita = X.NumParticlesFilt./AggNcells;
PPerCapita_Shell = X.NumParticlesFilt./ShellNcells;

% Binning by aggregate size, get average and std particles per capita:
PPerCapita_Rad = zeros(length(X.Group_Rad),4);
PPerCapita_Shell_Rad = zeros(length(X.Group_Rad), 4, length(DeltaR));
for bb = 1:length(X.Group_Rad)
    x = X.Group_Rad{bb}(:,1);
    n = polyval(Y.c_est, 4/3*pi*x.^3);
    p = X.Group_Rad{bb}(:,2);
    for rr = 1:length(DeltaR)
        xc = x - DeltaR(rr);
        if xc > 0 % the shell width is smaller than the radius of the group
            ns = polyval(Y.c_est, 4/3*pi*(x.^3 - xc.^3));
        else % the shell width is the same size or larger than the group radius
            ns = nan; % we just won't plot the thing when dr > r
        end
        ps = p./ns; % estimate number of particles per capita, shell only
        PPerCapita_Shell_Rad(bb,:,rr) = [mean(x), std(x), mean(ps), std(ps)];
    end
    p = p./n; % estimated number of particles per capita
    PPerCapita_Rad(bb,:) = [mean(x), std(x), mean(p), std(p)];
end
ix = find(PPerCapita_Shell_Rad(:,1,1) >= 4);
for rr = 1:length(DeltaR)
    ppcrs(:,:,rr) = PPerCapita_Shell_Rad(ix, :, rr);
end
PPerCapita_Rad = PPerCapita_Rad(PPerCapita_Rad(:,1) >= 4, :);

% Estimating particles per capita in the regression:
Reg_AggVol = 4/3*pi*X.X_Rad(:,2).^3;
Reg_Ncells = polyval(Y.c_est, Reg_AggVol);
PPerCapita_Reg = X.R_Rad./Reg_Ncells;

% FIGURE ------------------------------------------------------------------

figure('units','centimeters','position',[3,3,lil_width,10]);

% PLOT A: -----------------------------------------------------------------
MinRad = min(X.AggregateRadiusFilt);
MaxRad = max(X.AggregateRadiusFilt);

skip = 2;
ax_a = axes('position',[.10,.48,.88,.49]);
hold on; box on; set(gca,'linewidth',1);
plot(X.AggregateRadiusFilt, PPerCapita,'o','markersize',2,'markerfacecolor',[.65,.65,.65],'markeredgecolor','none');
% plot(PPerCapita_Rad(:,1), PPerCapita_Rad(:,3), 'd','markersize',5,'markerfacecolor','k','markeredgecolor','none');
errorbar(PPerCapita_Rad(1:skip:end,1), PPerCapita_Rad(1:skip:end,3), PPerCapita_Rad(1:skip:end,4), PPerCapita_Rad(1:skip:end,4), PPerCapita_Rad(1:skip:end,2), PPerCapita_Rad(1:skip:end,2),...
    'kd','markersize',5,'markerfacecolor','k','markeredgecolor','none','linewidth',1,'capsize',4);
plot([MinRad, MaxRad],[0,0],'k:','linewidth',1);
xlim([MinRad-1, 40]);
ylim([-1e-3, 1e-2]);

% Labels:
xlabel('Radius [um]');
ylabel('Microbeads per cell');
set(gca,'xscale','linear');
set(gca,'yscale','linear');
set(gca,'fontsize',7);
set(gca,'layer','top');

% Inset, what if we are just talking about feeding the shell:
skip = 1;
Colormap_a = flipud([254,217,142; 254,153,41; 217,95,14; 153,52,4]/255);
ax_inset = axes('position',[.55,.68,.40,.25]);
hold on; box on; set(gca,'linewidth',1);
plot([1,40],[0,0],'k:','linewidth',1);
shadedErrorBar(PPerCapita_Rad(1:skip:end,1), PPerCapita_Rad(1:skip:end,3), PPerCapita_Rad(1:skip:end,4),...
        'lineProps',{'-o','markerfacecolor','none','markeredgecolor','none','linewidth',1,'color','k'});
shadedErrorBar(PPerCapita_Shell_Rad(1:skip:end,1,1), PPerCapita_Shell_Rad(1:skip:end,3,1), PPerCapita_Shell_Rad(1:skip:end,4,1),...
        'lineProps',{'-o','markerfacecolor','none','markeredgecolor','none','linewidth',1,'color',[232, 89, 12]/255});
set(gca,'xscale','log');
set(gca,'yscale','linear');
xlabel('Radius [\mum]');
ylabel('P/N');
xlim([MinRad, 35]);
ylim([-1e-3, 2e-2])
set(gca,'fontsize',7);
set(gca,'layer','top');

ax_legend = axes('position',[.33,.73,.13,.25]);
imshow(DR_Ill);
annotation("textbox",[.56,.80,.10,.05],'string','1\mum','edgecolor','none','fontsize',7,'color',[232, 89, 12]/255);
annotation("textbox",[.56,.69,.10,.05],'string','Full','edgecolor','none','fontsize',7,'color','k');

% PLOT B: -----------------------------------------------------------------
skip = 1;
XFit = linspace(min(PPerCapita_Rad(:,1)), max(PPerCapita_Rad(:,1)), 20);
ax_b = axes('position',[.10,.08,.39,.32]);
hold on; box on; set(gca,'linewidth',1);

yyaxis left;
p1=plot(PPerCapita_Rad(:,1), PPerCapita_Rad(:,3),'d','markersize',5,'markerfacecolor','k','markeredgecolor','none');
ylabel('Mean Beads/Capita');
ylim([3e-4,2e-2]);
set(gca,'xscale','log');
set(gca,'yscale','log');

yyaxis right;
p2=plot(PPerCapita_Rad(:,1), (PPerCapita_Rad(:,4)).^(-1), 's','markersize',5,'markerfacecolor',[.2,.4,.8],'markeredgecolor','none');
ylim([1e2,1e4]);
set(gca,'yscale','log');
xlim([MinRad,MaxRad]);
yticklabels({});

T = annotation('textbox','string','Consistency (1/\sigma)');
T.Position = [.15,.35,.50,.04];
T.FontSize = 7;
T.Color = [.2,.4,.8];
T.EdgeColor = 'none';

% Labels
xlabel('Radius [\mum]');
set(gca,'fontsize',7);
set(gca,'ticklength',[.03,.03]);
set(gca,'layer','top');

% Set axes colors:
ax_b.YAxis(1).Color = [0,0,0];
ax_b.YAxis(2).Color = [.2,.4,.8];

% PLOT C: -----------------------------------------------------------------
%{
skip = 1;
XFit = linspace(min(PPerCapita_Rad(:,1)), max(PPerCapita_Rad(:,1)), 20);
ax_c = axes('position',[.55,.39,.39,.22]);
hold on; box on; set(gca,'linewidth',1);
yyaxis left;
for rr = 1:length(DeltaR)
    % labels(rr) = scatter(PPerCapita_Shell_Rad(1:skip:end,1,rr), PPerCapita_Shell_Rad(1:skip:end,3,rr), 25, Colormap_a(rr,:), 'd');
    plot(PPerCapita_Shell_Rad(1:skip:end,1,rr), PPerCapita_Shell_Rad(1:skip:end,3,rr), '-','linewidth',1.5,'color',Colormap_a(rr,:));
end
plot(PPerCapita_Rad(:,1), PPerCapita_Rad(:,3),'k-','linewidth',1);
ylim([3e-4,2e-2]);
set(gca,'xscale','log');
set(gca,'yscale','log');
yticklabels({});

yyaxis right;
for rr = 1:length(DeltaR)
    % labels(rr) = scatter(PPerCapita_Shell_Rad(1:skip:end,1,rr), PPerCapita_Shell_Rad(1:skip:end,4,rr).^(-1), 25, Colormap_a(rr,:), 's');
    plot(PPerCapita_Shell_Rad(1:skip:end,1,rr), PPerCapita_Shell_Rad(1:skip:end,4,rr).^(-1), '-.','linewidth',1.5,'color',Colormap_a(rr,:))
end
plot(PPerCapita_Rad(:,1), PPerCapita_Rad(:,4).^(-1),'k-.','linewidth',1);
ylim([1e2,1e4]);
set(gca,'yscale','log');
xlim([MinRad,MaxRad]);

% Labels
xlabel('Radius [\mum]');
set(gca,'fontsize',7);
set(gca,'ticklength',[.03,.03]);
set(gca,'layer','top');

% Set axes colors:
ax_c.YAxis(1).Color = [0,0,0];
ax_c.YAxis(2).Color = [0,0,0];
%}

% PLOT D: -----------------------------------------------------------------
B_sweep = 1:.2:4;
B_colors = parula(length(B_sweep));
MarkerStyle = {'o','d','s','pentagram'};
ax_c = axes('position',[.57,.08,.39,.32]);
% ax_d = axes('position',[.10,.06,.88,.27]);
hold on; box on; set(gca,'linewidth',1);
% for rr = 1:length(DeltaR)
%     scatter(PPerCapita_Shell_Rad(1:skip:end,3,rr), PPerCapita_Shell_Rad(1:skip:end,4,rr).^(-1), 25,Colormap_a(rr,:),MarkerStyle{rr});
% end
plot(PPerCapita_Rad(1:end,3), (PPerCapita_Rad(1:end,4)).^(-1),'kx','linewidth',1,'markersize',6);
MaxY = 8000;
% plot([7e-4, max(PPerCapita_Rad(:,3))], [MaxY, 0], 'k--','linewidth',1);
xlabel('Mean Beads/Capita');
ylabel('Consistency (1/\sigma)');
set(gca,'fontsize',7);
set(gca,'layer','top');
set(gca,'xscale','linear');
set(gca,'yscale','linear');
xlim([0,5e-3]);
ylim([0,1e4]);
yticks([]);
xticks([]);

% Inset:
%{
ax_di = axes('position',[.50,.15,.40,.15]);
hold on; box on; set(gca,'linewidth',1);
for rr = 1:length(DeltaR)
    scatter(PPerCapita_Shell_Rad(1:skip:end,3,rr), PPerCapita_Shell_Rad(1:skip:end,4,rr).^(-1), 25,Colormap_a(rr,:),MarkerStyle{rr});
end
plot(PPerCapita_Rad(1:end,3), (PPerCapita_Rad(1:end,4)).^(-1),'kx','linewidth',1,'markersize',6);
xlabel('Mean Beads/Capita');
ylabel('Consistency (1/\sigma)');
set(gca,'fontsize',7);
set(gca,'layer','top');
set(gca,'xscale','log');
set(gca,'yscale','log');
%}

% -------------------------------------------------------------------------
% Diagram on top illustrating the shell width:
% DeltaR_Ill = imread([drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\Imgs\','DeltaR_Ill.png']);
% ax_ill = axes('position',[.25,.50,.15,.11]);
% hold on; box on; set(gca,'linewidth',1);
% imshow(DeltaR_Ill);
% set(gca,'layer','top');

% Panel Labels:
annotation('textbox',[.01,.92,.03,.10],'string','a','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.12,.54,.03,.10],'string','b','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.57,.54,.03,.10],'string','c','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.13,.25,.03,.10],'string','d','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Figure_04'],'-dpng',res);
    exportgraphics(gcf,'Figure_04.pdf');
    close(gcf);
end

%% Figure 5
% A model shows that there is a fitness trade-off between mean encounter
% rate and consistency of encounters.

% Load data:
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;
fig_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\07_Sims\'];
sims_ill = imread([fig_folder,'Sim_Ill_3.png']);
% X=load([fig_folder, 'g0_alpha_sims\', 'sims_FS_17-Jun-2025.mat']);
X = load([drive_folder,':\My Drive\02_Scripts\EncounterRates\States_Competition_Model\g0_alpha_sims\', 'sims_FS_23-Jan-2026.mat']);
Z=load([fig_folder, 'swimming_speed_sims\', 'sims_02-Mar-2025.mat']);
% W3=load([fig_folder, 'g0_alpha_sims\', 'sims_MF_17-Jun-2025.mat'],'G_0','alpha',...
    % 'n_frac','b_frac','FinalC_FS');
W3 = load([drive_folder,':\My Drive\02_Scripts\EncounterRates\States_Competition_Model\g0_alpha_sims\', 'sims_MF_23-Jan-2026.mat']);
M=load([fig_folder, 'growth_ratio_sims\', 'sims_11-Feb-2025.mat']);
V=load([fig_folder, 'packing_frac_sims\','sims_03-Mar-2025.mat']);

% Create a colormap for all heatmaps of the figure to share:
cmap = redwhiteblue(-.5, .5, 25);

% Converting swimming strength to swimming speed:
V0 = sqrt(3*Z.Ds);

% FIGURE ------------------------------------------------------------------
figure('units','centimeters','position',[3,3,lil_width,11.5]);

% PLOT A ------------------------------------------------------------------
% An illustration of the simulations
ax_a = axes('position',[0,.70,1,.30]);
imshow(sims_ill);

% PLOT B ------------------------------------------------------------------
% Show example runs, first n_frac...
ax_bi = axes('position',[.11,.40,.22,.14]);
hold on; box on; set(gca,'linewidth',1);
p1 = plot(W3.n_frac{4,9,1},'-','linewidth',1,'color',[.6,.6,.6]);
% plot(W3.n_frac{4,6,2},'-','linewidth',1,'color',[.6,.6,.6]);
% plot(W3.n_frac{4,6,3},'-','linewidth',1,'color',[.6,.6,.6]);
p2 = plot(X.n_frac{4,9,1},'-','linewidth',1,'color','k');
% plot(X.n_frac{4,9,2},'-','linewidth',1,'color','k');
% plot(X.n_frac{4,9,3},'-','linewidth',1,'color','k');

% Labels:
xlabel('Time');
ylabel('# Fraction');
set(gca,'fontsize',7);
xlim([0,50]);
ylim([0,1]);
set(gca,'layer','top');
xticks([0,10,20,30,40,50]);
yticks([0,.25,.5,.75]);
set(gca,'ticklength',[.03,.03]);

% Legend:
L2a = annotation('textbox','position',[.14,.50,.22,.04],'String','No Fluc',...
    'fontsize',6,'edgecolor','none');
L2b = annotation('line',[.13,.15],[.517,.517],'linewidth',1,'color',[.6,.6,.6]);
L1a = annotation('textbox','position',[.14,.47,.22,.04],'String','Fluc',...
    'fontsize',6,'edgecolor','none');
L1b = annotation('line',[.13,.15],[.487,.487],'linewidth',1,'color',[0,0,0]);

% Then in b_frac...
ax_bii = axes('position',[.11,.54,.22,.14]);
hold on; box on; set(gca,'linewidth',1);
p1 = plot(W3.b_frac{4,9,1},'-','linewidth',1,'color',[.6,.6,.6]);
% plot(W3.b_frac{4,6,2},'-','linewidth',1,'color',[.6,.6,.6]);
% plot(W3.b_frac{4,6,3},'-','linewidth',1,'color',[.6,.6,.6]);
p2 = plot(X.b_frac{4,9,1},'-','linewidth',1,'color',[0,0,0]);
% plot(X.b_frac{4,9,2},'-','linewidth',1,'color',[0,0,0]);
% plot(X.b_frac{4,9,3},'-','linewidth',1,'color',[0,0,0]);

% Labels:
yL = ylabel('Biomass');
set(gca,'fontsize',7);
xlim([0,50]);
ylim([0,1]);
set(gca,'layer','top');
xticks([0,10,20,30,40,50]);
xticklabels({'','','','','',''});
yticks([.25,.5,.75,1]);
set(gca,'ticklength',[.03,.03]);


% PLOT C ------------------------------------------------------------------
% Heatmap of concentration changes in the mean field:
ax_c = axes('position',[.45,.45,.27,.23]); 
hold on; box on; set(gca,'linewidth',1);
imagesc(W3.alpha, 1:length(W3.G_0), mean(W3.FinalC_FS,3));
colormap(ax_c, cmap);
xlim([2,4]);
ylim([1,length(W3.G_0)]);
xticks([2.5,3,3.5]);
yticks([1,6,11]);
yticklabels({'3e3','3e4','3e5'});
xlabel('\lambda');
ylabel('Food Conc. [#/mL]');
clim([-.5,.5]);
set(ax_c,'fontsize',7);
set(ax_c,'layer','top');

% PLOT D ------------------------------------------------------------------
% Heatmap of actual concentration changes with a new state strategy
ax_d = axes('position',[.72,.45,.27,.23]);
hold on; box on; set(gca,'linewidth',1);
imagesc(X.alpha, 1:length(X.G_0), mean(X.FinalC_FS,3));
colormap(ax_d, cmap);
xlim([2,4]);
ylim([1,length(X.G_0)]);
xticks([2.5,3,3.5]);
yticks([1,6,11]);
yticklabels({'','',''});
xlabel('\lambda');
set(ax_c,'ticklength',[.03,.03]);
set(ax_d,'ticklength',[.03,.03]);
% ylabel('Food Concentration [#/mL]');
clim([-.5,.5]);
set(ax_d,'fontsize',7);
set(ax_d,'layer','top');
cb = colorbar('southoutside');
set(cb,'fontsize',6);
set(cb,'linewidth',1);
% set(ax_d,'position',[.38,.07,.27,.46]);
set(cb,'position',[.46,.36,.52,.02]);
set(cb,'ticks',[-.5,-.25,0,.25,.5]);

T1 = annotation('textbox','String','No Fluc');
T1.Position = [.45,.64,.27,.04];
T1.FontSize = 7;
T1.EdgeColor = 'none';
T2 = annotation('textbox','string','Fluc');
T2.Position = [.72,.64,.27,.04];
T2.FontSize = 7;
T2.EdgeColor = 'none';
T3 = annotation('textbox','string','\DeltaB');
T3.Position = [.37,.355,.15,.04];
T3.FontSize = 7;
T3.EdgeColor = 'none';
T3.FontAngle='italic';

T4 = annotation('rectangle');
T4.Position = [.57,.505,.015,.025];
T4.Color = 'k';
T4.LineWidth = 2;
T5 = annotation('rectangle');
T5.Position = [.84,.505,.015,.025];
T5.Color = [0,0,0];
T5.LineWidth = 2;


% PLOT E ------------------------------------------------------------------
% Plot with increasing swimming speeds:
ax_e = axes('position',[.11,.103,.29,.19]);
hold on; box on; set(gca,'linewidth',1);
% plot(V0(:,2), zeros(size(Z.Ds(:,2))), 'k--','linewidth',1);
% errorbar(V0(:,2), mean(Z.FinalC_FS,2), 2*std(Z.FinalC_FS,0,2), 'k-','linewidth',1,'capsize',4);
% xlabel('Speed [\mum/s]');
plot(Z.Ds(:,2), zeros(size(Z.Ds(:,2))), 'k--','linewidth',1);
errorbar(Z.Ds(:,2), mean(Z.FinalC_FS,2), 2*std(Z.FinalC_FS,0,2), 'k-','linewidth',1,'capsize',4);
xlabel('D_s [\mum^2/s]');
ylabel('\DeltaB');
set(gca,'xscale','log');
set(ax_e,'fontsize',7);
set(ax_e,'layer','top');
% xlim([min(V0(:,2)), max(V0(:,2))]);
% xticks([.1,1,10,100]);
xlim([min(Z.Ds(:,2)),max(Z.Ds(:,2))]);
xticks([1e-2,1e-1,1e0,1e1,1e2,1e3])
ylim([-.5,.5]);
set(ax_e,'ticklength',[.04,.04]);
title('Swimming');

% PLOT F ------------------------------------------------------------------
% Plot with differences in growth rates between the states:
ax_f = axes('position',[.40,.103,.29,.19]);
hold on; box on; set(gca,'linewidth',1);
plot(M.MuRatio, zeros(size(M.MuRatio)), 'k--','linewidth',1);
errorbar(M.MuRatio, mean(M.FinalC_FS,2), 2*std(M.FinalC_FS,0,2), 'k-','linewidth',1,'capsize',4);
xlabel('\mu_{mc} / \mu_{sc}');
set(ax_f,'fontsize',7);
set(ax_f,'layer','top');
xlim([min(M.MuRatio), max(M.MuRatio)]);
xticks([.2,.5,.8]);
yticks([-.5,0,.5]);
yticklabels({'','',''});
set(ax_f,'ticklength',[.04,.04]);
title('Growth Penalty');


% PLOT G ------------------------------------------------------------------
ax_g = axes('position',[.69,.103,.29,.19]);
hold on; box on; set(gca,'linewidth',1);
plot(V.Phi, zeros(size(V.Phi)), 'k--','linewidth',1);
errorbar(V.Phi, mean(V.FinalC_FS,2), 2*std(V.FinalC_FS,0,2), 'k-','linewidth',1,'capsize',4);
xlabel('\phi');
% ylabel('\DeltaB');
set(ax_g,'fontsize',7);
set(ax_g,'layer','top');
xlim([0,1]);
xticks([.2,.5,.8]);
ylim([-.5,.5]);
yticks([-.5,0,.5]);
yticklabels({'','',''});
set(ax_g,'ticklength',[.04,.04]);
title('Cell Packing');

% What level of packing fraction changes can counterbalance growth
% penalties?
x1 = V.Phi - V.Phi(end); x2 = mean(V.FinalC_FS,2); x2 = x2 - x2(end);
y1 = M.MuRatio; y2 = mean(M.FinalC_FS,2);
q = -y2;
vq = interp1(x2,x1,q,'linear','extrap');
ix = find(vq >= -1);
Range = y1(ix);
min(Range)

% -------------------------------------------------------------------------
% Panel Labels:
annotation('textbox',[.01,.92,.03,.10],'string','a','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.01,.64,.03,.10],'string','b','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.36,.64,.03,.10],'string','c','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.01,.24,.03,.10],'string','d','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
% annotation('textbox',[.45,.64,.03,.10],'string','e','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
% annotation('textbox',[.45,.64,.03,.10],'string','f','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Figure_05'],'-dpng',res);
    exportgraphics(gcf,'Figure_05.pdf');
    close(gcf);
end


%% Figure 6: Including phage
% LOAD DATA
% First, load from same simulations as above ------------------------------
% datafolder = 'E:\2025\01_ProjectEncounter\202511XX_SimsRevisions\Phage_Food_Sweeps_Sims\';
datafolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\07_Sims\'];
X = load([datafolder,'sims_PHAGE_2026-01-26.mat']);
G0_F = X.G0_F;
G0_P = X.G0_P;
BS = X.Phage_Burst_Size;

% Calculate integrated biomass fraction:
for ff = 1:length(G0_F)
    for pp = 1:length(G0_P)
        for bb = 1:length(BS)
            for mm = 1:10
                IntegratedBiomassFraction(ff,bb,pp,mm) = sum(X.b_frac{bb,pp,mm});
            end
        end
    end
end
IBF = mean(IntegratedBiomassFraction,4);

% Second, load data from Gillespie simulations ----------------------------
gill_folder = [drive_folder,':\Shared drives\Schwartzman Lab\Data\tom\Sweeping_Size_BurstSize_PhageConc_EMP\'];
names = dir(gill_folder);
folder_names = names([names.isdir]);
folder_names = folder_names(3:end);

group_sizes     = [1e1,3e1,1e2,3e2,1e3,3e3,1e4];
phage_conc      = [1e4,3e4,1e5,3e5,1e6,3e6,1e7];
% burst_sizes     = [0,33,67,100,133,167,200];
% burst_sizes     = [0,50,100,150,200];
burst_sizes     = [0,100,200];
Msims           = [1,1,1,3,10,30,100];

% Calculate some important quantities:
T_avg_2  = zeros(length(folder_names),1);
T_err_2  = T_avg_2;
T_calc_2 = T_avg_2;
Infection_Array = cell(length(folder_names),1);
Infection_Array_Ordered = Infection_Array;

% Loop through all parameter combinations:
for nn = 1:length(folder_names)
    disp(folder_names(nn).name);
    filelist = dir([gill_folder,folder_names(nn).name,'/*.mat']);
    msims_this_size = length(filelist);

    % Set up the measurement arrays:
    T_inf  = zeros(length(filelist),1);
    T_calc = T_inf;
    T_err  = T_inf;

    % Loop through the files in the folder:
    for ff = 1:length(filelist)
        disp(filelist(ff).name);
        Z = load([gill_folder,folder_names(nn).name,'\',filelist(ff).name]);
        
        % Record time of infection for all members in the file:
        ix_2 = find(Z.LOG(end).STATE == 2); % get all the MC states
        BD = Z.LOG(end).BANK_DEATH;
        Infection_Times = BD(ix_2);

        % Order the infection times into an array:
        Infection_Array_Ordered{nn} = [Infection_Array_Ordered{nn}; sort(Infection_Times)];

        % Get the average infection time, and also append the whole list to
        % an array:
        T_inf(ff) = mean(Infection_Times);
        if length(Infection_Times) > 1
            T_err(ff) = std(Infection_Times);
        end
        Infection_Array{nn} = [Infection_Array{nn}; Infection_Times];

        % Calculate the expected infection time:
        R = computePhageAcquisitionRate(Z.initial_biomass, Z.phage_conc);
        T_calc(ff) = 1./R;
        
    end

    % Get averages from across the different sample runs:
    T_avg_2(nn) = mean(T_inf,1);
    T_err_2(nn) = mean(T_err,1);
    T_calc_2(nn) = mean(T_calc);
end

% Reshape all arrays to match the parameter sweeps:
name_list = cell(length(folder_names),1);
for ii = 1:length(folder_names)
    name_list{ii} = folder_names(ii).name;
end
name_list = reshape(name_list, length(phage_conc), length(burst_sizes), length(group_sizes));
Infection_Array = reshape(Infection_Array, length(phage_conc), length(burst_sizes), length(group_sizes));
Infection_Array_Ordered = reshape(Infection_Array_Ordered, length(phage_conc), length(burst_sizes), length(group_sizes));

% Record an infection matrix:
IM_bar = zeros(size(Infection_Array));
IM_err = IM_bar;
for nn = 1:length(group_sizes)
    for bb = 1:length(burst_sizes)
        for pp = 1:length(phage_conc)

            % Remove infs from dataset:
            % finite_data = Infection_Array{pp,bb,nn}(~isinf(Infection_Array{pp,bb,nn}) & ~isnan(Infection_Array{pp,bb,nn}));

            % Make Infs the max simulation time before averaging:
            data = Infection_Array{pp,bb,nn};
            data(data == Inf) = Z.tmax;
            finite_data = data;

            % Average:
            IM_bar(pp,bb,nn) = mean(finite_data);
            IM_err(pp,bb,nn) = std(finite_data);
        end
    end
end
% avg_correction_factor = mean(IM_bar(:,4)./IM_bar(:,1));

% -------------------------------------------------------------------------
%% FIGURE
figure('units','centimeters','position',[3,3,big_width,10]);

% A & B: ------------------------------------------------------------------

% Sweep burst size holding phage and food constant:
rowa_x = [.10,.32,.54,.76];
rowa_y = 0.78;
rowa_w = 0.22;
rowa_h = 0.20;
% ff = 3;
pp = 2;
for bb = 1:length(BS)
    ax = axes('position',[rowa_x(bb), rowa_y, rowa_w, rowa_h]);
    hold on; box on; set(gca,'linewidth',1);
    for mm = 1:10
        plot(X.b_frac{bb,pp,mm},'-','color',[.8,.2,.2]);
        plot(X.n_frac{bb,pp,mm},'-','color','k');
    end
    xlim([0,218]);
    ylim([0,1]);
    text(150,.90,['BS=',num2str(BS(bb))],'fontsize',6);
    xticks([50,100,150,200]);
    yticks([.25,.5,.75,1]);
    if bb == 1
        % ylabel('Biomass Fraction');
        text(50,.75,'Biomass','color',[.8,.2,.2],'fontsize',7);
        text(50,.10,'Number','color','k','fontsize',7);
    else
        yticklabels({'','','','',''});
    end
    set(gca,'fontsize',7);
    set(gca,'ticklength',[.02,.02]);
    xlabel('Time (# Rounds)');
end

% Sweep phage holding food and burst size constant:
rowb_y = .58;
rowb_w = .88/4;
rowb_h = .20;
rowb_x = rowb_w*(0:6) + 0.10;

bb = 3;
for pp = 1:4%length(G0_P)
    ax = axes('position',[rowb_x(pp), rowb_y, rowb_w, rowb_h]);
    hold on; box on; set(gca,'linewidth',1);
    for mm = 1:10
        plot(X.b_frac{bb,pp,mm},'-','color',[.8,.2,.2]);
        plot(X.n_frac{bb,pp,mm},'k-');
    end
    xlim([0,218]);
    ylim([0,1]);
    xticks([50,100,150,200]);
    yticks([.25,.5,.75,1]);
    text(155,.9,['P=',num2str(G0_P(pp),'%1.0e')],'fontsize',6);
    if pp == 1
        % ylabel('Biomass Fraction');
    else
        yticklabels({'','','','',''});
    end
    set(gca,'fontsize',7);
    set(gca,'ticklength',[.02,.02]);
    xlabel('Time (# Rounds)');
end

% Panel C: ----------------------------------------------------------------
%{
% Left:
axC1 = axes('position',[.07,.43,.17,.15]);
hold on; box on; set(gca,'linewidth',1);
imagesc(IBF(:,:,3));
xlabel('Phage Conc.');
ylabel('Food Conc.');
xticks([]);
yticks([]);
xlim([.5,7.5]);
ylim([.5,7.5]);
clim([0,125]);
set(gca,'layer','top');
colormap(axC1,'pink');
set(gca,'fontsize',7);

% Right:
axC2 = axes('position',[.28,.43,.17,.15]);
hold on; box on; set(gca,'linewidth',1);
imagesc(squeeze(IBF(3,:,:)));
xticks([1,2,3,4]);
xticklabels({'0','50','100','200'});
yticks([]);
xlabel('Burst Size');
ylabel('Phage Conc.');
xlim([.5,4.5]);
ylim([.5,7.5]);
clim([0,125]);
cb = colorbar('eastoutside');
set(gca,'layer','top');
cb.Position = [.46,.43,.015,.15];
cb.LineWidth = 1;
colormap(axC2,'pink');
set(gca,'fontsize',7);
%}

% PANEL C: Timescales -----------------------------------------------------
%{
axC = axes('position',[.05,.05,.43,.50]);
hold on; box on; set(gca,'linewidth',1);
set(gca,'fontsize',7);

% Sweeping two variables:
n_cells = logspace(0,4,100); % number of cells in a group
n_p = logspace(4,7,50); % phage concentration
radius_phage = 1e-5; % cm
correction_factor = 1; %avg_correction_factor; %0.27; % this comes from the above analysis

% Two models: Diffusion only, or turbulence only
for nn = 1:length(n_cells)
    for pp = 1:length(n_p)
        Gamma_D(nn,pp) = computePhageAcquisitionRate(n_cells(nn), n_p(pp));
        Gamma_T(nn,pp) = computeFoodAcquisitionRate(n_cells(nn), radius_phage, n_p(pp));
    end
end

% Consider the timescale to get X cells from a single cell - is this possible?
Colormap = [0,60,48; 
            1,102,94;
            53,151,143;
            128,205,193;
            199,234,229;
            245,245,245;
            246,232,195;
            223,194,125;
            191,129,45;
            140,81,10;
            84,48,5]/255;
doubling_time = 30*60; % secs
assembly_time = log2(n_cells) * doubling_time;
Assembly_Time = repmat(assembly_time', [1,length(n_p)]);
A = correction_factor*(1./Gamma_D);
PossibilityMatrix = A./3600;%./Assembly_Time;

hold on; box on; set(gca,'linewidth',1);
imagesc(n_p, n_cells, PossibilityMatrix);
colormap(axC, Colormap);
% [C,h] = contour(n_p, n_cells, PossibilityMatrix,[0.1,1,10],'w-','linewidth',1.5,'ShowText',true,'LabelFormat','%1.1f DT');
[C,h] = contour(n_p, n_cells, PossibilityMatrix,[20,60,180,600]/60,'w-','linewidth',1.5,'ShowText',true,'LabelFormat','%2.2f [hrs]');
clabel(C,h,'fontsize',8,'color','w');
text(2e4,3e0,'Can assemble','fontsize',8,'color','w');
text(1.5e6, 2e3,{'Cannot','assemble'},'fontsize',8,'color','w');
% clim([0,10]);
clim([0,10]);
xlim([1e4,1e7]);
ylim([1,1e4]);
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'layer','top');
set(gca,'fontsize',7);
xlabel('Infect. Phage Conc [#/mL]');
ylabel('Group Biomass [# Cells]');
%}

% Panel C, Timescales with turbulence: ------------------------------------
axC = axes('position',[.05,.08,.24,.40]);
hold on; box on; set(gca,'linewidth',1);
set(gca,'fontsize',7);

% Sweeping two variables:
n_cells = logspace(0,4,100); % number of cells in a group
n_p = logspace(4,7,50); % phage concentration
radius_phage = 1e-5; % cm
correction_factor = 1; %avg_correction_factor; %0.27; % this comes from the above analysis

% Two models: Diffusion only, or turbulence only
for nn = 1:length(n_cells)
    for pp = 1:length(n_p)
        Gamma_D(nn,pp) = computePhageAcquisitionRate(n_cells(nn), n_p(pp));
        Gamma_T(nn,pp) = computeFoodAcquisitionRate(n_cells(nn), radius_phage, n_p(pp));
    end
end
% Gamma_N = max(cat(3,Gamma_D, Gamma_T),[],3);
Gamma_N = Gamma_T;

% Consider the timescale to get X cells from a single cell - is this possible?
Colormap = [0,60,48; 
            1,102,94;
            53,151,143;
            128,205,193;
            199,234,229;
            245,245,245;
            246,232,195;
            223,194,125;
            191,129,45;
            140,81,10;
            84,48,5]/255;
doubling_time = 30*60; % secs
assembly_time = (log2(n_cells) * doubling_time)/3600; % in hrs
Assembly_Time = repmat(assembly_time', [1,length(n_p)]);
A = correction_factor*(1./Gamma_N);
PossibilityMatrix = A./3600;%./Assembly_Time;
AssemblyMatrix = PossibilityMatrix./Assembly_Time;

hold on; box on; set(gca,'linewidth',1);
imagesc(n_p, n_cells, PossibilityMatrix);
colormap(axC, Colormap);
% [C,h] = contour(n_p, n_cells, PossibilityMatrix,[0.1,1,10],'w-','linewidth',1.5,'ShowText',true,'LabelFormat','%1.1f DT');
[C,h] = contour(n_p, n_cells, PossibilityMatrix,[20,60,180,720,1440]/60,'w-','linewidth',1.5,'ShowText',true,'LabelFormat','%2.2f [hrs]');
clabel(C,h,'fontsize',8,'color','w');
% text(2e4,3e0,'Can assemble','fontsize',8,'color','w');
% text(1.5e6, 2e3,{'Cannot','assemble'},'fontsize',8,'color','w');
% clim([0,10]);
clim([1,12]);
xlim([1e4,1e7]);
ylim([1,1e4]);
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'layer','top');
set(gca,'fontsize',7);
xlabel('Infect. Phage Conc [#/mL]');
ylabel('Group Biomass [# Cells]');

% % Inset: Assembly time 100 cells ----------
% axCi = axes('position',[.32,.38,.15,.15]);
% hold on; box on; set(gca,'linewidth',1);
% imagesc(n_p, n_cells, PossibilityMatrix > 2.5);
% colormap(axCi, Colormap);
% % [C,h] = contour(n_p, n_cells, PossibilityMatrix,[0.1,1,10],'w-','linewidth',1.5,'ShowText',true,'LabelFormat','%1.1f DT');
% % [C,h] = contour(n_p, n_cells, AssemblyMatrix,linspace(.01,30,100),'k-','linewidth',1.5);%,'ShowText',true,'LabelFormat','%2.2f');
% % clabel(C,h,'fontsize',8,'color','w');
% clim([0,6]);
% xlim([1e4,1e7]);
% ylim([1,1e4]);
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% set(gca,'layer','top');
% set(gca,'fontsize',7);
% % xlabel('Infect. Phage Conc [#/mL]');
% % ylabel('Group Biomass [# Cells]');

% PANEL D: Gillespie simulations ------------------------------------------
rowd_h = .40;
rowd_w = 0.21;
% rowd_x = [.56,.77,.56,.77];
% rowd_y = [.31,.31,.07,.07];
rowd_y = [.078,.08,.08];
rowd_x = [.35,.56,.77];

Colors_BS = abyss(7);
conc_list = [1,4,7];
% burst_list = [1,3,5];
burst_list = [1,2,3];
for pp = 1:3%length(conc_list)
    ax = axes('position',[rowd_x(pp), rowd_y(pp), rowd_w, rowd_h]);
    hold on; box on; set(gca,'linewidth',1);
    for bb = 1:length(burst_list)
        % burst_sizes(burst_list(bb))
        errorbar(group_sizes, squeeze(IM_bar(conc_list(pp),burst_list(bb),:)./3600), squeeze(IM_err(conc_list(pp),burst_list(bb),:)./3600), 'o-','linewidth',1,'color',Colors_BS(burst_list(bb),:));
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([5,1.7e4]);
    set(gca,'fontsize',7);
    xticks([1e1,1e2,1e3,1e4]);
    xlabel('Group biomass [# Cells]');
    ylim([.01,120]);
    yticks([1e-1,1e0,1e1,1e2]);
    if ismember(pp,[1])
        ylabel('\Delta t [hrs]');
    else
        yticklabels({'','','',''});
    end
    set(gca,'ticklength',[.03,.03]);
    set(gca,'layer','top');
    text(1e3, 70, ['\phi=',num2str(phage_conc(conc_list(pp)), '%1.0e')],'fontsize',7)
end
legend({'BS=0','100','200'},'location','northwest')

% Annotations:
% Box1a = annotation('rectangle',[.10,.84,.88,.14],'color',[.2,.4,.8],'linewidth',2);
% Box1b = annotation('rectangle',[.28,.495,.17,.02],'color',[.2,.4,.8],'linewidth',2);
% Box2a = annotation('rectangle',[.10,.65,.88,.14],'color',[.7,.3,.7],'linewidth',2);
% Box2b = annotation('rectangle',[.07,.475,.17,.02],'color',[.7,.3,.7],'linewidth',2);
an = annotation('textbox', [.07,.71,.15,.05], 'String', 'S2 Fraction', 'FitBoxToText', 'on','fontsize',7,'edgecolor','none');
an.Rotation = 90; % Rotates 90 degrees counter-clockwise
annotation('textbox',[.01,.92,.03,.10],'string','a','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.01,.72,.03,.10],'string','b','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.00,.44,.03,.10],'string','c','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.30,.44,.03,.10],'string','d','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Figure_06'],'-dpng',res);
    exportgraphics(gcf,'Figure_06.pdf');
    close(gcf);
end



%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% SUPPLEMENTAL FIGURES
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% Adding kernels together in quadrature rather than simple superposition:
% Inspired by Saffman-Turner 1956.

% 3D Encounter kernel (theory): -------------------------------------------
kT    = 1.38e-23 * 300;          % thermal energy in [J]
nu    = 1e-3;                    % water dynamic viscosity in [Pa*s]
mu    = 1e-6;                    % water kinematic viscosity (nu/rho)
DRhor = 25;                      % difference in density between bacteria and water [kg/m^3]
g     = 10;                      % acceleration due to gravity [m/s^2]
eps   = 1e-5;                    % ocean energy dissipation rate, [W/kg]
rb    = logspace(-10,-3,1e3);    % bacteria agg radius in [m]
rr    = logspace(-10,-3,1e3);    % resource patch radius in [m]
Conv2UM = (1e6)^3;               % Conversion factor from cubic meters to cubic microns
Conv2CM = (1e2)^3;               % Conversion factor from cubic meters to cubic centimeters (mL)

% Make patch size matrix:
[Rb, Rr] = meshgrid(rb, rr);

% Diffusive Kernel:
Db = kT./(6*pi*nu*Rb);
Dr = kT./(6*pi*nu*Rr);
G_D = 4*pi*(Db + Dr) .* (Rb + Rr);

% Buoyancy Kernel:
Ub = 2 * DRhor * g * Rb.^2 / (9*nu);
Ur = 2 * DRhor * g * Rr.^2 / (9*nu);
G_B = pi * (Rb + Rr).^2 .* abs(Ub-Ur);

% Turbulence Kernel, chosen for a particular epsilon value, and for well below the Kolmogorov limit:
G_T = 1.3 * (Rb + Rr).^3 .* sqrt(eps./mu);

% Divide each kernel by (Rb+Rr).^2:
G_D_vel = G_D./((Rb+Rr).^2);
G_B_vel = G_B./((Rb+Rr).^2);
G_T_vel = G_T./((Rb+Rr).^2);

% Add in quadrature:
G_N_vel = sqrt(G_D_vel.^2 + G_B_vel.^2 + G_T_vel.^2);

% Re-multiply (Rb+Rr)^2
G_N = (Rb+Rr).^2 .* G_N_vel;
G_N_orig = G_D + G_B + G_T;

% Diffusion:
%{
figure;
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','log');
set(gca,'yscale','log');
imagesc(1e6*rb, 1e6*rr, log10(G_D));
colormap(turbo(25));

% Diffusion vel:
figure;
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','log');
set(gca,'yscale','log');
imagesc(1e6*rb, 1e6*rr, log10(G_D_vel));
colormap(turbo(25));

% Bouyancy:
figure;
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','log');
set(gca,'yscale','log');
imagesc(1e6*rb, 1e6*rr, log10(G_B));
colormap(turbo(25));

% Bouyancy vel:
figure;
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','log');
set(gca,'yscale','log');
imagesc(1e6*rb, 1e6*rr, log10(G_B_vel));
colormap(turbo(25));

% Turb:
figure;
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','log');
set(gca,'yscale','log');
imagesc(1e6*rb, 1e6*rr, log10(G_T));
colormap(turbo(25));

% Turb vel:
figure;
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','log');
set(gca,'yscale','log');
imagesc(1e6*rb, 1e6*rr, log10(G_T_vel));
colormap(turbo(25));
%}

% Net:
figure('units','centimeters','position',[3,3,10,8]);
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','log');
set(gca,'yscale','log');
imagesc(1e6*rb, 1e6*rr, log10(Conv2CM * G_N));
colormap(turbo(25));
colorbar;
xlim([1e-1,1e3])
ylim([1e-4,1e3]);
xlabel('Microbe Size [\mum]');
ylabel('Resource Size [\mum]');
set(gca,'layer','top');
clim([-11, 0]);
set(gca,'fontsize',7);
% figure;
% hold on; box on; set(gca,'linewidth',1);
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% imagesc(1e6*rb, 1e6*rr, log10(G_N_orig));
% colormap(turbo(25));

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_Quadrature'],'-dpng',res);
    close(gcf);
end

%% Changing the densities of sinking objects
% 3D Encounter kernel (theory): -------------------------------------------
kT    = 1.38e-23 * 300;          % thermal energy in [J]
nu    = 1e-3;                    % water dynamic viscosity in [Pa*s]
mu    = 1e-6;                    % water kinematic viscosity (nu/rho)
DRhob = linspace(0,200,201);                      % difference in density between particles and water [kg/m^3]
DRhor = linspace(0,200,201);                      % difference in density between patches and water [kg/m^3]
% DRhor = 25;
g     = 10;                      % acceleration due to gravity [m/s^2]
eps   = 1e-5;                    % ocean energy dissipation rate, [W/kg]
rb    = 1e-5;                    % bacteria agg radius in [m]
rr    = 1e-6;                    % resource patch radius in [m]
Conv2UM = (1e6)^3;               % Conversion factor from cubic meters to cubic microns
Conv2CM = (1e2)^3;               % Conversion factor from cubic meters to cubic centimeters (mL)

% Make density matrix:
[rho_b, rho_r] = meshgrid(DRhob, DRhor);

% Diffusive Kernel:
Db = kT./(6*pi*nu*rb);
Dr = kT./(6*pi*nu*rr);
G_D = 4*pi*(Db + Dr) .* (rb + rr);

% Buoyancy Kernel:
Ub = 2 * rho_b * g * rb.^2 / (9*nu);
Ur = 2 * rho_r * g * rr.^2 / (9*nu);
G_B = pi * (rb + rr).^2 .* abs(Ub-Ur);

% Turbulence Kernel, chosen for a particular epsilon value, and for well below the Kolmogorov limit:
G_T = 1.3 * (rb + rr).^3 .* sqrt(eps./mu);
% G_T2 = 1.37*pi.*eps ^ (1/3) * (Rb+Rr).^(7/3);
% Lk = 0.5*(mu^3./eps).^(1/4);
% ix1 = Rb+Rr <= Lk;
% ix2 = Rb+Rr > Lk;
% G_T = G_T1.*ix1 + G_T2.*ix2;

% Make a swimming kernel, too: only one item swims
% G_S = 4*pi*SwimStrength .* (Rb + Rr);

% Net kernel:
% G_N = Conv2UM * (G_D + G_B + G_T);
G_N = Conv2CM * (G_D + G_B + G_T);

% Show figure:
figure('units','centimeters','position',[3,3,10,10]); 
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','linear');
set(gca,'yscale','linear');
imagesc(DRhob, DRhor, log10(G_N));
% plot(DRhob, G_N,'-','linewidth',2,'color','k');
axis equal;
xlim([0,200]);
ylim([0,200]);
set(gca,'layer','top');
xlabel('\Delta\rho Bacteria [kg/m^3]');
ylabel('\Delta\rho Resource [kg/m^3]');


% Colors:
MinRange = min(log10(G_N),[],'all');
MaxRange = max(log10(G_N),[],'all');
Ncolors  = 10;
CMap = bone(Ncolors);
clim([MinRange, MaxRange]);
colormap(CMap);
cb = colorbar; %('AxisLocationMode','manual','position',[.48,.48,.015,.50]);
set(cb,'linewidth',1);
cb.Label.String = 'log10(Encounter Kernel [cm^3/s])';

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_VaryingDensity'],'-dpng',res);
    close(gcf);
end

%% Supplemental Timetracks:
% Show the average number of particles attached to aggregates of a specific
% size class, over time.

clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

% Inputs:
FilterFLPerArea = 5e3;
MinSize = 5;

% Load data:
datafolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\02_Timelapse\A\'];
datalist = dir([datafolder,'*.mat']);
Time = [0,30,60,90,120,245];
n_0 = 1e7;
X = LOAD_DATA_BIN_CALIBRATE_TIMETRACKS([datafolder, datalist(3).name], FilterFLPerArea, MinSize, 1);

% Track a few bins:
% ix_track = 1:8;
ix_track = 8:28;
SizeClass = X.Avg_Rad(ix_track,1);
SizeClassErr = X.Avg_Rad(ix_track,2);
SizeClassLabels = cell(1,1);
for bb = 1:length(SizeClass)
    SizeClassLabels{bb} = [num2str(SizeClass(bb),2),' \mum'];
end

% Within tracked bins, get the average number of particles attached.
P = zeros(length(ix_track), length(Time));
S = zeros(length(ix_track), length(Time));
for tt = 2:length(Time)
    X = LOAD_DATA_BIN_CALIBRATE_TIMETRACKS([datafolder,datalist(tt).name], FilterFLPerArea, MinSize, 1);
    P(:,tt) = X.Avg_Rad(ix_track, 3);
    S(:,tt) = X.Avg_Rad(ix_track, 4)./sqrt(X.N_Rad(ix_track));
end

% Set up the theoretical result:
modelfun_allvary   = @(parameters,t) ( parameters(1)/parameters(2) * 1e7 * (1 - exp(-parameters(2)*t)) );

% Do a regression given the theory result:
beta_0  = 1e-2;
gamma_0 = 1e-8*linspace(1,20,length(ix_track));

% Fit the model, letting all terms vary:
fit_params = zeros(2, length(ix_track));
err_params = fit_params;
for ii = 1:length(ix_track)
    parameters0 = [gamma_0(ii), beta_0];
    mdl = fitnlm(Time, P(ii,:), modelfun_allvary, parameters0);
    fit_params(:,ii) = mdl.Coefficients.Estimate;
    err_params(:,ii) = mdl.Coefficients.SE;
    tst_params(:,ii) = mdl.Coefficients.tStat;
    pval_params(:,ii) = mdl.Coefficients.pValue;
end

% Get average values for beta:
beta_bar  = mean(fit_params(2,:));

% Define functions where beta is not allowed to vary:
modelfun_gammavary = @(parameters,t) ( parameters(1)/.025 * 1e7 * (1 - exp(-.019*t)) );
modelfun_betavary  = @(parameters,t) ( 2.4e-8/parameters(1) * 1e7 * (1 - exp(-parameters(1)*t)) );

% Fit the model with restrictions on beta:
for ii = 1:length(ix_track)

    % Beta restricted:
    parameters0 = [gamma_0(ii)];
    mdl_gammavary = fitnlm(Time, P(ii,:), modelfun_gammavary, parameters0);
    fit_params_gammavary(:,ii) = mdl_gammavary.Coefficients.Estimate;

    % Gamma restricted:
    parameters0 = 0.01;
    mdl_betavary = fitnlm(Time, P(ii,:), modelfun_betavary, parameters0);
    fit_params_betavary(:,ii) = mdl_betavary.Coefficients.Estimate;

end

% Generate residuals for each restriction:
residuals = zeros(3,length(ix_track)); % Pre-allocate space for residuals
for ii = 1:length(ix_track)
    % All vary:
    y_reg = modelfun_allvary(fit_params(:,ii), Time);
    residuals(1,ii) = sum(sqrt((P(ii,2:end) - y_reg(2:end)).^2)./y_reg(2:end));

    % Beta restriction:
    y_reg = modelfun_gammavary(fit_params_gammavary(:,ii), Time);
    residuals(2,ii) = sum(sqrt((P(ii,2:end) - y_reg(2:end)).^2)./y_reg(2:end));

    % Gamma restriction:
    y_reg = modelfun_betavary(fit_params_betavary(:,ii), Time);
    residuals(3,ii) = sum(sqrt((P(ii,2:end) - y_reg(2:end)).^2)./y_reg(2:end));
end

% Statistics of the fit parameters: ---------------------------------------
% Calculate a Pearson's r-correlation for trends of Gamma and Beta:
c_gamma = cov(SizeClass, fit_params(1,:)');
rho_gamma = c_gamma(1,2)./sqrt(c_gamma(1,1).*c_gamma(2,2));
c_beta  = cov(SizeClass, fit_params(2,:)');
rho_beta = c_beta(1,2)./sqrt(c_beta(1,1)*c_beta(2,2));

% Do a Mann-Kendall test for trends:
[H_gamma, p_gamma] = Mann_Kendall(fit_params(1,:)', 0.05);
[H_beta, p_beta] = Mann_Kendall(fit_params(2,:)', 0.05);
fprintf(['Mann-Kendall test says the p-values are: \n',num2str(p_gamma),' & ',num2str(p_beta),'\n'])

% Make a table of the beta fit parameters and their t-scores and p-values:
beta = fit_params(2,:)';
se = err_params(2,:)';
t_score = tst_params(2,:)';
p_val = pval_params(2,:)';
T1 = table(SizeClass,beta,se,t_score,p_val);
writetable(T1,[printfolder,'Supplemental_BetaFitParams.csv']);

% Let's take the average of beta from 10um up, and compare all other beta
% fits to those ones:
n = 6; % number of datapoints in each sample
ix = find(SizeClass > 10);
avg_beta = mean(beta(ix));
delta_beta = beta - avg_beta;
t_score_avg = delta_beta./se;
p_val_avg = 2*(1-tcdf(abs(t_score_avg),n-1));
T2 = table(SizeClass, beta, se, t_score_avg, p_val_avg);

% Show figures: --------------------------------------------------------
figure('units','centimeters','position',[3,3,15,15]);

% All parameters vary:
ax = axes('position',[.05,.55,.40,.44]);
hold on; box on; set(gca,'linewidth',1);
Colors = turbo(length(ix_track));
for ii = 1:length(ix_track)
    y_reg = modelfun_allvary(fit_params(:,ii), Time);
    t_plot = linspace(0,max(Time),1e3);
    y_plot = modelfun_allvary(fit_params(:,ii), t_plot);
    residuals(1,ii) = sum(sqrt((P(ii,2:end) - y_reg(2:end)).^2)./y_reg(2:end));
    errorbar(Time, P(ii,:), S(ii,:), '.-','markersize',12,'linewidth',1,'color',Colors(ii,:),'capsize',6);
    plot(t_plot, y_plot, '-','linewidth',1,'color',Colors(ii,:));
end
set(gca,'layer','top');
xlabel('Time [min]');
ylabel('No. Microbeads');
set(gca,'fontsize',7);
xlim([0,250]);
% ylim([0,30]);
set(gca,'yscale','linear');
set(gca,'xscale','linear');
colormap(Colors);

% Restrict beta:
ax = axes('position',[.52,.55,.40,.44]);
hold on; box on; set(gca,'linewidth',1);
Colors = turbo(length(ix_track));
for ii = 1:length(ix_track)
    y_reg = modelfun_gammavary(fit_params_gammavary(:,ii), Time);
    t_plot = linspace(0,max(Time),1e3);
    y_plot = modelfun_gammavary(fit_params_gammavary(:,ii), t_plot);
    residuals(2,ii) = sum(sqrt((P(ii,2:end) - y_reg(2:end)).^2)./y_reg(2:end));
    errorbar(Time, P(ii,:), S(ii,:), '.-','markersize',12,'linewidth',1,'color',Colors(ii,:),'capsize',6);
    plot(t_plot, y_plot, '-','linewidth',1,'color',Colors(ii,:));
end
set(gca,'layer','top');
xlabel('Time [min]');
ylabel('No. Microbeads');
set(gca,'fontsize',7);
xlim([0,250]);
% ylim([0,30]);
set(gca,'yscale','linear');
set(gca,'xscale','linear');
colormap(Colors);

% Colorbar:
cb = colorbar('location','manual');
cb.Position = [.93,.55,.015,.44];
cb.LineWidth = 1;
clim([1,length(ix_track)]);
cb.Ticks = 1:2:length(ix_track);
cb.TickLabels = SizeClassLabels(1:2:end);

% Show the fit table over size:
ax = axes('position',[.05,.05,.40,.40]);
hold on; box on; set(gca,'linewidth',1);
yyaxis left;
errorbar(SizeClass, fit_params(1,:), 2*err_params(1,:),'.-','linewidth',1,'markersize',8);
ylabel('\Gamma [mL/min]');
ylim([0,6e-8]);
yyaxis right;
errorbar(SizeClass, fit_params(2,:), 2*err_params(2,:),'.-','linewidth',1,'markersize',8);
plot(SizeClass, avg_beta*ones(size(SizeClass)), '--','linewidth',1);
ylabel('\beta [1/min]');
ylim([0,0.07]);
xlabel('Radius [\mum]');
set(gca,'fontsize',7);

% Show the residuals:
ax = axes('position',[.57,.05,.40,.40]);
hold on; box on; set(gca,'linewidth',1);
plot(SizeClass, residuals(1,:),'o-','linewidth',1,'color','k');
plot(SizeClass, residuals(2,:), 's-','linewidth',1,'color',[.6,.6,.6]);
plot(SizeClass, residuals(3,:), '^-','linewidth',1,'color',[.8,.4,.08]);
xlabel('Radius [\mum]');
ylabel('Residuals');
set(gca,'yscale','linear');
legend({'All','\beta rest.','\Gamma rest.'},'location','northeast'); legend('boxoff');
set(gca,'fontsize',7);

annotation('textbox',[.065,.94,.02,.05],'string','a','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.54,.94,.02,.05],'string','b','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.07,.40,.02,.05],'string','c','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')
annotation('textbox',[.58,.40,.02,.05],'string','d','fontsize',9,'fontweight','bold','horizontalalignment','center','verticalalignment','middle','edgecolor','none')

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_Timetracks'],'-dpng',res);
    close(gcf);
end

% % Restrict gamma:
% figure;
% % ax = axes('position',[.52,.55,.40,.44]);
% hold on; box on; set(gca,'linewidth',1);
% Colors = turbo(length(ix_track));
% for ii = 1:length(ix_track)
%     y_reg = modelfun_betavary(fit_params_betavary(:,ii), Time);
%     t_plot = linspace(0,max(Time),1e3);
%     y_plot = modelfun_betavary(fit_params_betavary(:,ii), t_plot);
%     residuals(3,ii) = sum(sqrt((P(ii,2:end) - y_reg(2:end)).^2)./y_reg(2:end));
%     errorbar(Time, P(ii,:), S(ii,:), '.-','markersize',12,'linewidth',1,'color',Colors(ii,:),'capsize',6);
%     plot(t_plot, y_plot, '-','linewidth',1,'color',Colors(ii,:));
% end
% set(gca,'layer','top');
% xlabel('Time [min]');
% ylabel('No. Particles');
% set(gca,'fontsize',7);
% xlim([0,250]);
% % ylim([0,30]);
% set(gca,'yscale','linear');
% set(gca,'xscale','linear');
% colormap(Colors);




%% Show linear version of the regression:
% Inputs:
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

FilterFLPerArea = 5e3;
MinSize         = 4;
MaxSize         = 50;
FilterSize      = [MinSize, MaxSize];
MicrobeadDiam   = 1;
fig_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\'];

% Loading and analysis ----------------------------------------------------
% Load a single example:
X = LOAD_DATA_BIN_CALIBRATE([fig_folder, '02_Timelapse\A\', 'A_090mins_001_EDF.nd2_TD_c_Probabilities.mat'], FilterFLPerArea, FilterSize, MicrobeadDiam);


% FIGURE ------------------------------------------------------------------
figure('units','centimeters','position',[3,3,8,8]);
hold on; box on; set(gca,'linewidth',1);
plot(X.AggregateRadiusFilt, X.NumParticlesFilt, 'o','markersize',3,'markerfacecolor',[.65,.65,.65],'markeredgecolor','none');
errorbar(X.Avg_Particles(:,1), X.Bins_Particles, 0.6, 0.6, X.Avg_Particles(:,2), X.Avg_Particles(:,2), 'd','markersize',6,'markerfacecolor','k','markeredgecolor','none','linewidth',1,'color','k');
p3=plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',1,'color',[.8,.2,.2]);
set(gca,'xscale','linear');
set(gca,'yscale','linear');
ylim([0, 80]);
xlim([3,30]);
set(gca,'fontsize',7);
set(gca,'TickLength',[.02,.02]);
set(gca,'layer','top');
xlabel('Group Radius [\mum]');
ylabel('# Microbeads');

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_Regression_Linear'],'-dpng',res);
    close(gcf);
end

%% Show timelapse of microbeads with regressions:
% Do regressions of one treatment over time. Track the way that the
% particle attachments scale.
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

datafolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\02_Timelapse\A\'];

% Get data:
datalist   = dir([datafolder,'*.mat']);
chosen_dataset = 1:6;

% Inputs:
FilterFLPerArea = 5e3;
FilterSize = [4, Inf];
MicroBeadDiam = 1;
Time = [0,30,60,90,120,245];

% Set up figure:
figure('units','centimeters','position',[3,3,18,4]);
axs = tight_subplot(1,length(chosen_dataset), 0, 0.17, 0.05);

for tt = 1:length(chosen_dataset)

    % Load data, bin it, and do regressions:
    datalist(chosen_dataset(tt)).name
    X = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(chosen_dataset(tt)).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
    
    % Make figure:
    axes(axs(tt));
    hold on; box on; set(gca,'linewidth',1);
    
    % Plot raw data:
    plot(X.AggregateRadiusFilt, X.NumParticlesFilt, '.','markersize',8,'color',[.7,.7,.7]);

    % Plot data without outliers:
    for bb = 1:length(X.Group_Particles_Out)
        plot(X.Group_Particles_Out{bb}(:,1), X.Group_Particles_Out{bb}(:,2), '.','markersize',8,'color',[.4,.4,.4]);
    end
    
    % Plot bin averaging:
    errorbar(X.Avg_Particles(:,1), X.Bins_Particles,...
        0.6*ones(size(X.Bins_Particles)),  0.6*ones(size(X.Bins_Particles)), ...
        X.Avg_Particles(:,2), X.Avg_Particles(:,2), ...
        'd','markerfacecolor','k','markeredgecolor','none','linewidth',1,'color','k');
    
    % Plot regressions:
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',2,'color',[.8,.2,.2]);
    % plot(X.X_Rad(:,2), X.R_Rad, '-','linewidth',2,'color',[.2,.4,.8]);
    % plot(X.X_ParticlesW(:,2), X.R_ParticlesW, '-','linewidth',2,'color',[.2,.7,.2]);
    
    % Show slope:
    text(6, 300, ['\lambda = ',num2str(X.B_Particles(2),'%0.2f'),' \pm ',num2str(2*X.m_std_error,'%0.2f')],'fontsize',8);

    % Aesthetics:
    xlim([5,50]);
    ylim([5e-1,5e2]);
    xticks([10,20,30,40]);
    yticks([1,3,10,30,100]);
    xticklabels({'10','20','30','40'});
    if tt == 1
        yticklabels({'1','3','10','30','100'});
        ylabel('No. microbeads');
    end
    set(gca,'ticklength',[.03,.03]);
    xlabel('Radius [\mum]');
    % title([num2str(X.B_Particles(2)),'; ', num2str(X.B_Rad(2)), '; ', num2str(X.B_ParticlesW(2))]);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'fontsize',7);
    title(['Time: ',num2str(Time(tt)),' min']);
    disp(tt);
    disp(['No. groups counted = ',num2str(length(X.AggregateRadiusFilt))]);

end

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_Timelapse_A'],'-dpng',res);
    close(gcf);
end

%% Time tracks of the two lower concentrations:
% Do regressions of one treatment over time. Track the way that the
% particle attachments scale.
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

datafolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\02_Timelapse\Alt\'];

% Get data:
datalist = dir([datafolder,'*.mat']);

% Inputs:
FilterFLPerArea = 5e3;
FilterSize = [4, Inf];
MicroBeadDiam = 1;
Time = [30,60,90,120,245];
Con  = [5e6, 1e6];

% Set up figure:
figure('units','centimeters','position',[3,3,18,7]);
axs = tight_subplot(length(Con),length(Time), 0, 0.10, 0.05);
axis_counter = 0;

for cc = 1:length(Con)
for tt = 1:length(Time)

    axis_counter = axis_counter + 1;

    % Load data, bin it, and do regressions:
    datalist(axis_counter).name
    X = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(axis_counter).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
    
    % Make figure:
    axes(axs(axis_counter));
    hold on; box on; set(gca,'linewidth',1);
    
    % Plot raw data:
    plot(X.AggregateRadiusFilt, X.NumParticlesFilt, '.','markersize',8,'color',[.7,.7,.7]);

    % Plot data without outliers:
    for bb = 1:length(X.Group_Particles_Out)
        plot(X.Group_Particles_Out{bb}(:,1), X.Group_Particles_Out{bb}(:,2), '.','markersize',8,'color',[.4,.4,.4]);
    end
    
    % Plot bin averaging:
    errorbar(X.Avg_Particles(:,1), X.Bins_Particles,...
        0.6*ones(size(X.Bins_Particles)),  0.6*ones(size(X.Bins_Particles)), ...
        X.Avg_Particles(:,2), X.Avg_Particles(:,2), ...
        'd','markerfacecolor','k','markeredgecolor','none','linewidth',1,'color','k');
    
    % Plot regressions:
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',2,'color',[.8,.2,.2]);
    % plot(X.X_Rad(:,2), X.R_Rad, '-','linewidth',2,'color',[.2,.4,.8]);
    % plot(X.X_ParticlesW(:,2), X.R_ParticlesW, '-','linewidth',2,'color',[.2,.7,.2]);
    
    % Show slope:
    text(6, 50, ['\lambda = ',num2str(X.B_Particles(2),'%0.2f')],'fontsize',8);

    % Aesthetics:
    xlim([5,50]);
    ylim([5e-1,1e2]);
    xticks([10,20,30,40]);
    yticks([1,3,10,30,100]);
    if cc == 2
        xticklabels({'10','20','30','40'});
        xlabel('Radius [\mum]');
    end
    if tt == 1
        yticklabels({'1','3','10','30','100'});
        ylabel('No. particles');
    end
    set(gca,'ticklength',[.03,.03]);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'fontsize',7);
    if cc == 1
        title(['Time: ',num2str(Time(tt)),' min']);
    end

end
end

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_Timelapse_Alt'],'-dpng',res);
    close(gcf);
end

%% Show changing concentrations, over time. Divide by the expected concentration difference to show the overlay

clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

% Get data:
datafolder_con = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\01_ChangingConcentration\New_Data\'];
datalist = dir([datafolder_con,'*.mat']);
dataset1 = [11,6,1]; % 3 different concentrations, each at 100rpm, each at X mins of shaking, each 1um beads
dataset2 = [12,7,2]; % 60 mins
dataset3 = [13,8,3]; % 90 mins
dataset4 = [14,9,4]; % 120 mins
dataset5 = [15,10,5]; % 245 mins
dataset = [dataset1; dataset2; dataset3; dataset4; dataset5]; % concatenate

Colors = abyss(3);
titles = {'30 min','60 min','90 min','120 min','245 min'};

% Inputs:
FilterFLPerArea = 5e3; % [au]
FilterSize      = [4, 100];   % [um]
MicroBeadDiam   = 1;   % [um]

% Generate figure backbone:
figure('units','centimeters','position',[3,3,18,4]);
axs = tight_subplot(1,size(dataset,1), 0.02, 0.17, 0.02);

for tt = 1:size(dataset,1)

    % Make a subfigure:
    axes(axs(tt));
    hold on; box on; set(gca,'linewidth',1);

    % Division parameters:
    DivideBy = [1,5,10];
    disp([titles{tt}]);
    for cc = 1:size(dataset,2)
    
        % Load data, bin it, and do regressions:
        Y = LOAD_DATA_BIN_CALIBRATE([datafolder_con, datalist(dataset(tt,cc)).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
    
        % % Plot the bins:
        errorbar(Y.Avg_Particles(:,1), Y.Bins_Particles/DivideBy(cc),...
            0.6*ones(size(Y.Bins_Particles)),  0.6*ones(size(Y.Bins_Particles)), ...
            Y.Avg_Particles(:,2), Y.Avg_Particles(:,2), ...
            'o','markerfacecolor','none','markeredgecolor',Colors(cc,:),'color',Colors(cc,:),'linewidth',0.5,'markersize',4,'capsize',4);
        % scatter(Y.Avg_Particles(:,1), Y.Bins_Particles/DivideBy(cc), 6,'o','linewidth',0.5,'markeredgecolor',Colors(cc,:));
        
        disp(['C = ',num2str(DivideBy(cc))]);
        disp(['No. groups counted: ',num2str(length(Y.AggregateRadiusFilt))]);

        % Regression:
        % plot(Y.X_Particles(:,2), Y.R_Particles/DivideBy(cc), '-','linewidth',2,'color',Colors(cc,:));
        % plot(Y.X_Rad(:,2), Y.R_Rad, '-','linewidth',2,'color',[.2,.4,.8]);
        
    end

    % Set aesthetics:
    xlim([5,100]);
    xticks([5,10,20,50]);
    xticklabels({'5','10','20','50'});
    ylim([5e-2,5e1]);
    xlabel('Radius [\mum]');
    ylabel('P/C');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'fontsize',7);
    set(gca,'layer','top');
    set(gca,'ticklength',[.03,.03]);
    title(titles{tt});

end

% Print figure:
if PRINT_FIGURES == 1
    print([printfolder,'Supplemental_Concentration'],'-dpng',res);
    close all;
end


%% Showing shaking speed progression regressions:

clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

% Do a regression for each shaking speed tested:
% Get data:
datafolder_SS = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\04_ChangingShakingSpeed\Newest_Data\Mins_060\'];
% datafolder_SS = 'G:\My Drive\01_Data\ProjectEncounter\Microbeads\202412_Nikon_Microbeads\04_ChangingShakingSpeed\Old_Data\Combo\';
datalist_SS = dir([datafolder_SS,'*.mat']);

% Generate figure backbone:
figure;
set(gcf,'units','centimeters','position',[3,2,15,6]);
axs = tight_subplot(1,length(datalist_SS), 0, 0.13, 0.06);

% Inputs:
FilterFLPerArea = 5e3;
FilterSize      = [4, Inf];
MicroBeadDiam   = 1;
SS              = [50,100,150];
% SS = [0,22,47,100,150];

for ss = 1:length(datalist_SS)

    % Load data, bin it, and do regressions:
    datalist_SS(ss).name
    X = LOAD_DATA_BIN_CALIBRATE([datafolder_SS, datalist_SS(ss).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
    
    % Plot this set on the same figure:
    axes(axs(ss));
    hold on; box on; set(gca,'linewidth',1);

    % Plot the raw data:
    plot(X.AggregateRadiusFilt, X.NumParticlesFilt, '.','markersize',8,'color',[.8,.8,.8]);
    for bb = 1:length(X.Group_Particles_Out)
        plot(X.Group_Particles_Out{bb}(:,1), X.Group_Particles_Out{bb}(:,2), '.','markersize',8,'color',[.4,.4,.4]);
    end

    % Plot the bins:
    errorbar(X.Avg_Particles(:,1), X.Bins_Particles,...
        0.6*ones(size(X.Bins_Particles)),  0.6*ones(size(X.Bins_Particles)), ...
        X.Avg_Particles(:,2), X.Avg_Particles(:,2), ...
        'd','markerfacecolor','k','markeredgecolor','none','linewidth',1,'color','k');
    
    % Regression:
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',2,'color',[.8,.2,.2]);

    % Aesthetics:
    xlim([5,50]);
    ylim([5e-1,1e2]);
    xticks([10,20,30,40]);
    yticks([1,3,10,30,100]);
    xticklabels({'10','20','30','40'});
    if ss == 1
        yticklabels({'1','3','10','30','100'});
        ylabel('No. microbeads');
    end
    set(gca,'ticklength',[.03,.03]);
    xlabel('Radius [\mum]');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'fontsize',7);
    title(['Shaking speed: ',num2str(SS(ss)),' rpm']);
    text(6,70,['\lambda=',num2str(X.B_Particles(2)),' \pm ',num2str(2*X.m_std_error,'%0.2f')],'fontsize',7);
    disp(['No. groups counted: ',num2str(length(X.AggregateRadiusFilt))]);

end

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_ShakingSpeed'],'-dpng',res);
    close(gcf);
end

%% Changing shaker:

clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

datafolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\04_ChangingShakingSpeed\ChangingShaker\'];
datalist = dir([datafolder,'*.mat']);

% Properties:
Colors = [0,0,0; .2,.4,.8; .8,.2,.2];
titles = {'Rotation_Well','Rotation_Flask','Side2Side'};

% Inputs:
FilterFLPerArea = 3e3; % [au]
MinSize         = 3.5;   % [um]
FilterSize      = [MinSize, Inf];
MicroBeadDiam   = 1;   % [um]
FIGSAVE         = 0;

% Generate figure backbone:
figure('units','centimeters','position',[3,3,18,6]);
axs = tight_subplot(1,size(datalist,1), 0.03, 0.12, 0.05);

for ii = 1:size(datalist)
    
    % Make a subfigure:
    axes(axs(ii));
    hold on; box on; set(gca,'linewidth',1);

    % Load data, bin it, and do regressions:
    Y = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(ii).name], FilterFLPerArea, FilterSize, MicroBeadDiam);

    % Plot the raw data:
    plot(Y.AggregateRadiusFilt, Y.NumParticlesFilt,'.','markersize',3,'color',[.6,.6,.6]);

    % Plot the bins:
    errorbar(Y.Avg_Particles(:,1), Y.Bins_Particles,...
        0.6*ones(size(Y.Bins_Particles)),  0.6*ones(size(Y.Bins_Particles)), ...
        Y.Avg_Particles(:,2), Y.Avg_Particles(:,2), ...
        'd','markerfacecolor',Colors(ii,:),'markeredgecolor','none','linewidth',1,'color',Colors(ii,:));
    
    % Regression:
    plot(Y.X_Particles(:,2), Y.R_Particles, '-','linewidth',2,'color',Colors(ii,:));

    % Set aesthetics:
    xlim([2,40]);
    xticks([5,10,20,40]);
    xticklabels({'5','10','20','40'});
    ylim([5e-1,2e2]);
    xlabel('Radius [\mum]');
    ylabel('No. microbeads');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'fontsize',7);
    title(titles{ii},'interpreter','none');
    legend({['\lambda = ',num2str(Y.B_Particles(2),'%0.2f'),' \pm ',num2str(2*Y.m_std_error,'%0.2f')]},'location','northwest'); legend('boxoff');

    % Label y-ticks for the first one:
    if ii == 1
        yticks([1,3,10,30,100]);
        yticklabels({'1','3','10','30','100'})
    end

    % More aesthetics:
    set(gca,'ticklength',[.02,.02]);
    disp(['No. groups counted: ',num2str(length(Y.AggregateRadiusFilt))]);
end

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_ChangingShakerMean'],'-dpng',res);
    close(gcf);
end

%% AMINE BEADS experiments
% Load data:
datafolder = 'F:\2025\01_ProjectEncounter\20251203_Revisions_Experiments\EDFs\MATs\combo_all\';
datalist = dir([datafolder, '*comb.mat']);

FL_per_area = 5e3;
FilterSize = [3.5,50];
FilterNoParts = [0, Inf];
FilterCirc = [0.30,1];

% Organize the data:
TimePts = [30,60,90,150]'; % mins
AmineIX = [2,4,6,8];
CarbIX = [1,3,5,7];

% Amine separate:
figure('units','centimeters','position',[3,3,big_width,14]);

% Panel A: Amine beads ----------------------------------------------------
rowa_w = .23;
rowa_h = 0.25;
rowa_x = linspace(0,.92,5) + 0.07;
rowa_y = 0.72;
for tt = 1:length(TimePts)

    ax = axes('position',[rowa_x(tt), rowa_y, rowa_w, rowa_h]);
    amine_ixs = AmineIX(tt);

    % Load the data for the current amine index
    disp(datalist(amine_ixs).name);
    X = LOAD_DATA_BIN_CALIBRATE_Revisions_2([datafolder, datalist(amine_ixs).name], FL_per_area, FilterSize, FilterNoParts, FilterCirc, 1);

    % Plot:
    hold on; box on; set(gca,'linewidth',1);
    plot(X.AggregateRadiusFilt, X.NumParticlesFilt, 'o', 'markersize', 4,'color',[.8,.4,.08]);
    errorbar(X.Avg_Particles(:,1), X.Bins_Particles, 0.6*ones(length(X.Avg_Particles),1), 0.6*ones(length(X.Avg_Particles),1), X.Avg_Particles(:,2), X.Avg_Particles(:,2), 'kd','markersize',4,'markerfacecolor','k','linewidth',1);
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',1,'color',[.8,.2,.2]);    
    
    xlabel('Radius [\mum]');
    yticks([1,10,100]);
    if tt == 1
        ylabel('No. Particles');
    else
        yticklabels({'','',''});
    end
    set(gca,'xscale','log');    
    set(gca,'yscale','log');
    xlim([3,30]);
    ylim([.5,200]);
    title([num2str(X.B_Particles(2),3),' \pm ',num2str(2*X.B_Particles(3),3)])
    set(gca,'fontsize',7);
    set(gca,'ticklength',[.03,.03]);
end

% Panel B - Carboxyl separate: --------------------------------------------
rowb_w = .23;
rowb_h = 0.25;
rowb_x = linspace(0,.92,5) + 0.07;
rowb_y = 0.39;
for tt = 1:length(TimePts)
    axB = axes('position',[rowb_x(tt), rowb_y, rowb_w, rowb_h]);
    carb_ixs = CarbIX(tt);

    % Load the data for the current carboxyl index
    disp(datalist(carb_ixs).name);
    X = LOAD_DATA_BIN_CALIBRATE_Revisions_2([datafolder, datalist(carb_ixs).name], FL_per_area, FilterSize, FilterNoParts, FilterCirc, 1);

    % Plot:
    hold on; box on; set(gca,'linewidth',1);
    plot(X.AggregateRadiusFilt, X.NumParticlesFilt, 'o', 'markersize', 4,'color',[.2,.4,.8]);
    errorbar(X.Avg_Particles(:,1), X.Bins_Particles, 0.6*ones(length(X.Avg_Particles),1), 0.6*ones(length(X.Avg_Particles),1), X.Avg_Particles(:,2), X.Avg_Particles(:,2), 'kd','markersize',4,'markerfacecolor','k','linewidth',1);
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',1,'color',[.8,.2,.2]);    
    
    xlabel('Radius [\mum]');
    yticks([1,10,100]);
    if tt == 1
        ylabel('No. Particles');
    else
        yticklabels({'','',''});
    end
    set(gca,'xscale','log');    
    set(gca,'yscale','log');
    xlim([3,30]);
    ylim([.5,200]);
    title([num2str(X.B_Particles(2),3),' \pm ',num2str(2*X.B_Particles(3),3)])
    set(gca,'fontsize',7);
    set(gca,'ticklength',[.03,.03]);
end

% Panel C -----------------------------------------------------------------
% Put them on the same plot:
% Scale by relative concentration:
% RelativeConcentrationRatio = 1;
RelativeConcentrationRatio = 0.92;
% RelativeConcentrationRatio = 0.85;
rowc_w = .23;
rowc_h = 0.25;
rowc_x = linspace(0,.92,5) + 0.07;
rowc_y = 0.07;
for tt = 1:length(TimePts)
    amine_ixs = AmineIX(tt);
    carb_ixs = CarbIX(tt);

    % Load the data for the current amine index
    disp(datalist(carb_ixs).name);
    X = LOAD_DATA_BIN_CALIBRATE_Revisions_2([datafolder, datalist(carb_ixs).name], FL_per_area, FilterSize, FilterNoParts, FilterCirc, 1);
    Y = LOAD_DATA_BIN_CALIBRATE_Revisions_2([datafolder, datalist(amine_ixs).name], FL_per_area, FilterSize, FilterNoParts, FilterCirc, 1);

    % Plot:
    axC = axes('position',[rowc_x(tt), rowc_y, rowc_w, rowc_h]);
    hold on; box on; set(gca,'linewidth',1);
    % plot(X.AggregateRadiusFilt, X.NumParticlesFilt, 'o', 'markersize', 4,'color',[.2,.4,.8]);
    pC=errorbar(X.Avg_Particles(:,1), X.Bins_Particles, 0.6*ones(length(X.Avg_Particles),1), 0.6*ones(length(X.Avg_Particles),1), X.Avg_Particles(:,2), X.Avg_Particles(:,2), 'd','markersize',4,'markerfacecolor',[.2,.4,.8],'linewidth',1,'markeredgecolor','k','color',[.2,.4,.8]);
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',1.5,'color','k');
    pA=errorbar(Y.Avg_Particles(:,1), Y.Bins_Particles./RelativeConcentrationRatio, 0.6*ones(length(Y.Avg_Particles),1), 0.6*ones(length(Y.Avg_Particles),1), Y.Avg_Particles(:,2), Y.Avg_Particles(:,2), 's','markersize',4,'markerfacecolor',[.8,.4,.1],'linewidth',1,'markeredgecolor','k','color',[.8,.4,.1]);
    plot(Y.X_Particles(:,2), Y.R_Particles./RelativeConcentrationRatio, '-','linewidth',1.5,'color','k'); 
    

    xlabel('Radius [\mum]');
    yticks([1,10,100]);
    if tt == 1
        ylabel('No. Particles');
    else
        yticklabels({'','',''});
    end
    set(gca,'xscale','log');    
    set(gca,'yscale','log');
    xlim([3,30]);
    ylim([.5,200]);
    legend_label_C = ['Carb. m = ',num2str(X.B_Particles(2),2),' \pm ',num2str(2*X.B_Particles(3),2)];
    legend_label_A = ['Amine m = ',num2str(Y.B_Particles(2),2),' \pm ',num2str(2*Y.B_Particles(3),2)];
    legend([pC,pA],{legend_label_C, legend_label_A},'location','northwest','fontsize',5);
    legend('boxoff');
    % title([num2str(X.B_Particles(2),3),' \pm ',num2str(2*X.B_Particles(3),3)])
    set(gca,'fontsize',7);
    set(gca,'ticklength',[.03,.03]);
end

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_AmineBeads'],'-dpng',res);
    close(gcf);
end



%% Poisson predictions, make a full panel of one example:
% Show a histogram plot where each bin is shown separately:
% Inputs:
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

FilterFLPerArea = 5e3;
FilterSize      = [4,50];
MicroBeadDiam   = 1;

fig_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\02_Timelapse\A\'];
X = LOAD_DATA_BIN_CALIBRATE([fig_folder, 'A_090mins_001_EDF.nd2_TD_c_Probabilities.mat'], ...
    FilterFLPerArea, FilterSize, MicroBeadDiam);

% Generate Poisson predictions:
f       = zeros(length(X.Bins_Rad),51);
PCPoiss = cell(length(X.Bins_Rad),1);
MeanR   = zeros(length(X.Bins_Rad),1);
chi2    = zeros(length(X.Bins_Rad),3);
counter = 0;
for bb = 1:length(X.Bins_Rad)
    
    % Get data:
    rr_bb = X.Group_Rad{bb}(:,1); % aggregate radius
    pp_bb = X.Group_Rad{bb}(:,2); % number of particles attached
    MeanR(bb,1) = mean(rr_bb);

    % Get histogram if there are enough datapoints:
    if length(rr_bb) > 10
        counter = counter + 1;

        bin_edges   = -0.5:1:50.5;
        bin_centers = 0:50;
        [f(bb,:), ~] = histcounts(pp_bb, bin_edges, 'normalization','pdf');    
    
        % Find Poisson distribution with the same avg:
        x_poiss = bin_centers;
        y_poiss = poisspdf(x_poiss, X.Avg_Rad(bb,3));
        PCPoiss{bb,1} = [x_poiss; y_poiss];
    
        % Goodness of fit test:
        [Ctemp, ~]  = histcounts(pp_bb, bin_edges,'normalization','count');
        Ytemp       = poisspdf(bin_centers, X.Avg_Rad(bb,3));
        Cnet        = sum(Ctemp); % total number of counts
        Etemp       = Cnet * Ytemp;
        chi2(bb,1)  = sum((Ctemp - Etemp).^2 ./ Etemp); % chi^2 statistic
        chi2(bb,2)  = length(rr_bb) - 1; % degrees of freedom
        chi2(bb,3)  = 1-chi2cdf(chi2(bb,1), chi2(bb,2)); % p-values        
    end
end

% Make a panel-by-panel plot:
Colors = turbo(length(PCPoiss));
NetFigs = counter;
NW = 5;
NT = ceil(NetFigs/NW);
figure('units','centimeters','position',[3,3,big_width,big_width]);
axs = tight_subplot(NT, NW, .01,.01,.01);
for bb = 1:length(PCPoiss)
    if ~isempty(PCPoiss{bb})
        axes(axs(bb)); hold on; box on; set(gca,'linewidth',1);
        plot(0:50, f(bb,:), '.-','linewidth',2,'markersize',12,'color',Colors(bb,:));
        plot(PCPoiss{bb}(1,:), PCPoiss{bb}(2,:), '-','linewidth',1,'markersize',12,'color','k');
        xlim([0,50]);
        ylim([0,1.5*max(f(bb,:))])
        text(5, 1.3*max(f(bb,:)), ['R=',num2str(MeanR(bb),3),'\mum']);
        text(25, 1.1*max(f(bb,:)), ['N=',num2str(chi2(bb,2)+1)],'fontsize',7)
    end
end

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_PoissonPanel'],'-dpng',res);
    close(gcf);
end

%% Show Poisson agreement/disagreement when changing time, concentration, shaking speed and 
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

datafolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\01_ChangingConcentration\New_Data\'];
datalist = dir([datafolder,'*.mat']);

% Inputs:
FilterFLPerArea = 5e3; % [au]
MinSize         = 4;   % [um]
FilterSize      = [MinSize, 50];
MicroBeadDiam   = 1;   % [um]

% Show across time: -------------------------------------------------------
% t1 = 60 mins
% t2 = 120 mins
% t3 = 245 mins
T1 = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(2).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
T2 = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(4).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
T3 = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(5).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
ix1 = find(T1.N_Rad > 5);
ix2 = find(T2.N_Rad > 5);
ix3 = find(T3.N_Rad > 5);
t1s = T1.Avg_Rad(ix1,4).^2./T1.Avg_Rad(ix1,3);
t2s = T2.Avg_Rad(ix2,4).^2./T2.Avg_Rad(ix2,3);
t3s = T3.Avg_Rad(ix3,4).^2./T3.Avg_Rad(ix3,3);

% Figure: -----------------------------------------------------------------
figure('units','centimeters','position',[3,3,big_width,6]);
ax1 = axes('position',[.06,.15,.31,.78]);
hold on; box on; set(gca,'linewidth',1);
plot([0,50],[.5,.5], 'k--','linewidth',1);
plot([0,50],[2,2], 'k--','linewidth',1);
p1=plot(T1.Avg_Rad(ix1, 1), t1s, 'o','markersize',4,'markeredgecolor',[0,0,0],'linewidth',1);
p2=plot(T2.Avg_Rad(ix2, 1), t2s, 's','markersize',5,'markeredgecolor',[.4,.4,.4],'linewidth',1);
p3=plot(T3.Avg_Rad(ix3, 1), t3s, '^','markersize',4,'markeredgecolor',[.7,.7,.7],'linewidth',1);
xlim([0,30]);
ylim([0,9]);
xlabel('Radius [\mum]');
ylabel('Var(P)/E(P)');
legend([p1,p2,p3],{['60min, N=',num2str(length(T1.AggregateRadiusFilt))],['120min, N=',num2str(length(T2.AggregateRadiusFilt))],['245min, N=',num2str(length(T3.AggregateRadiusFilt))]},'location','northwest');
set(gca,'fontsize',7);
xticks([5:5:30]);
yticks(1:2:10);
title('Time');

% Show across concentration: ----------------------------------------------
% c1 = lo
% c2 = med
% c3 = hi
T1 = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(13).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
T2 = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(8).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
T3 = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(3).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
ix1 = find(T1.N_Rad > 5);
ix2 = find(T2.N_Rad > 5);
ix3 = find(T3.N_Rad > 5);
t1s = T1.Avg_Rad(ix1,4).^2./T1.Avg_Rad(ix1,3);
t2s = T2.Avg_Rad(ix2,4).^2./T2.Avg_Rad(ix2,3);
t3s = T3.Avg_Rad(ix3,4).^2./T3.Avg_Rad(ix3,3);

ax2 = axes('position',[.37,.15,.31,.78]);
hold on; box on; set(gca,'linewidth',1);
plot([0,50],[.5,.5], 'k--','linewidth',1);
plot([0,50],[2,2], 'k--','linewidth',1);
p1=plot(T1.Avg_Rad(ix1, 1), t1s, 'o','markersize',4,'markeredgecolor',[0,0,0],'linewidth',1);
p2=plot(T2.Avg_Rad(ix2, 1), t2s, 's','markersize',5,'markeredgecolor',[.4,.4,.4],'linewidth',1);
p3=plot(T3.Avg_Rad(ix3, 1), t3s, '^','markersize',4,'markeredgecolor',[.7,.7,.7],'linewidth',1);
xlim([0,30]);
ylim([0,9]);
xlabel('Radius [\mum]');
legend([p1,p2,p3],{['10^6 [#/mL], N=',num2str(length(T1.AggregateRadiusFilt))],['5*10^6, N=',num2str(length(T2.AggregateRadiusFilt))],['10^7, N=',num2str(length(T3.AggregateRadiusFilt))]},'location','northwest');
set(gca,'fontsize',7);
xticks([5:5:30]);
yticks(1:2:10);
yticklabels({});
title('Concentration');

% Show across shaking speed: ----------------------------------------------
% s1 = lo
% s2 = med
% s3 = hi
datafolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\04_ChangingShakingSpeed\Newest_Data\Mins_060\'];
datalist = dir([datafolder,'*.mat']);
T1 = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(1).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
T2 = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(2).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
T3 = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(3).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
ix1 = find(T1.N_Rad > 5);
ix2 = find(T2.N_Rad > 5);
ix3 = find(T3.N_Rad > 5);
t1s = T1.Avg_Rad(ix1,4).^2./T1.Avg_Rad(ix1,3);
t2s = T2.Avg_Rad(ix2,4).^2./T2.Avg_Rad(ix2,3);
t3s = T3.Avg_Rad(ix3,4).^2./T3.Avg_Rad(ix3,3);

ax1 = axes('position',[.68,.15,.31,.78]);
hold on; box on; set(gca,'linewidth',1);
plot([0,50],[.5,.5], 'k--','linewidth',1);
plot([0,50],[2,2], 'k--','linewidth',1);
p1=plot(T1.Avg_Rad(ix1, 1), t1s, 'o','markersize',4,'markeredgecolor',[0,0,0],'linewidth',1);
p2=plot(T2.Avg_Rad(ix2, 1), t2s, 's','markersize',5,'markeredgecolor',[.4,.4,.4],'linewidth',1);
p3=plot(T3.Avg_Rad(ix3, 1), t3s, '^','markersize',4,'markeredgecolor',[.7,.7,.7],'linewidth',1);
xlim([0,30]);
ylim([0,9]);
xlabel('Radius [\mum]');
legend([p1,p2,p3],{['50 rpm, N=',num2str(length(T1.AggregateRadiusFilt))],['100, N=',num2str(length(T2.AggregateRadiusFilt))],['150, N=',num2str(length(T3.AggregateRadiusFilt))]},'location','northwest');
set(gca,'fontsize',7);
xticks([5:5:30]);
yticks(1:2:10);
yticklabels({});
title('Shaking Speed');

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_TestStatistic'],'-dpng',res);
    close(gcf);
end

%% Show Poisson agreement/disagreement when changing shaking modality:

clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

datafolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\04_ChangingShakingSpeed\ChangingShaker\'];
datalist = dir([datafolder,'*.mat']);

% Properties:
Colors = [0,0,0; .2,.4,.8; .8,.2,.2];
titles = {'Rotation_Well','Rotation_Flask','Side2Side'};

% Inputs:
FilterFLPerArea = 5e3; % [au]
MinSize         = 4;   % [um]
FilterSize      = [MinSize, 50];
MicroBeadDiam   = 1;   % [um]
FIGSAVE         = 0;

% Generate figure backbone:
figure('units','centimeters','position',[3,3,18,6]);
axs = tight_subplot(1,size(datalist,1), 0.03, 0.12, 0.03);
ts = cell(3,1);
for ii = 1:size(datalist)
    
    % Make a subfigure:
    axes(axs(ii));
    hold on; box on; set(gca,'linewidth',1);

    % Load data, bin it, and do regressions:
    Y = LOAD_DATA_BIN_CALIBRATE([datafolder, datalist(ii).name], FilterFLPerArea, FilterSize, MicroBeadDiam);
    [f,PC,~] = DO_POISSON_STUFF(Y, 3, 0);
    plot(2*Y.Avg_Rad(:,1), Y.Avg_Rad(:,3), 'd','markersize',5,'markeredgecolor','none','markerfacecolor','k');
    plot(2*Y.Avg_Rad(:,1), Y.Avg_Rad(:,4).^2, 's','markersize',5,'markeredgecolor','none','markerfacecolor',[.2,.4,.8]);
    set(gca,'xscale','linear');
    set(gca,'yscale','log');
    title(titles{ii},'interpreter','none');
    xlabel('Diameter [\mum]')
    xlim([4,50]);
    xticks([5,10,25,50]);
    xticklabels({'5','10','25','50'});
    ylim([1e-1,3e2]);
    yticks([1,10,30]);
    yticklabels({'1','10','30'});
    set(gca,'ticklength',[.03,.03]);
    set(gca,'fontsize',7);
    L = legend({'Avg.','Variance'},'location','northwest'); legend('boxoff');

    % Make the test statistic:
    ix = find(Y.N_Rad > 5);
    ts{ii}(:,1) = Y.Avg_Rad(ix,1);
    ts{ii}(:,2) = Y.Avg_Rad(ix,4).^2./Y.Avg_Rad(ix,3);
end

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_Poisson_DifferentShakers'],'-dpng',res);
    close(gcf);
end

%%
%{
Colors = [0,0,0; .4,.4,.4; .7,.7,.7];
figure;
hold on; box on; set(gca,'linewidth',1);
plot([0,50],[.5,.5], 'k--','linewidth',1);
plot([0,50],[2,2], 'k--','linewidth',1);
for ii = 1:length(ts)
    p(ii) = plot(ts{ii}(:,1), ts{ii}(:,2), 'o','markersize',4,'markeredgecolor',Colors(ii,:),'linewidth',1);
end
xlim([0,30]);
ylim([0,9]);
xlabel('Radius [\mum]');
ylabel('Var(P)/E(P)');
set(gca,'fontsize',7);
xticks([5:5:30]);
yticks(1:2:10);
title('Changing shaking modality');
legend(p, titles);

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_TestStatisticShakers'],'-dpng',res);
    close(gcf);
end
%}

%% Supplement: What is the relationship between aggregate size and shape
% Thomas C. Day
% April 2025
% We want to know what the discrepancy in the Poisson distribution mean and
% variance at large sizes is caused by. We don't know if this is caused by
% differential stickiness between aggregates, maybe due to aggregates being
% in differential developmental states, or some other reason. One
% alternative hypothesis is that larger aggregates are combinations of
% smaller aggregates, and therefore have a wide range of "lumpinesses".
% This could lead to a larger variance in geometric encounters. Let's test
% this by measuring the size and shape of each aggregate.

% Inputs:
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

datafolder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\01_ChangingConcentration\New_Data\'];
filename = 'A_090mins_001_EDF.nd2_TD_c_Probabilities.mat';
FilterFLPerArea = 5e3;
FilterSize = [4, 50];
MicroBeadDiam = 1; % microns

% T = LOAD_DATA_BIN_CALIBRATE([datafolder, filename], FilterFLPerArea, FilterSize, MicroBeadDiam);

% Load data:
X = load([datafolder, filename]);

% Calibrate:
AggregateRadius = X.XYres * sqrt(X.AggregateArea/pi);

% Pre-allocate for some measurements:
CCListFilt = X.CC.PixelIdxList;
AggregateRadiusFilt1 = AggregateRadius;
AggregateAreaFilt = X.AggregateArea;
NumParticlesFilt = X.NumParticles;

% Filter by fluorescence per unit area:
FL_per_Area = X.IntegratedFL./X.AggregateArea;
ix_filterFL = find(FL_per_Area > FilterFLPerArea);
CCListFilt(ix_filterFL) = [];
AggregateRadiusFilt1(ix_filterFL) = [];
AggregateAreaFilt(ix_filterFL) = [];
NumParticlesFilt(ix_filterFL) = [];

% Filter by size:
ix_filterS = find(AggregateRadiusFilt1 < FilterSize(1) | AggregateRadiusFilt1 > FilterSize(2));
AggregateRadiusFilt1(ix_filterS) = [];
AggregateAreaFilt(ix_filterS) = [];
CCListFilt(ix_filterS) = [];
NumParticlesFilt(ix_filterS) = [];
NumParticlesFilt1 = NumParticlesFilt;

% Make a new CC array:
CCFilt.Connectivity = X.CC.Connectivity;
CCFilt.ImageSize = X.CC.ImageSize;
CCFilt.NumObjects = length(CCListFilt);
CCFilt.PixelIdxList = CCListFilt;

% Now make size and shape measurements on the CC:
Measurements = regionprops(CCFilt, 'Area','EquivDiameter','Circularity','MajorAxisLength','MinorAxisLength','Perimeter');
AggregateRadiusFilt = X.XYres * [Measurements.EquivDiameter]';
AggregateCircularity = [Measurements.Circularity]';
% figure; plot(AggregateRadiusFilt1, AggregateRadiusFilt/2,'.');

% Binning by aggregate radius: --------------------------------------------
nbins_guess = calcnbins(AggregateRadiusFilt,'all');
Nbins = nbins_guess.scott;
Edges_Rad = logspace(log10(min(AggregateRadiusFilt)), log10(max(AggregateRadiusFilt)), Nbins+1);
Bins_Rad  = mean([Edges_Rad(1:end-1); Edges_Rad(2:end)])';
[~,~,BinIx_Rad] = histcounts(AggregateRadiusFilt, Edges_Rad);

% Preallocate for space:
Group_Rad = cell(length(Bins_Rad), 1);
Avg_Rad   = zeros(length(Bins_Rad), 6);
N_Rad     = zeros(length(Bins_Rad), 1);

for bb = 1:length(Bins_Rad)
    ix_Rad = find(BinIx_Rad == bb);
    Group_Rad{bb}(:,1) = AggregateRadiusFilt(ix_Rad);
    Group_Rad{bb}(:,2) = NumParticlesFilt(ix_Rad);
    Group_Rad{bb}(:,3) = AggregateCircularity(ix_Rad);

    % Make averages:
    if length(ix_Rad) > 5
        NotOutliers = find(~isoutlier(Group_Rad{bb}(:,2),'mean'));

        Avg_Rad(bb,1) = mean(Group_Rad{bb}(NotOutliers,1)); % aggregate radius
        Avg_Rad(bb,2) = std(Group_Rad{bb}(NotOutliers,1));
        Avg_Rad(bb,3) = mean(Group_Rad{bb}(NotOutliers,2)); % number of particles
        Avg_Rad(bb,4) = std(Group_Rad{bb}(NotOutliers,2));
        Avg_Rad(bb,5) = mean(Group_Rad{bb}(NotOutliers,3)); % aggregate circularity
        Avg_Rad(bb,6) = std(Group_Rad{bb}(NotOutliers,3));
    else
        Avg_Rad(bb,1) = mean(Group_Rad{bb}(:,1));
        Avg_Rad(bb,2) = std(Group_Rad{bb}(:,1));
        Avg_Rad(bb,3) = mean(Group_Rad{bb}(:,2));
        Avg_Rad(bb,4) = std(Group_Rad{bb}(:,2));
        Avg_Rad(bb,5) = mean(Group_Rad{bb}(:,3));
        Avg_Rad(bb,6) = std(Group_Rad{bb}(:,3));
    end
end
Avg_Rad_Plot = Avg_Rad;
Num_Par_Plot = NumParticlesFilt;

% Try filtering things below a threshold circularity
CircThresh_List = 0:.1:.6;
figure('units','centimeters','position',[3,3,10,10]);
axs = tight_subplot(3,2,.05,.05,.05);
for ii = 1:length(CircThresh_List)
    CircThresh = CircThresh_List(ii);
    ix_filterC = find(AggregateCircularity < CircThresh);
    AggregateRadiusFilt(ix_filterC) = [];
    AggregateCircularity(ix_filterC) = [];
    NumParticlesFilt(ix_filterC) = [];
    
    % Binning by aggregate radius: --------------------------------------------
    nbins_guess = calcnbins(AggregateRadiusFilt,'all');
    Nbins = nbins_guess.scott;
    Edges_Rad = logspace(log10(min(AggregateRadiusFilt)), log10(max(AggregateRadiusFilt)), Nbins+1);
    Bins_Rad  = mean([Edges_Rad(1:end-1); Edges_Rad(2:end)])';
    [~,~,BinIx_Rad] = histcounts(AggregateRadiusFilt, Edges_Rad);
    
    % Preallocate for space:
    Group_Rad = cell(length(Bins_Rad), 1);
    Avg_Rad   = zeros(length(Bins_Rad), 6);
    N_Rad     = zeros(length(Bins_Rad), 1);
    
    for bb = 1:length(Bins_Rad)
        ix_Rad = find(BinIx_Rad == bb);
        Group_Rad{bb}(:,1) = AggregateRadiusFilt(ix_Rad);
        Group_Rad{bb}(:,2) = NumParticlesFilt(ix_Rad);
        Group_Rad{bb}(:,3) = AggregateCircularity(ix_Rad);
    
        % Make averages:
        N_Rad(bb) = length(ix_Rad);
        if length(ix_Rad) > 5
            NotOutliers = find(~isoutlier(Group_Rad{bb}(:,2),'mean'));
    
            Avg_Rad(bb,1) = mean(Group_Rad{bb}(NotOutliers,1)); % aggregate radius
            Avg_Rad(bb,2) = std(Group_Rad{bb}(NotOutliers,1));
            Avg_Rad(bb,3) = mean(Group_Rad{bb}(NotOutliers,2)); % number of particles
            Avg_Rad(bb,4) = std(Group_Rad{bb}(NotOutliers,2));
            Avg_Rad(bb,5) = mean(Group_Rad{bb}(NotOutliers,3)); % aggregate circularity
            Avg_Rad(bb,6) = std(Group_Rad{bb}(NotOutliers,3));
        else
            Avg_Rad(bb,1) = mean(Group_Rad{bb}(:,1));
            Avg_Rad(bb,2) = std(Group_Rad{bb}(:,1));
            Avg_Rad(bb,3) = mean(Group_Rad{bb}(:,2));
            Avg_Rad(bb,4) = std(Group_Rad{bb}(:,2));
            Avg_Rad(bb,5) = mean(Group_Rad{bb}(:,3));
            Avg_Rad(bb,6) = std(Group_Rad{bb}(:,3));
        end
    end

    % Get the residuals:
    residuals = Avg_Rad(:,4).^2 - Avg_Rad(:,3);
    residuals(isnan(residuals)) = [];
    sum_res(ii) = sqrt(sum(residuals.^2))./length(residuals);

    % Make a plot of this:
    if ii < 7
        axes(axs(ii)); hold on; box on; set(gca,'linewidth',1);
        ts = Avg_Rad(:,4).^2./Avg_Rad(:,3);
        % plot(Avg_Rad(:,1), Avg_Rad(:,3), 'kd');
        % plot(Avg_Rad(:,1), Avg_Rad(:,4).^2, 'bs');
        plot(Avg_Rad(:,1), ts, 'ko','markersize',4,'linewidth',1);
        plot([0,75],[.5,.5],'k--','linewidth',1);
        plot([0,75],[2,2],'k--','linewidth',1);
        xlim([0,75]);
        ylim([0,10]);
        title(['Circ thresh = ',num2str(CircThresh)]);
    end

end

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_SizeShapeRelationship_TestStatistic'],'-dpng',res);
    close(gcf);
end

% Plots: ------------------------------------------------------------------
% Show the average and variance of the circularity:
AggDiam = X.XYres * [Measurements.EquivDiameter]';
AggCirc = [Measurements.Circularity]';
figure('units','centimeters','position',[3,3,18,6]); 

% Plot A: Circ vs Radius
ax1 = axes('position',[.08,.15,.25,.80]); hold on; box on; set(gca,'linewidth',1);
plot(AggDiam, AggCirc, '.','markersize',10,'color',[.6,.6,.6]);
errorbar(Avg_Rad_Plot(:,1), Avg_Rad_Plot(:,5), Avg_Rad_Plot(:,6),'>','color',[.8,.4,.04],'markersize',5,'markerfacecolor',[.8,.4,.1],'markeredgecolor','none','linewidth',1,'capsize',6);
xlabel('Group Diameter [\mum]');
ylabel('Circularity');

% PLot B: Num Particles vs. Circ
ax2 = axes('position',[.41,.15,.25,.80]); hold on; box on; set(gca,'linewidth',1);
plot(AggCirc, Num_Par_Plot, '.','markersize',5,'color',[.6,.6,.6]);
% errorbar(Avg_Rad_Plot(:,5), Avg_Rad_Plot(:,3), Avg_Rad_Plot(:,4), '>','markersize',5,'markerfacecolor',[.8,.4,.1],'markeredgecolor','none','color',[.8,.4,.1],'capsize',.6);
xlabel('Circularity');
ylabel('No. Microbeads');

% Plot C: Circ thresholds:
ax3 = axes('position',[.74,.15,.25,.80]); hold on; box on; set(gca,'linewidth',1);
plot(CircThresh_List, sum_res, '>-','linewidth',1,'color',[.8,.4,.1],'markerfacecolor',[.8,.4,.1],'markeredgecolor','none','markersize',5);
xlabel('Circularity Threshold');
ylabel('RMS residual per datapoint');

set(ax1,'fontsize',7);
set(ax2,'fontsize',7);
set(ax3,'fontsize',7);
set(ax1,'layer','top');
set(ax2,'layer','top');
set(ax3,'layer','top');

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_SizeShapeRelationship'],'-dpng',res);
    close(gcf);
end

% Reviewer comments: ------------------------------------------------------
% Bin by circularity class - do we see higher variance by class?
Measurements = regionprops(CCFilt, 'Area','EquivDiameter','Circularity','MajorAxisLength','MinorAxisLength','Perimeter');
AggregateRadiusFilt = 0.5 * X.XYres * [Measurements.EquivDiameter]';
AggregateCircularity = [Measurements.Circularity]';
AggregateElongation = [Measurements.MajorAxisLength]'./[Measurements.MinorAxisLength]';

nclasses = 5;
Circ_Classes = linspace(0,1,nclasses+1);
Measurements_Classes = cell(nclasses,1);
ts = cell(nclasses,1);
Avg_Rad = cell(nclasses,1);
figure('units','centimeters','position',[3,3,10,10]); 
Colors = bone(nclasses+1);
% tiledlayout('flow');
for cc = 1:(length(Circ_Classes)-1)
    ix = find(AggregateCircularity > Circ_Classes(cc) & AggregateCircularity <= Circ_Classes(cc+1));
    if ~isempty(ix)
        Measurements_Classes{cc} = [AggregateRadiusFilt(ix), AggregateCircularity(ix), AggregateElongation(ix), NumParticlesFilt1(ix)];

        % Within each circ bin, bin by aggregate radius and calc mean and
        % variance:
        aggrad = Measurements_Classes{cc}(:,1);
        cirrad = Measurements_Classes{cc}(:,2);
        numrad = Measurements_Classes{cc}(:,4);
        nbins_guess = calcnbins(aggrad,'all');
        Nbins = nbins_guess.scott;
        % Edges_Rad = logspace(log10(min(aggrad)), log10(max(aggrad)), Nbins+1);
        Edges_Rad = linspace(min(aggrad), max(aggrad), Nbins+1);
        Bins_Rad  = mean([Edges_Rad(1:end-1); Edges_Rad(2:end)])';
        [~,~,BinIx_Rad] = histcounts(aggrad, Edges_Rad);
        
        % Preallocate for space:
        Group_Rad   = cell(length(Bins_Rad), 1);
        Avg_Rad{cc} = zeros(length(Bins_Rad), 6);
        N_Rad       = zeros(length(Bins_Rad), 1);
        
        for bb = 1:length(Bins_Rad)
            ix_Rad = find(BinIx_Rad == bb);
            Group_Rad{bb}(:,1) = aggrad(ix_Rad);
            Group_Rad{bb}(:,2) = numrad(ix_Rad);
            Group_Rad{bb}(:,3) = cirrad(ix_Rad);
        
            % Make averages:
            N_Rad(bb) = length(ix_Rad);
            if length(ix_Rad) > 2
                NotOutliers = find(~isoutlier(Group_Rad{bb}(:,2),'mean'));
        
                Avg_Rad{cc}(bb,1) = mean(Group_Rad{bb}(NotOutliers,1)); % aggregate radius
                Avg_Rad{cc}(bb,2) = std(Group_Rad{bb}(NotOutliers,1));
                Avg_Rad{cc}(bb,3) = mean(Group_Rad{bb}(NotOutliers,2)); % number of particles
                Avg_Rad{cc}(bb,4) = std(Group_Rad{bb}(NotOutliers,2));
                Avg_Rad{cc}(bb,5) = mean(Group_Rad{bb}(NotOutliers,3)); % aggregate circularity
                Avg_Rad{cc}(bb,6) = std(Group_Rad{bb}(NotOutliers,3));
            else
                %{
                Avg_Rad(bb,1) = mean(Group_Rad{bb}(:,1));
                Avg_Rad(bb,2) = std(Group_Rad{bb}(:,1));
                Avg_Rad(bb,3) = mean(Group_Rad{bb}(:,2));
                Avg_Rad(bb,4) = std(Group_Rad{bb}(:,2));
                Avg_Rad(bb,5) = mean(Group_Rad{bb}(:,3));
                Avg_Rad(bb,6) = std(Group_Rad{bb}(:,3));
                %}
            end
        end

        % nexttile; 
        hold on; box on; set(gca,'linewidth',1);
        ts{cc} = Avg_Rad{cc}(:,4).^2./Avg_Rad{cc}(:,3);
        plot(Avg_Rad{cc}(:,1), ts{cc}, '.-','color',Colors(cc,:));
    end
end
plot(0:40, 1*ones(1,41),'k:','linewidth',1);
plot(0:40, 2*ones(1,41),'k:','linewidth',1);
xlabel('Radius [\mum]');
ylabel('test statistic');
xticks([5,10,15,20,30]);

figure('units','centimeters','position',[3,3,17.8,12]);
tiledlayout('flow');
for cc = 1:nclasses
    nexttile;
    hold on;
    plot(Measurements_Classes{cc}(:,1), Measurements_Classes{cc}(:,4), '.');
    errorbar(Avg_Rad{cc}(:,1), Avg_Rad{cc}(:,3), Avg_Rad{cc}(:,4), Avg_Rad{cc}(:,4), Avg_Rad{cc}(:,2), Avg_Rad{cc}(:,2), 'kd');
    xlabel('Radius [\mum]');
    ylabel('No. Particles')
    yyaxis right;
    plot(Measurements_Classes{cc}(:,1), Measurements_Classes{cc}(:,2), '.'); ylim([0,1]);
end
        
if PRINT_FIGURES    
    print([printfolder,'Supplemental_BinningByCirc'],'-dpng',res);
    close(gcf);
end



%% Confocal segmentations:
% Show one example segmentation process, then some plots.
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;

fig_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\06_ConfocalSegmentations\'];
Seg = load([fig_folder,'20240801_00-0CM_W_1-5mL_150rpm_1e6_021hr_ Series027_Seg3D_progress.mat']);
load([fig_folder,'DATA_Figure_04_inset.mat']);
XYres = 0.135;
MinSize = 0.05;

% Generate "weird stuff" list to color red in plot:
% Record Cell Volumes and Aspect Ratios fit from Ellipsoids:
CellVol = zeros(length(Seg.Part4.Ellipsoids),1);
CellAsp = CellVol;
for ii = 1:length(Seg.Part4.Ellipsoids)
    CellVol(ii) = XYres^3 * 4/3 * pi * Seg.Part4.Ellipsoids(ii).S(1,1) * Seg.Part4.Ellipsoids(ii).S(2,2) * Seg.Part4.Ellipsoids(ii).S(3,3);
    CellAsp(ii) = Seg.Part4.Ellipsoids(ii).S(1,1) / Seg.Part4.Ellipsoids(ii).S(3,3);
end

% Record cell volumes from a convex hull of the surface:
V = zeros(length(Seg.Part4.Surfaces),1);
for kk = 1:length(Seg.Part4.Surfaces)
    [~,V(kk,1)] = convhull(XYres * Seg.Part4.Surfaces(kk).vertices);
end
ix = find(V < MinSize);

% Weird fits filter:
ix2 = find(CellVol < 0);
ix3 = find(CellAsp < 0);
ix4 = find(CellAsp > 10);
ix_net = unique([ix; ix2; ix3; ix4]);

% Generate colors based on weird ellipsoid fits:
Colors = repmat(.5, length(Seg.Part4.Ellipsoids), 3);
for ii = 1:length(ix_net)
    Colors(ix_net(ii),:) = [1,0,0];
end

% Generate an estimated volume from the max-projection cross-section:
FilteredEllipsoids = Seg.Part4.Ellipsoids;
FilteredEllipsoids(ix_net) = [];
CellCenters = zeros(length(FilteredEllipsoids),3);
CellLengths = zeros(length(FilteredEllipsoids),3);
for kk = 1:length(FilteredEllipsoids)
    CellCenters(kk,:) = XYres * [FilteredEllipsoids(kk).T(1,end), FilteredEllipsoids(kk).T(2,end), FilteredEllipsoids(kk).T(3,end)];
    CellLengths(kk,:) = XYres * [FilteredEllipsoids(kk).S(1,1), FilteredEllipsoids(kk).S(2,2), FilteredEllipsoids(kk).S(3,3)];
end
GroupCOM        = mean(CellCenters);
CellCenters     = CellCenters - GroupCOM;
rCOM            = (CellCenters)./vecnorm(CellCenters, 2, 2);
CellCenters_Sh  = CellCenters + mean(CellLengths(:,1)).*rCOM;
CellCenters_Zproj = CellCenters_Sh(:,1:2);
[K, CrossSectionalArea] = convhull(CellCenters_Zproj(:,1), CellCenters_Zproj(:,2));
estimated_radius = sqrt(CrossSectionalArea/pi);
circle_x = estimated_radius * cos(linspace(0,2*pi,1e3));
circle_y = estimated_radius * sin(linspace(0,2*pi,1e3));

% Generate a regression for estimated volume and number of cells:
c_est = polyfit(GroupVEst, GroupN, 1);
c_est(2) = 0;
x_reg_est = linspace(min(GroupVEst), max(GroupVEst), 100);
y_reg_est = c_est(1) * x_reg_est + c_est(2);
y_est = polyval(c_est, GroupVEst);
residuals = GroupN - y_est;
residuals_std = std(residuals);
ss_res = sum((residuals.^2));
ss_tot = sum((GroupN - mean(GroupN)).^2);
r_squared = 1 - (ss_res./ss_tot)
n = length(GroupN);
x_mean = mean(GroupVEst);
x_var = var(GroupVEst);
m_std_error = residuals_std/sqrt(x_var*(n-1));
b_std_error = sqrt( 1/n * (var(GroupN) * sum(GroupVEst.^2))./ ((n-1)*x_var) );

% Make a panel figure: ----------------------------------------------------
figure('units','centimeters','position',[3,3,18,5]);

% Panel A: Raw image slice
x_location = 65; y_location = 65; Scalebar_length = 2; PixelSize = 0.135;
slice = 50;
axA = axes('position',[.00,.05,.22,.90]);
hold on; box on; set(gca,'linewidth',1);
imagesc(Seg.Part1.Orig(:,:,slice));
axis equal;
xlim([60,170]);
ylim([60,170]);
xticks([]);
yticks([]);
quiver(x_location, y_location, Scalebar_length/PixelSize, 0, 'ShowArrowHead', 'off','color','w','linewidth',3);

% Panel B: 3d surfaces
axB = axes('position',[.22,.05,.22,.90]);
hold on; box on; set(gca,'linewidth',1);
for kk = 1:length(Seg.Part4.Surfaces)
    patch('faces',Seg.Part4.Surfaces(kk).faces,'vertices', 0.135 * Seg.Part4.Surfaces(kk).vertices, 'facecolor',Colors(kk,:), 'edgecolor','none','facealpha',0.7);
end
plot(GroupAS{5}, 'facecolor',[.5,.5,.5],'edgecolor','k','facealpha',0.2);
view(3); axis equal;
camlight; material dull;
xticks([]);
yticks([]);
zticks([]);

% Panel C: Z-projection:
axC = axes('position',[.44,.05,.22,.90]);
hold on; box on; set(gca,'linewidth',1);
plot(CellCenters_Zproj(:,1), CellCenters_Zproj(:,2), 'k.','markersize',12);
plot(CellCenters_Zproj(K,1), CellCenters_Zproj(K,2),'k-','linewidth',2);
patch(circle_x, circle_y, [.5,.5,.5],'edgecolor','none','facealpha',0.4);
plot([0,estimated_radius],[0,0], 'k-','linewidth',1);
axis equal;
xticks([]);
yticks([]);

% Panel D: Measurements of number of cells vs. estimated volume
axD = axes('position',[.73,.16,.25,.80]);
regression_text = ['y = (',num2str(c_est(1)),' \pm ',num2str(m_std_error),') x'];
hold on; box on; set(gca,'linewidth',1);
plot(GroupVEst, GroupN, 'k.','markersize',12);
plot(x_reg_est, y_reg_est, 'k--','linewidth',1);
text(500, max(y_reg_est), regression_text,'fontsize',7);
xlabel('V [\mum^3]');
ylabel('Ncells');
set(gca,'fontsize',7);
% set(gca,'xscale','log');
% set(gca,'yscale','log');

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_ConfocalSegmentations'],'-dpng',res);
    close(gcf);
end


%% Supplemental figure: Effective diffusion constant as a function of run speeds and characteristic runtimes:

v0 = logspace(-1,2,100);
tau = logspace(-1,1,100);
[V0,Tau] = meshgrid(v0, tau);

EffectiveDiffusionConstant = 1/3 * Tau .* V0.^2;
figure('units','centimeters','position',[3,3,10,8]);
hold on; box on; set(gca,'linewidth',1);
set(gca,'xscale','log');
set(gca,'yscale','log');
imagesc(v0, tau, log10(EffectiveDiffusionConstant));
contour(v0, tau, log10(EffectiveDiffusionConstant), [-2,-1,0,1,2,3],'ShowText',true,'color','k','linewidth',1);
set(gca,'layer','top');
xticks([.1,1,10,100,1000]);
yticks([.1,1,10]);
set(gca,'ticklength',[.03,.03]);
xlabel('Swimming Speed, v_0 [\mum/s]');
ylabel('Characteristic run time, \tau [s]');
cb = colorbar;
cb.LineWidth = 1;
colormap('hot');

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_SwimmingStrengths'],'-dpng',res);
    close(gcf);
end

%% What is the kolmogorov length scale for varying energy dissipation rates?
eps = logspace(-11,-3,100);
mu = 1e-6;
Lk = 0.5*(mu^3./eps).^(1/4);
figure; plot(eps, 1e6*Lk); 
set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('\epsilon [W/kg]');
ylabel('Kolmogorov length [um]')

%% Supplement: What happens at very low food concentrations? For the simulations
% Load data:
clearvars -except PRINT_FIGURES drive_folder printfolder res big_width lil_width med_width;
fig_folder = [drive_folder,':\My Drive\04_Manuscripts\02_ProjectEncounter\2025_Feb\07_Sims\'];
X = load([fig_folder,'g0_alpha_sims\sims_MF_17-Jun-2025.mat'],'NumInd','G_0'...
    ,'alpha','n_frac','b_frac','FinalC_FS');
Y = load([fig_folder,'g0_alpha_sims\sims_FS_17-Jun-2025.mat']);
lambda = 7;

% Mean field:
% Plot 1: Number of individuals over time:
Colors = parula(length(X.alpha));
figure('units','centimeters','position',[3,3,20,7]);
ax1 = axes('position',[.05,.15,.25,.80]);
hold on; box on; set(gca,'linewidth',1);
for ll = 1:length(X.alpha)
for gg = 1%:length(X.G_0)
    plot(X.NumInd{gg,ll,1},'-','linewidth',1,'color',Colors(ll,:));
    plot(X.NumInd{gg,ll,2},'-','linewidth',1,'color',Colors(ll,:));
    plot(X.NumInd{gg,ll,3},'-','linewidth',1,'color',Colors(ll,:));
end
end
xlabel('Time');
ylabel('N');
set(gca,'yscale','log');

% Plot 2: Fraction of individuals that are state 2:
ax2 = axes('position',[.375,.15,.25,.80]);
hold on; box on; set(gca,'linewidth',1);
for aa = 1:length(X.alpha)
for gg = 1%:length(X.G_0)
    plot(X.n_frac{gg,aa,1},'-','linewidth',1,'color',Colors(aa,:));
    plot(X.n_frac{gg,aa,2},'-','linewidth',1,'color',Colors(aa,:));
    plot(X.n_frac{gg,aa,3},'-','linewidth',1,'color',Colors(aa,:));
end
end
xlabel('Time');
ylabel('Frac.');
ylim([0,1]);

% Plot 3: Biomass fraction of individuals that are state 2:
ax2 = axes('position',[.70,.15,.25,.80]);
hold on; box on; set(gca,'linewidth',1);
for aa = 1:length(X.alpha)
for gg = 1%:length(X.G_0)
    plot(X.b_frac{gg,aa,1},'-','linewidth',1,'color',Colors(aa,:));
    plot(X.b_frac{gg,aa,2},'-','linewidth',1,'color',Colors(aa,:));
    plot(X.b_frac{gg,aa,3},'-','linewidth',1,'color',Colors(aa,:));
end
end
xlabel('Time');
ylabel('Biomass');
ylim([0,1]);

% FS:
% Plot 1: Number of individuals over time:
Colors = parula(length(X.alpha));
figure('units','centimeters','position',[3,3,20,7]);
ax1 = axes('position',[.05,.15,.25,.80]);
hold on; box on; set(gca,'linewidth',1);
for aa = 1:length(X.alpha)
for gg = 1%:length(Y.G_0)
    plot(Y.NumInd{gg,aa,1},'-','linewidth',1,'color',Colors(aa,:));
    plot(Y.NumInd{gg,aa,2},'-','linewidth',1,'color',Colors(aa,:));
    plot(Y.NumInd{gg,aa,3},'-','linewidth',1,'color',Colors(aa,:));
end
end
xlabel('Time');
ylabel('N');
set(gca,'yscale','log');

% Plot 2: Fraction of individuals that are state 2:
ax2 = axes('position',[.375,.15,.25,.80]);
hold on; box on; set(gca,'linewidth',1);
for aa = 1:length(X.alpha)
for gg = 1%:length(Y.G_0)
    plot(Y.n_frac{gg,aa,1},'-','linewidth',1,'color',Colors(aa,:));
    plot(Y.n_frac{gg,aa,2},'-','linewidth',1,'color',Colors(aa,:));
    plot(Y.n_frac{gg,aa,3},'-','linewidth',1,'color',Colors(aa,:));
end
end
xlabel('Time');
ylabel('Frac.');
ylim([0,1]);

% Plot 3: Biomass fraction of individuals that are state 2:
ax2 = axes('position',[.70,.15,.25,.80]);
hold on; box on; set(gca,'linewidth',1);
for aa = 1:length(X.alpha)
for gg = 1%:length(Y.G_0)
    plot(Y.b_frac{gg,aa,1},'-','linewidth',1,'color',Colors(aa,:));
    plot(Y.b_frac{gg,aa,2},'-','linewidth',1,'color',Colors(aa,:));
    plot(Y.b_frac{gg,aa,3},'-','linewidth',1,'color',Colors(aa,:));
end
end
xlabel('Time');
ylabel('Biomass');
ylim([0,1]);

%% Supplement: How fast do microbial aggregates approach their terminal velocity settling speed when sinking?
%{
% I solved the equations of motion. Make a plot of this:
g       = 10; % [m/s^2]
rho_w   = 1e3; % [kg/m^3]
rho_a   = 1.1e3; % [kg/m^3]
r       = 50*1e-6; % [m]
nu      = 1e-3; % water dynamic viscosity in [Pa*s]
epsilon = 1e-4; % energy dissipation rate in [W/kg]

% Solving for v and plotting:
t = linspace(0, 0.1, 1e4); % time [seconds]
v = (2*g*(rho_a-rho_w)*r^2)/(9*nu) * ( 1 - exp(-(9*nu)/(2*r^2*rho_a)*t) ); % speed
figure; hold on; box on; set(gca,'linewidth',1);
plot(t,v,'k-','linewidth',1);

% Compare timescale to kolmogorov time scale for shaking conditions:
t_sinking = (2*r^2*rho_a)/(9*nu);
t_kolmogorov = sqrt(nu/(rho_w*epsilon));
disp(['Sinking timescale : ',num2str(t_sinking,4)]);
disp(['Kolmogorov timescale : ',num2str(t_kolmogorov,4)]);
%}

%% Supplement: What is the relationship between maximum occupancy density and attachment rate?
p = (0:300)';
r = linspace(1,100,1e2);
[P,R] = meshgrid(p,r);
S = P .* (4*pi*R.^2).^(-1);

PHI = [.01,.1,.3,.9];
% PHI = 0.02;

figure('units','centimeters','position',[3,3,10,10]);
tiledlayout('flow');
for ii = 1:length(PHI)
    phi0 = PHI(ii);
    a0 = pi*0.5^2; % 1 micron diameter disk area
    sigma0(ii) = phi0/a0;
    
    Stemp = S./sigma0(ii);
    Stemp(Stemp > 1) = 1;
    Stemp(5,1:10);
    Stemp(10,1:10);
    Stemp(20,1:25);

    nexttile;
    hold on; box on; set(gca,'linewidth',1);
    contourf(R,P,Stemp,[0,.01,.1,.25,1],'k--','linewidth',1,'showtext','on');
    xlabel('Group radius [\mum]');
    ylabel('p');
    xlim([0,50]);
    ylim([0,100]);
    xticks([0,10,20,30,40,50]);
    yticks(0:10:300)
    set(gca,'layer','top');
    colormap('sky');
    title(['\phi = ',num2str(PHI(ii))]);
end

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_Langmuir'],'-dpng',res);
    close(gcf);
end



%% Supplement: Theory diagrams for phage encounters via turbulence and diffusion separately:
% Sweeping two variables:
n_cells = logspace(0,4,100); % number of cells in a group
n_p = logspace(4,7,50); % phage concentration
radius_phage = 1e-5; % cm
correction_factor = 1;

% Two models: Diffusion only, or turbulence only
for nn = 1:length(n_cells)
    for pp = 1:length(n_p)
        Gamma_D(nn,pp) = computePhageAcquisitionRate(n_cells(nn), n_p(pp));
        Gamma_T(nn,pp) = computeFoodAcquisitionRate(n_cells(nn), radius_phage, n_p(pp));
    end
end
Gamma_N = Gamma_T + Gamma_D;

% Consider the timescale to get X cells from a single cell - is this possible?
Colormap = [0,60,48; 
            1,102,94;
            53,151,143;
            128,205,193;
            199,234,229;
            245,245,245;
            246,232,195;
            223,194,125;
            191,129,45;
            140,81,10;
            84,48,5]/255;
doubling_time = 30*60; % secs
assembly_time = (log2(n_cells) * doubling_time)/3600; % in hrs
Assembly_Time = repmat(assembly_time', [1,length(n_p)]);
A = correction_factor*(1./Gamma_D);
B = correction_factor*(1./Gamma_N);
TimeMatrix_Hrs_D = A./3600;
TimeMatrix_Hrs_N = B./3600;
AssemblyMatrix = PossibilityMatrix./Assembly_Time;

figure('units','centimeters','position',[3,3,15,8]);
ax1 = axes('position',[.08,.10,.45,.88]);
hold on; box on; set(gca,'linewidth',1);
imagesc(n_p, n_cells, TimeMatrix_Hrs_D);
colormap(ax1, Colormap);
[C,h] = contour(n_p, n_cells, TimeMatrix_Hrs_D,[20,60,180,600]/60,'w-','linewidth',1.5,'ShowText',true,'LabelFormat','%2.2f [hrs]');
clabel(C,h,'fontsize',8,'color','w');
% text(2e4,3e0,'Can assemble','fontsize',8,'color','w');
% text(1.5e6, 2e3,{'Cannot','assemble'},'fontsize',8,'color','w');
% clim([0,10]);
clim([0,10]);
xlim([1e4,1e7]);
ylim([1,1e4]);
yticks([1,10,100,1e3,1e4]);
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'layer','top');
set(gca,'fontsize',7);
xlabel('Infect. Phage Conc [#/mL]');
ylabel('Group Biomass [# Cells]');
set(gca,'ticklength',[.03,.03]);

ax2 = axes('position',[.53,.10,.45,.88]);
hold on; box on; set(gca,'linewidth',1);
imagesc(n_p, n_cells, TimeMatrix_Hrs_N);
colormap(ax2, Colormap);
[C,h] = contour(n_p, n_cells, TimeMatrix_Hrs_N,[20,60,180,600]/60,'w-','linewidth',1.5,'ShowText',true,'LabelFormat','%2.2f [hrs]');
clabel(C,h,'fontsize',8,'color','w');
% text(2e4,3e0,'Can assemble','fontsize',8,'color','w');
% text(1.5e6, 2e3,{'Cannot','assemble'},'fontsize',8,'color','w');
% clim([0,10]);
clim([0,10]);
xlim([1e4,1e7]);
ylim([1,1e4]);
yticks([1,10,100,1e3,1e4]);
yticklabels({'','','','',''});
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'layer','top');
set(gca,'fontsize',7);
xlabel('Infect. Phage Conc [#/mL]');
set(gca,'ticklength',[.03,.03]);

% Print figure:
if PRINT_FIGURES    
    print([printfolder,'Supplemental_InfectionTimescales_DiffTurb'],'-dpng',res);
    close(gcf);
end










%{
%% Supplement: AMINE BEADS concentrations of beads
datafolder = 'E:\2025\01_ProjectEncounter\20251119_MicrobeadsRevisions\Bead_Concentration_Tests\Nikon\10x\';
img_folder = 'E:\2025\01_ProjectEncounter\20251119_MicrobeadsRevisions\Bead_Concentration_Tests\Confocal\';
load([datafolder,'BeadConcentration_Nikon10x.mat']);
Ni = ni(ixs);
Nbar = mean(Ni,2);
Nerr = std(Ni,0,2);
N_type_stats = [Nbar, Nerr];
x_type = 2:3:length(filelist);
y_type = repelem(x_type,3);
r_type = repmat([-.5,0,.5],1,8);
y_type = y_type + r_type;

figure('units','centimeters','position',[3,3, lil_width, 11])
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on; box on; set(gca,'linewidth',1);
bar(x_type(1:4), Nbar(1:4), 0.85, 'facecolor',[.8,.4,.1],'edgecolor','k','facealpha',1,'linewidth',1);
bar(x_type(5:end), Nbar(5:end), 0.85, 'facecolor',[.2,.4,.8],'edgecolor','k','facealpha',1,'linewidth',1);
errorbar(x_type, Nbar, Nerr, 'k.','markersize',1,'linewidth',1,'capsize',6);
plot(y_type, STATS(:,1), 'ko','linewidth',1);
xticks(x_type);
xticklabels({'T=0, Alg','T=90, Alg','T=0, noc','T=90, noc','T=0, alg','T=90, alg','T=0, noc','T=90, noc'});
ylabel('No. objects');
xlim([0,25]);
set(gca,'layer','top');
set(gca,'fontsize',7);

nexttile; hold on; box on; set(gca,'linewidth',1);
errorbar(y_type(1:12), STATS(1:12,2), STATS(1:12,3), 'pentagram','linewidth',1,'capsize',6,'color',[.8,.4,.1]);
errorbar(y_type(13:end), STATS(13:end,2), STATS(13:end,3), 'o','linewidth',1,'capsize',6,'color',[.2,.4,.8]);
xticks(x_type);
xticklabels({'T=0, Alg','T=90, Alg','T=0, noc','T=90, noc','T=0, alg','T=90, alg','T=0, noc','T=90, noc'});
ylabel('Area [px]');
xlim([0,25]);
set(gca,'fontsize',7);

% Determine what the relative difference in concentration is:
RelativeConcentrationRatio = mean(Nbar(1:4))./mean(Nbar(5:end));
disp(['Relative concentration ratio = ',num2str(RelativeConcentrationRatio)]);

% Show the amine beads under confocal:
ax = axes('position',[.55,.26,.38,.20]); 
% hold on; box on; set(gca,'linewidth',1);
img = imread([img_folder,'Montage_Examples_Amine_Alg_90mins_5um_scalebar.tif']);
imshow(img); colormap('gray');
xticks([]);
yticks([]);




%% Supplement: AMINE BEADS supplemental panels
datafolder = 'E:\2025\01_ProjectEncounter\20251119_MicrobeadsRevisions\MATs\';
datalist = dir([datafolder, '*comb.mat']);

FL_per_area = 4e3;
FilterSize = [4,100];
FilterNoParts = [0, 0.12];

% Organize the data:
TimePts = [30,60,90,120,250]'; % mins
AmineIX = [1,3,5,7,9];
CarbIX = [2,4,6,8,10];

% Amine separate:
figure('units','centimeters','position',[3,3,big_width,7]);
tiledlayout(2,5,'TileSpacing','none','Padding','compact');
for tt = 1:length(TimePts)
    amine_ixs = AmineIX(tt);

    % Load the data for the current amine index
    disp(datalist(amine_ixs).name);
    X = LOAD_DATA_BIN_CALIBRATE_Revisions([datafolder, datalist(amine_ixs).name], FL_per_area, FilterSize, FilterNoParts, 1);

    % Plot:
    nexttile; 
    hold on; box on; set(gca,'linewidth',1);
    plot(X.AggregateRadiusFilt, X.NumParticlesFilt, 'pentagram', 'markersize', 4,'color',[.8,.4,.08]);
    errorbar(X.Avg_Particles(:,1), X.Bins_Particles, 0.6*ones(length(X.Avg_Particles),1), 0.6*ones(length(X.Avg_Particles),1), X.Avg_Particles(:,2), X.Avg_Particles(:,2), 'kd','markersize',4,'markerfacecolor','k','linewidth',1);
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',1,'color',[.8,.2,.2]);    
    text(5,50,['\lambda=',num2str(X.B_Particles(2),2),'\pm',num2str(2*X.B_Particles(3),2)],'fontsize',7);

    if tt == 1
        ylabel('No. microbeads');
        yticks([1,10,100]);
        yticklabels({1,10,100});
    else
        yticks([1,10,100]);
        yticklabels({'','',''});
        xticks([10,100]);
        xticklabels({'',''});
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'fontsize',7);
    xlim([FilterSize(1), FilterSize(2)]);
    ylim([5e-1,1e2]);
    % title([num2str(X.B_Particles(2)),' \pm ',num2str(X.B_Particles(3))])
end

% Carboxyl separate:
for tt = 1:length(TimePts)
    carb_ixs = CarbIX(tt);

    % Load the data for the current amine index
    disp(datalist(carb_ixs).name);
    X = LOAD_DATA_BIN_CALIBRATE_Revisions([datafolder, datalist(carb_ixs).name], FL_per_area, FilterSize, FilterNoParts, 1);

    % Plot:
    nexttile; hold on; box on; set(gca,'linewidth',1);
    plot(X.AggregateRadiusFilt, X.NumParticlesFilt, 'o', 'markersize', 4,'color',[.2,.4,.8]);
    errorbar(X.Avg_Particles(:,1), X.Bins_Particles, 0.6*ones(length(X.Avg_Particles),1), 0.6*ones(length(X.Avg_Particles),1), X.Avg_Particles(:,2), X.Avg_Particles(:,2), 'kd','markersize',4,'markerfacecolor','k','linewidth',1);
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',1,'color',[.8,.2,.2]);    
    text(5,60,['\lambda=',num2str(X.B_Particles(2),2),'\pm',num2str(2*X.B_Particles(3),2)],'fontsize',7);
    
    if tt == 1
        ylabel('No. microbeads');
        yticks([1,10,100]);
        yticklabels({1,10,100});
    else
        yticks([1,10,100]);
        yticklabels({'','',''});
    end
    xlabel('Radius [\mum]');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'fontsize',7);
    xlim([FilterSize(1), FilterSize(2)]);
    ylim([5e-1,1e2]);
    % title([num2str(X.B_Particles(2)),' \pm ',num2str(X.B_Particles(3))])
end

% Put avgs on the same plot:
RelativeConcentrationRatio = 1;
figure('units','centimeters','position',[3,3,big_width,7]);
tiledlayout(2,5,'TileSpacing','none','Padding','compact');
for tt = 1:length(TimePts)
    amine_ixs = AmineIX(tt);
    carb_ixs = CarbIX(tt);

    % Load the data for the current amine index
    disp(datalist(carb_ixs).name);
    X = LOAD_DATA_BIN_CALIBRATE_Revisions([datafolder, datalist(carb_ixs).name], FL_per_area, FilterSize, FilterNoParts, 1);
    Y = LOAD_DATA_BIN_CALIBRATE_Revisions([datafolder, datalist(amine_ixs).name], FL_per_area, FilterSize, FilterNoParts, 1);

    % Plot:
    nexttile; hold on; box on; set(gca,'linewidth',1);
    errorbar(X.Avg_Particles(:,1), X.Bins_Particles, 0.6*ones(length(X.Avg_Particles),1), 0.6*ones(length(X.Avg_Particles),1), X.Avg_Particles(:,2), X.Avg_Particles(:,2), 'o','markersize',4,'markerfacecolor',[.2,.4,.8],'linewidth',1,'markeredgecolor','k','color',[.2,.4,.8]);
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',1.5,'color','k');
    errorbar(Y.Avg_Particles(:,1), Y.Bins_Particles./RelativeConcentrationRatio, 0.6*ones(length(Y.Avg_Particles),1), 0.6*ones(length(Y.Avg_Particles),1), Y.Avg_Particles(:,2), Y.Avg_Particles(:,2), 'pentagram','markersize',6,'markerfacecolor',[.8,.4,.1],'linewidth',1,'markeredgecolor','k','color',[.8,.4,.1]);
    plot(Y.X_Particles(:,2), Y.R_Particles./RelativeConcentrationRatio, '-','linewidth',1.5,'color','k'); 
    

    if tt == 1
        ylabel('No. microbeads');
        yticks([1,10,100]);
        yticklabels({1,10,100});
    else
        yticks([1,10,100]);
        yticklabels({'','',''});
    end
    % xlabel('Radius [\mum]');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'fontsize',7);
    xlim([FilterSize(1), 50]);
    ylim([5e-1,1e2]);
    xticks([5,10,25]);
    xticklabels({'','',''});
end

% Scale by relative concentration:
RelativeConcentrationRatio = 0.60444;
for tt = 1:length(TimePts)
    amine_ixs = AmineIX(tt);
    carb_ixs = CarbIX(tt);

    % Load the data for the current amine index
    disp(datalist(carb_ixs).name);
    X = LOAD_DATA_BIN_CALIBRATE_Revisions([datafolder, datalist(carb_ixs).name], FL_per_area, FilterSize, FilterNoParts, 1);
    Y = LOAD_DATA_BIN_CALIBRATE_Revisions([datafolder, datalist(amine_ixs).name], FL_per_area, FilterSize, FilterNoParts, 1);

    % Plot:
    nexttile; hold on; box on; set(gca,'linewidth',1);
    errorbar(X.Avg_Particles(:,1), X.Bins_Particles, 0.6*ones(length(X.Avg_Particles),1), 0.6*ones(length(X.Avg_Particles),1), X.Avg_Particles(:,2), X.Avg_Particles(:,2), 'o','markersize',4,'markerfacecolor',[.2,.4,.8],'linewidth',1,'markeredgecolor','k','color',[.2,.4,.8]);
    plot(X.X_Particles(:,2), X.R_Particles, '-','linewidth',1.5,'color','k');
    errorbar(Y.Avg_Particles(:,1), Y.Bins_Particles./RelativeConcentrationRatio, 0.6*ones(length(Y.Avg_Particles),1), 0.6*ones(length(Y.Avg_Particles),1), Y.Avg_Particles(:,2), Y.Avg_Particles(:,2), 'pentagram','markersize',6,'markerfacecolor',[.8,.4,.1],'linewidth',1,'markeredgecolor','k','color',[.8,.4,.1]);
    plot(Y.X_Particles(:,2), Y.R_Particles./RelativeConcentrationRatio, '-','linewidth',1.5,'color','k'); 
    

    if tt == 1
        ylabel('No. microbeads');
        yticks([1,10,100]);
        yticklabels({1,10,100});
    else
        yticks([1,10,100]);
        yticklabels({'','',''});
    end
    xlabel('Radius [\mum]');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    set(gca,'fontsize',7);
    xlim([FilterSize(1), 50]);
    ylim([5e-1,1e2]);
    xticks([5,10,25]);
end



%% Supplement: AMINE BEADS, directly compare power law exponents:
datafolder = 'E:\2025\01_ProjectEncounter\20251119_MicrobeadsRevisions\MATs\';
t_list = dir([datafolder,'*comb.mat']);

FilterFLPerArea = 4e3;
FilterSize = [4,100];
FilterNoParts = [1e-4, 0.12];

TimePts = [30,60,90,120,250]'; % mins
AmineIX = [1,3,5,7,9];
CarbIX = [2,4,6,8,10];

% Load timetracks of carboxyl beads: --------------------------------------
for tt = 1:length(CarbIX)
    % Load data for carboxyl beads
    W = LOAD_DATA_BIN_CALIBRATE_Revisions([datafolder, t_list(CarbIX(tt)).name], FilterFLPerArea, FilterSize, FilterNoParts, 1);
    W_Carboxyl(tt) = W;
end

% Load timetracks of amine beads: -----------------------------------------
for tt = 1:length(AmineIX)
    W = LOAD_DATA_BIN_CALIBRATE_Revisions([datafolder, t_list(AmineIX(tt)).name], FilterFLPerArea, FilterSize, FilterNoParts, 1);
    W_Amine(tt) = W;
end

% Specifically compare the regression slopes:
Colors = [.2,.4,.8; .8,.4,.1];
figure('units','centimeters','position',[3,3,lil_width,lil_width]); 
hold on; box on; set(gca,'linewidth',1);
B_c = [W_Carboxyl.B_Particles];
B_a = [W_Amine.B_Particles];
errorbar(TimePts, B_c(2,:), 2*B_c(3,:), 'o','markersize',7,'markerfacecolor',Colors(1,:),'color',Colors(1,:),'linewidth',1,'markeredgecolor','k');
errorbar(TimePts, B_a(2,:), 2*B_a(3,:), 'pentagram','markersize',9,'markerfacecolor',Colors(2,:),'color',Colors(2,:),'linewidth',1,'markeredgecolor','k');
xlabel('Time [min]');
ylabel('\lambda');

% Welch's t-test to compare the slopes:
s_delta = sqrt(B_c(3,:).^2 + B_a(3,:).^2); % combined std devs
t = (B_a(2,:) - B_c(2,:))./s_delta;
%}













%% Functions
function [X] = LOAD_DATA_BIN_CALIBRATE(filename, FilterFLPerArea, FilterSize, MicrobeadDiam)
    % Loading, filtering, calibrating: ------------------------------------
    % Load data:
    X = load(filename);

    % Calibrate:
    AggregateRadius = X.XYres * sqrt(X.AggregateArea/pi);

    % Filter by fluorescence per unit area:
    FL_per_Area = X.IntegratedFL ./ X.AggregateArea;
    ix_filterFL = find(FL_per_Area < FilterFLPerArea);
    X.AggregateRadiusFilt = AggregateRadius(ix_filterFL);
    X.NumParticlesFilt    = X.NumParticles(ix_filterFL);
    X.IntegratedFLFilt    = X.IntegratedFL(ix_filterFL);

    % Filter by size:
    ix_filterS = find(X.AggregateRadiusFilt > FilterSize(1) & X.AggregateRadiusFilt < FilterSize(2));
    X.AggregateRadiusFilt = X.AggregateRadiusFilt(ix_filterS);
    X.NumParticlesFilt = X.NumParticlesFilt(ix_filterS);
    X.IntegratedFLFilt = X.IntegratedFLFilt(ix_filterS);

    % Show how many aggregate were filtered this way:
    disp(['Total aggregates filtered: ',num2str(length(X.AggregateArea) - length(X.AggregateRadiusFilt)), ' / ',num2str(length(X.AggregateArea))])

    % Binning: ------------------------------------------------------------
    % Binning by number of particles encountered:
    Edges_Particles = -.5:1:max(X.NumParticlesFilt)+.5;
    Bins_Particles  = (0:max(X.NumParticlesFilt))';
    [~, ~, BinIx_Particles] = histcounts(X.NumParticlesFilt, Edges_Particles);

    % Pre-allocate space:
    Group_Particles     = cell(length(Bins_Particles), 1);
    Group_Particles_Out = Group_Particles;
    Avg_Particles       = zeros(length(Bins_Particles), 2);
    Fluo_Particles      = zeros(length(Bins_Particles), 2);
    N_Particles         = zeros(length(Bins_Particles), 1);

    % Loop through particle bins:
    for bb = 1:length(Bins_Particles)
        ix_Particles = find(BinIx_Particles == bb);
        Group_Particles{bb}(:,1) = X.AggregateRadiusFilt(ix_Particles);
        Group_Particles{bb}(:,2) = X.NumParticlesFilt(ix_Particles);
        Group_Particles{bb}(:,3) = X.IntegratedFLFilt(ix_Particles);

        % Find outliers within this bin:
        NotOutliers = find(~isoutlier(Group_Particles{bb}(:,1), 'mean'));
        Group_Particles_Out{bb} = Group_Particles{bb}(NotOutliers, :);

        % Get a histogram:
        N_Particles(bb) = length(NotOutliers);

        % Measure averages on non-outliers:
        Avg_Particles(bb,1) = mean(Group_Particles_Out{bb}(:,1));
        Avg_Particles(bb,2) = std(Group_Particles_Out{bb}(:,1));
        Fluo_Particles(bb,1) = mean(Group_Particles_Out{bb}(:,3));
        Fluo_Particles(bb,2) = std(Group_Particles_Out{bb}(:,3));
    end
    X.Group_Particles     = Group_Particles;
    X.Group_Particles_Out = Group_Particles_Out;
    X.N_Particles         = N_Particles;
    X.Avg_Particles       = Avg_Particles;
    X.Bins_Particles      = Bins_Particles;
    X.Fluo_Particles      = Fluo_Particles;

    % Binning by aggregate radius: ----------------------------------------
    nbins_guess = calcnbins(X.AggregateRadiusFilt,'all');
    Nbins = nbins_guess.scott;
    % Nbins = 10;
    % Edges_Rad = linspace(min(X.AggregateRadiusFilt), max(X.AggregateRadiusFilt), Nbins);
    Edges_Rad = logspace(log10(min(X.AggregateRadiusFilt)), log10(max(X.AggregateRadiusFilt)), Nbins);
    Bins_Rad  = mean([Edges_Rad(1:end-1); Edges_Rad(2:end)])';
    [~,~,BinIx_Rad] = histcounts(X.AggregateRadiusFilt, Edges_Rad);
    Group_Rad = cell(length(Bins_Rad),1);
    Avg_Rad   = zeros(length(Bins_Rad),2);
    Fluo_Rad  = zeros(length(Bins_Rad),2);
    N_Rad     = zeros(length(Bins_Rad),1);
    for bb = 1:length(Bins_Rad)
        ix_Rad = find(BinIx_Rad == bb);
        Group_Rad{bb}(:,1) = X.AggregateRadiusFilt(ix_Rad);
        Group_Rad{bb}(:,2) = X.NumParticlesFilt(ix_Rad);
        Group_Rad{bb}(:,3) = X.IntegratedFLFilt(ix_Rad);

        % In bins where there are sufficient datapoints, get averages for
        % Poisson statistics:
        if length(ix_Rad) > 5
            NotOutliers = find(~isoutlier(Group_Rad{bb}(:,2),'mean'));

            Avg_Rad(bb,1) = mean(Group_Rad{bb}(NotOutliers,1));
            Avg_Rad(bb,2) = std(Group_Rad{bb}(NotOutliers,1));
            Avg_Rad(bb,3) = mean(Group_Rad{bb}(NotOutliers,2));
            Avg_Rad(bb,4) = std(Group_Rad{bb}(NotOutliers,2));
            Fluo_Rad(bb,1) = mean(Group_Rad{bb}(NotOutliers,3));
            Fluo_Rad(bb,2) = std(Group_Rad{bb}(NotOutliers,3));
        else
            Avg_Rad(bb,1) = mean(Group_Rad{bb}(:,1));
            Avg_Rad(bb,2) = std(Group_Rad{bb}(:,1));
            Avg_Rad(bb,3) = mean(Group_Rad{bb}(:,2));
            Avg_Rad(bb,4) = std(Group_Rad{bb}(:,2));
            Fluo_Rad(bb,1) = mean(Group_Rad{bb}(:,3));
            Fluo_Rad(bb,2) = std(Group_Rad{bb}(:,3));
        end
        N_Rad(bb) = length(ix_Rad);
    end
    X.Group_Rad = Group_Rad;
    X.N_Rad     = N_Rad;
    X.Avg_Rad   = Avg_Rad;
    X.Bins_Rad  = Bins_Rad;
    X.Fluo_Rad  = Fluo_Rad;

    % Regression: ---------------------------------------------------------
    [X.X_Particles, X.Y_Particles, X.R_Particles, X.B_Particles, X.m_std_error] = REGRESSION_AVGS(X.Avg_Particles(:,1), X.Bins_Particles, MicrobeadDiam);
    [X.X_Rad, X.Y_Rad, X.R_Rad, X.B_Rad, ~] = REGRESSION_AVGS(X.Avg_Rad(:,1), X.Avg_Rad(:,3), MicrobeadDiam);

end

function [X] = LOAD_DATA_BIN_CALIBRATE_TIMETRACKS(filename, FilterFLPerArea, FilterSize, MicrobeadDiam)
    % Use a constant size binning algorithm to match bins across time.
    % Loading, filtering, calibrating: ------------------------------------
    % Load data:
    X = load(filename);

    % Calibrate:
    AggregateRadius = X.XYres * sqrt(X.AggregateArea/pi);

    % Filter by fluorescence per unit area:
    FL_per_Area = X.IntegratedFL ./ X.AggregateArea;
    ix_filterFL = find(FL_per_Area < FilterFLPerArea);
    X.AggregateRadiusFilt = AggregateRadius(ix_filterFL);
    X.NumParticlesFilt    = X.NumParticles(ix_filterFL);
    X.IntegratedFLFilt    = X.IntegratedFL(ix_filterFL);

    % Filter by size:
    ix_filterS = find(X.AggregateRadiusFilt > FilterSize);
    X.AggregateRadiusFilt = X.AggregateRadiusFilt(ix_filterS);
    X.NumParticlesFilt = X.NumParticlesFilt(ix_filterS);
    X.IntegratedFLFilt = X.IntegratedFLFilt(ix_filterS);

    % Show how many aggregate were filtered this way:
    disp(['Total aggregates filtered: ',num2str(length(X.AggregateArea) - length(X.AggregateRadiusFilt)), ' / ',num2str(length(X.AggregateArea))])

    % Binning: ------------------------------------------------------------
    % Binning by number of particles encountered:
    Edges_Particles = -.5:1:max(X.NumParticlesFilt)+.5;
    Bins_Particles  = (0:max(X.NumParticlesFilt))';
    [~, ~, BinIx_Particles] = histcounts(X.NumParticlesFilt, Edges_Particles);

    % Pre-allocate space:
    Group_Particles     = cell(length(Bins_Particles), 1);
    Group_Particles_Out = Group_Particles;
    Avg_Particles       = zeros(length(Bins_Particles), 2);
    Fluo_Particles      = zeros(length(Bins_Particles), 2);
    N_Particles         = zeros(length(Bins_Particles), 1);

    % Loop through particle bins:
    for bb = 1:length(Bins_Particles)
        ix_Particles = find(BinIx_Particles == bb);
        Group_Particles{bb}(:,1) = X.AggregateRadiusFilt(ix_Particles);
        Group_Particles{bb}(:,2) = X.NumParticlesFilt(ix_Particles);
        Group_Particles{bb}(:,3) = X.IntegratedFLFilt(ix_Particles);
        
        % Find outliers within this bin:
        NotOutliers = find(~isoutlier(Group_Particles{bb}(:,1), 'mean'));
        Group_Particles_Out{bb} = Group_Particles{bb}(NotOutliers, :);

        % Measure averages on non-outliers:
        Avg_Particles(bb,1) = mean(Group_Particles_Out{bb}(:,1));
        Avg_Particles(bb,2) = std(Group_Particles_Out{bb}(:,1));
        Fluo_Particles(bb,1) = mean(Group_Particles_Out{bb}(:,3));
        Fluo_Particles(bb,2) = std(Group_Particles_Out{bb}(:,3));
    end
    X.Group_Particles     = Group_Particles;
    X.Group_Particles_Out = Group_Particles_Out;
    X.N_Particles         = N_Particles;
    X.Avg_Particles       = Avg_Particles;
    X.Bins_Particles      = Bins_Particles;
    X.Fluo_Particles      = Fluo_Particles;

    % Binning by aggregate radius: ----------------------------------------
    nbins_guess = calcnbins(X.AggregateRadiusFilt,'all');
    Nbins = nbins_guess.scott;
    Nbins = 50;
    % Nbins = 10;
    % Edges_Rad = linspace(min(X.AggregateRadiusFilt), max(X.AggregateRadiusFilt), Nbins);
    % Edges_Rad = logspace(log10(min(X.AggregateRadiusFilt)), log10(max(X.AggregateRadiusFilt)), Nbins);
    Edges_Rad = logspace(log10(3), log10(100), Nbins+1);
    % Edges_Rad = linspace(3,100,Nbins+1);
    Bins_Rad  = mean([Edges_Rad(1:end-1); Edges_Rad(2:end)])';
    [~,~,BinIx_Rad] = histcounts(X.AggregateRadiusFilt, Edges_Rad);
    Group_Rad = cell(length(Bins_Rad),1);
    Avg_Rad   = zeros(length(Bins_Rad),2);
    Fluo_Rad  = zeros(length(Bins_Rad),2);
    N_Rad     = zeros(length(Bins_Rad),1);
    for bb = 1:length(Bins_Rad)
        ix_Rad = find(BinIx_Rad == bb);
        Group_Rad{bb}(:,1) = X.AggregateRadiusFilt(ix_Rad);
        Group_Rad{bb}(:,2) = X.NumParticlesFilt(ix_Rad);
        Group_Rad{bb}(:,3) = X.IntegratedFLFilt(ix_Rad);

        % In bins where there are sufficient datapoints, get averages for
        % Poisson statistics:
        if length(ix_Rad) > 5
            NotOutliers = find(~isoutlier(Group_Rad{bb}(:,2),'mean'));

            Avg_Rad(bb,1) = mean(Group_Rad{bb}(NotOutliers,1));
            Avg_Rad(bb,2) = std(Group_Rad{bb}(NotOutliers,1));
            Avg_Rad(bb,3) = mean(Group_Rad{bb}(NotOutliers,2));
            Avg_Rad(bb,4) = std(Group_Rad{bb}(NotOutliers,2));
            Fluo_Rad(bb,1) = mean(Group_Rad{bb}(NotOutliers,3));
            Fluo_Rad(bb,2) = std(Group_Rad{bb}(NotOutliers,3));
        else
            Avg_Rad(bb,1)  = nan;
            Avg_Rad(bb,2)  = nan;
            Avg_Rad(bb,3)  = nan;
            Avg_Rad(bb,4)  = nan;
            Fluo_Rad(bb,1) = nan;
            Fluo_Rad(bb,2) = nan;
        end
        N_Rad(bb) = length(ix_Rad);
    end
    X.Group_Rad = Group_Rad;
    X.N_Rad     = N_Rad;
    X.Avg_Rad   = Avg_Rad;
    X.Bins_Rad  = Bins_Rad;
    X.Fluo_Rad  = Fluo_Rad;

    % Regression: ---------------------------------------------------------
    [X.X_Particles, X.Y_Particles, X.R_Particles, X.B_Particles, m_std_error] = REGRESSION_AVGS(X.Avg_Particles(:,1), X.Bins_Particles, MicrobeadDiam);
    [X.X_Rad, X.Y_Rad, X.R_Rad, X.B_Rad, m_std_error] = REGRESSION_AVGS(X.Avg_Rad(:,1), X.Avg_Rad(:,3), MicrobeadDiam);
    X.Residuals = X.R_Particles - X.Y_Particles;
    
    % Wholesale regression:
    %{
    x1 = [];
    y1 = [];
    for bb = 1:length(X.Group_Particles_Out)
        x1 = [x1; X.Group_Particles_Out{bb}(:,1)];
        y1 = [y1; X.Group_Particles_Out{bb}(:,2)];
    end
    [X.X_ParticlesW, X.Y_ParticlesW, X.R_ParticlesW, X.B_ParticlesW] = REGRESSION_NET(x1, y1, MicrobeadDiam);
    %}
end

function [XData, YData, RegData, B, m_std_error] = REGRESSION_AVGS(XData, YData, MicrobeadDiam)
% Perform a power-law regression on the log-transformed dataset:

    % % Identify outliers and remove them before regression:
    % NoOutliers = find(~isoutlier(XData,'mean'));
    % XData = XData(NoOutliers);
    % YData = YData(NoOutliers);
    
    XData      = log10(XData + MicrobeadDiam/2);
    [XData,Ix] = sort(XData);
    YData      = log10(YData);
    YData      = YData(Ix);
    
    % Trim without nans:
    NoNans = find(~isnan(XData));
    XData  = XData(NoNans);
    YData  = YData(NoNans);

    % Trim without infinities:
    NoInf = find(~isinf(YData));
    XData = XData(NoInf);
    YData = YData(NoInf);

    % Make special form of XData for regression:
    XData  = [ones(length(XData),1), XData];

    % Regression:
    if length(YData) > 25
        B = XData(1:25,:)\YData(1:25);
    elseif length(YData) > 10
        B = XData(1:end-5,:)\YData(1:end-5);
    else
        B = XData \ YData;
    end
    RegData = XData*B;

    % Get standard error in slope:
    res = YData - RegData;
    xbar = XData(:,2) - mean(XData(:,2));
    m_std_error = sqrt( (length(res)-2).^(-1) * (sum(res.^2)/sum(xbar.^2)) );

    % Return fit:
    XData   = 10.^XData;
    YData   = 10.^YData;
    RegData = 10.^RegData;
end

function nbins = calcnbins(x, method, minimum, maximum)
% Calculate the "ideal" number of bins to use in a histogram, using a
% choice of methods.
% 
% NBINS = CALCNBINS(X, METHOD) calculates the "ideal" number of bins to use
% in a histogram, using a choice of methods.  The type of return value
% depends upon the method chosen.  Possible values for METHOD are:
% 'fd': A single integer is returned, and CALCNBINS uses the
% Freedman-Diaconis method,
% based upon the inter-quartile range and number of data.
% See Freedman, David; Diaconis, P. (1981). "On the histogram as a density
% estimator: L2 theory". Zeitschrift f�r Wahrscheinlichkeitstheorie und
% verwandte Gebiete 57 (4): 453-476.
% 'scott': A single integer is returned, and CALCNBINS uses Scott's method,
% based upon the sample standard deviation and number of data.
% See Scott, David W. (1979). "On optimal and data-based histograms".
% Biometrika 66 (3): 605-610.
% 
% 'sturges': A single integer is returned, and CALCNBINS uses Sturges'
% method, based upon the number of data.
% See Sturges, H. A. (1926). "The choice of a class interval". J. American
% Statistical Association: 65-66.
% 
% 'middle': A single integer is returned.  CALCNBINS uses all three
% methods, then picks the middle (median) value.
% 
% 'all': A structure is returned with fields 'fd', 'scott' and 'sturges',
% each containing the calculation from the respective method.
% 
% NBINS = CALCNBINS(X) works as NBINS = CALCNBINS(X, 'MIDDLE').
% 
% NBINS = CALCNBINS(X, [], MINIMUM), where MINIMUM is a numeric scalar,
% defines the smallest acceptable number of bins.
% 
% NBINS = CALCNBINS(X, [], MAXIMUM), where MAXIMUM is a numeric scalar,
% defines the largest acceptable number of bins.
% 
% Notes: 
% 1. If X is complex, any imaginary components will be ignored, with a
% warning.
% 
% 2. If X is an matrix or multidimensional array, it will be coerced to a
% vector, with a warning.
% 
% 3. Partial name matching is used on the method name, so 'st' matches
% sturges, etc.
% 
% 4. This function is inspired by code from the free software package R
% (http://www.r-project.org).  See 'Modern Applied Statistics with S' by
% Venables & Ripley (Springer, 2002, p112) for more information.
% 
% 5. The "ideal" number of depends on what you want to show, and none of
% the methods included are as good as the human eye.  It is recommended
% that you use this function as a starting point rather than a definitive
% guide.
% 
% 6. The wikipedia page on histograms currently gives a reasonable
% description of the algorithms used.
% See http://en.wikipedia.org/w/index.php?title=Histogram&oldid=232222820
% 
% Examples:     
% y = randn(10000,1);
% nb = calcnbins(y, 'all')
%    nb = 
%             fd: 66
%          scott: 51
%        sturges: 15
% calcnbins(y)
%    ans =
%        51
% subplot(3, 1, 1); hist(y, nb.fd);
% subplot(3, 1, 2); hist(y, nb.scott);
% subplot(3, 1, 3); hist(y, nb.sturges);
% y2 = rand(100,1);
% nb2 = calcnbins(y2, 'all')
%    nb2 = 
%             fd: 5
%          scott: 5
%        sturges: 8
% hist(y2, calcnbins(y2))
% 
% See also: HIST, HISTX
% 
% $ Author: Richard Cotton $		$ Date: 2008/10/24 $    $ Version 1.5 $
% Input checking
error(nargchk(1, 4, nargin));
if ~isnumeric(x) && ~islogical(x)
    error('calcnbins:invalidX', 'The X argument must be numeric or logical.')
end
if ~isreal(x)
   x = real(x);
   warning('calcnbins:complexX', 'Imaginary parts of X will be ignored.');
end
% Ignore dimensions of x.
if ~isvector(x)
   x = x(:);
   warning('calcnbins:nonvectorX', 'X will be coerced to a vector.');
end
nanx = isnan(x);
if any(nanx)
   x = x(~nanx);
   warning('calcnbins:nanX', 'Values of X equal to NaN will be ignored.');
end
if nargin < 2 || isempty(method)
   method = 'middle';
end
if ~ischar(method)
   error('calcnbins:invalidMethod', 'The method argument must be a char array.');
end
validmethods = {'fd'; 'scott'; 'sturges'; 'all'; 'middle'};
methodmatches = strmatch(lower(method), validmethods);
nmatches = length(methodmatches);
if nmatches~=1
   error('calnbins:unknownMethod', 'The method specified is unknown or ambiguous.');
end
method = validmethods{methodmatches};
if nargin < 3 || isempty(minimum)
   minimum = 1;
end
if nargin < 4 || isempty(maximum)
   maximum = Inf;
end
   
% Perform the calculation
switch(method)
   case 'fd'
      nbins = calcfd(x);
   case 'scott'
      nbins = calcscott(x);
    case 'sturges'
      nbins = calcsturges(x);
   case 'all'
      nbins.fd = calcfd(x);    
      nbins.scott = calcscott(x);
      nbins.sturges = calcsturges(x);
   case 'middle'
      nbins = median([calcfd(x) calcscott(x) calcsturges(x)]);
end
% Calculation details
   function nbins = calcfd(x)
      h = diff(prctile0(x, [25; 75])); %inter-quartile range
      if h == 0
         h = 2*median(abs(x-median(x))); %twice median absolute deviation
      end
      if h > 0
         nbins = ceil((max(x)-min(x))/(2*h*length(x)^(-1/3)));
      else
         nbins = 1;
      end
      nbins = confine2range(nbins, minimum, maximum);
   end
   function nbins = calcscott(x)
      h = 3.5*std(x)*length(x)^(-1/3);
      if h > 0 
         nbins = ceil((max(x)-min(x))/h);
      else 
         nbins = 1;
      end
      nbins = confine2range(nbins, minimum, maximum);
   end
   function nbins = calcsturges(x)
      nbins = ceil(log2(length(x)) + 1);
      nbins = confine2range(nbins, minimum, maximum);
   end
   function y = confine2range(x, lower, upper)
      y = ceil(max(x, lower));
      y = floor(min(y, upper));
   end
   function y = prctile0(x, prc)
      % Simple version of prctile that only operates on vectors, and skips
      % the input checking (In particluar, NaN values are now assumed to
      % have been removed.)
      lenx = length(x);
      if lenx == 0
         y = [];
         return
      end
      if lenx == 1
         y = x;
         return
      end
      
      function foo = makecolumnvector(foo)
         if size(foo, 2) > 1 
            foo = foo';
         end
      end
         
      sortx = makecolumnvector(sort(x));
      posn = prc.*lenx/100 + 0.5;
      posn = makecolumnvector(posn);
      posn = confine2range(posn, 1, lenx);
      y = interp1q((1:lenx)', sortx, posn);
   end
end

function [f, PCPoiss, chi2] = DO_POISSON_STUFF(X, Nmin, FigViz)
    % Thomas C. Day
    % This function loops through particular size classes and records the mean
    % number of particles encountered, then generates a Poisson prediction of
    % the distribution given that mean.

    %% Calculations:
    f       = cell(length(X.Bins_Rad),1);
    PCPoiss = cell(length(X.Bins_Rad),1);
    chi2    = zeros(length(X.Bins_Rad),3);
    counter = 0;
    for bb = 1:length(X.Bins_Rad)
        % Get data:
        rr_bb = X.Group_Rad{bb}(:,1); % aggregate radius
        pp_bb = X.Group_Rad{bb}(:,2); % number of particles attached
    
        % Get histogram if there are enough datapoints:
        if length(rr_bb) > Nmin
            counter = counter + 1;
    
            bin_edges   = -0.5:1:(max(pp_bb)+3);
            bin_centers = mean([bin_edges(1:end-1); bin_edges(2:end)]);
            [fx, ~] = histcounts(pp_bb, bin_edges, 'normalization','pdf');    
            f{bb} = [bin_centers', fx'];

            % Find Poisson distribution with the same avg:
            x_poiss = 0:50;
            y_poiss = poisspdf(x_poiss, X.Avg_Rad(bb,3));
            PCPoiss{bb,1} = [x_poiss; y_poiss];
        
            % Goodness of fit test:
            [Ctemp, ~]  = histcounts(pp_bb, bin_edges,'normalization','count');
            Ytemp       = poisspdf(bin_centers, X.Avg_Rad(bb,3));
            Cnet        = sum(Ctemp); % total number of counts
            Etemp       = Cnet * Ytemp;
            chi2(bb,1)  = sum((Ctemp - Etemp).^2 ./ Etemp); % chi^2 statistic
            chi2(bb,2)  = length(rr_bb) - 1; % degrees of freedom
            chi2(bb,3)  = 1-chi2cdf(chi2(bb,1), chi2(bb,2)); % p-values     

            % Different goodness of fit:
            ECDF = ecdf(pp_bb);
            PCDF = poisscdf(bin_centers, X.Avg_Rad(bb,3));

        end
    end


    %% Plots:
    if FigViz == 1
        % Show mean and variance on the same plot:
        figure('units','centimeters','position',[38,3,7,7]); 
        hold on; box on; set(gca,'linewidth',1);
        plot(X.Avg_Rad(:,1), X.Avg_Rad(:,3), 'd','markersize',5,'markeredgecolor','none','markerfacecolor','k');
        plot(X.Avg_Rad(:,1), X.Avg_Rad(:,4).^2, 's','markersize',5,'markeredgecolor','none','markerfacecolor',[.2,.4,.8]);
        text(X.Avg_Rad(:,1), X.Avg_Rad(:,3), num2str(chi2(:,2)));
        xlabel('Radius [\mu m]');
        ylabel('Measure');
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        
        % Show a histogram plot where each bin is shown separately:
        Colors = turbo(length(PCPoiss));
        NetFigs = length(X.Bins_Rad);
        NW = 6;
        NT = ceil(NetFigs/NW);
        figure('units','centimeters','position',[18,3,20,20]);
        axs = tight_subplot(NT, NW, .01,.01,.01);
        for bb = 1:length(PCPoiss)
            if ~isempty(PCPoiss{bb})
                axes(axs(bb)); hold on; box on; set(gca,'linewidth',1);
                plot(0:50, f(bb,:), '.-','linewidth',2,'markersize',12,'color',Colors(bb,:));
                plot(PCPoiss{bb}(1,:), PCPoiss{bb}(2,:), '-','linewidth',1,'color','k');
                % xlim([0,25]);
            end
        end
    elseif FigViz==2
        % Show stacked histogram plot:
        AddConstant = 0.05;
        Colors = turbo(length(PCPoiss));
        figure('units','centimeters','position',[3,3,15,20]); 
        hold on; box on; set(gca,'linewidth',1);
        for bb = 1:length(PCPoiss)
            if ~isempty(PCPoiss{bb})
                plot(0:50, f(bb,:) + AddConstant*bb, '.-','linewidth',2,'markersize',12,'color',Colors(bb,:));
                plot(PCPoiss{bb}(1,:), PCPoiss{bb}(2,:) + AddConstant*bb, '-','linewidth',1,'markersize',12,'color','k');
            end
        end
        xlabel('No. Beads');
        ylabel('PDF');
    elseif FigViz==3
        figure('units','centimeters','position',[3,20,7,7]); 
        hold on; box on; set(gca,'linewidth',1);
        yyaxis left;
        plot(X.Avg_Rad(:,1), chi2(:,2), '-');
        yyaxis right;
        plot(X.Avg_Rad(:,1), chi2(:,3), '-');
        xlabel('Radius [\mu m]');
    end
end

function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
    % tight_subplot creates "subplot" axes with adjustable gaps and margins
    %
    % [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
    %
    %   in:  Nh      number of axes in hight (vertical direction)
    %        Nw      number of axes in width (horizontaldirection)
    %        gap     gaps between the axes in normalized units (0...1)
    %                   or [gap_h gap_w] for different gaps in height and width 
    %        marg_h  margins in height in normalized units (0...1)
    %                   or [lower upper] for different lower and upper margins 
    %        marg_w  margins in width in normalized units (0...1)
    %                   or [left right] for different left and right margins 
    %
    %  out:  ha     array of handles of the axes objects
    %                   starting from upper left corner, going row-wise as in
    %                   subplot
    %        pos    positions of the axes objects
    %
    %  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
    %           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
    %           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
    % Pekka Kumpulainen 21.5.2012   @tut.fi
    % Tampere University of Technology / Automation Science and Engineering
    if nargin<3; gap = .02; end
    if nargin<4 || isempty(marg_h); marg_h = .05; end
    if nargin<5; marg_w = .05; end
    if numel(gap)==1; 
        gap = [gap gap];
    end
    if numel(marg_w)==1; 
        marg_w = [marg_w marg_w];
    end
    if numel(marg_h)==1; 
        marg_h = [marg_h marg_h];
    end
    axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
    axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
    py = 1-marg_h(2)-axh; 
    % ha = zeros(Nh*Nw,1);
    ii = 0;
    for ih = 1:Nh
        px = marg_w(1);
        
        for ix = 1:Nw
            ii = ii+1;
            ha(ii) = axes('Units','normalized', ...
                'Position',[px py axw axh], ...
                'XTickLabel','', ...
                'YTickLabel','');
            px = px+axw+gap(2);
        end
        py = py-axh-gap(1);
    end
    if nargout > 1
        pos = get(ha,'Position');
    end
    ha = ha(:);
end