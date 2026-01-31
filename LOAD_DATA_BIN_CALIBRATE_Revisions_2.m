function [X] = LOAD_DATA_BIN_CALIBRATE_Revisions_2(filename, FilterFLPerArea, FilterSize, FilterNoParts, FilterCirc, MicrobeadDiam)
%Thomas C. Day
% Loading data, binnning into histograms, etc.

    % Load:
    X = load(filename);

    % Calibrate:
    AggregateRadius = X.XYres * sqrt(X.AggregateArea/pi); % convert to microns

    % Filter by fluorescence:
    FL_per_Area = X.IntegratedFL ./ X.AggregateArea;
    ix_FilterFL = find(FL_per_Area < FilterFLPerArea);
    X.AggregateRadiusFilt = AggregateRadius(ix_FilterFL);
    X.NumParticlesFilt = X.NumParticles(ix_FilterFL);
    X.IntegratedFLFilt = X.IntegratedFL(ix_FilterFL);
    X.AggregateCircFilt = X.AggregateCirc(ix_FilterFL);

    % Filter by size:
    ix_FilterS = find(X.AggregateRadiusFilt > FilterSize(1) & X.AggregateRadiusFilt < FilterSize(2));
    X.AggregateRadiusFilt = X.AggregateRadiusFilt(ix_FilterS);
    X.NumParticlesFilt = X.NumParticlesFilt(ix_FilterS);
    X.IntegratedFLFilt = X.IntegratedFLFilt(ix_FilterS);
    X.AggregateCircFilt = X.AggregateCircFilt(ix_FilterS);

    % Filter by no. of particles:
    NperArea = X.NumParticlesFilt./(pi*X.AggregateRadiusFilt.^2);
    ix_FilterN = find(NperArea >= FilterNoParts(1) & NperArea < FilterNoParts(2));
    X.AggregateRadiusFilt = X.AggregateRadiusFilt(ix_FilterN);
    X.NumParticlesFilt = X.NumParticlesFilt(ix_FilterN);
    X.IntegratedFLFilt = X.IntegratedFLFilt(ix_FilterN);
    X.AggregateCircFilt = X.AggregateCircFilt(ix_FilterN);

    % Filter by circularity:
    ix_FilterC = find(X.AggregateCircFilt >= FilterCirc(1) & X.AggregateCircFilt < FilterCirc(2)); % Example threshold for circularity
    X.AggregateRadiusFilt = X.AggregateRadiusFilt(ix_FilterC);
    X.NumParticlesFilt = X.NumParticlesFilt(ix_FilterC);
    X.IntegratedFLFilt = X.IntegratedFLFilt(ix_FilterC);
    X.AggregateCircFilt = X.AggregateCircFilt(ix_FilterC);

    % Display number of aggregate filtered this way:
    disp(['Total aggregates filtered = ',num2str(length(X.AggregateArea) - length(X.AggregateRadiusFilt))]);

    % Binning -------------------------------------------------------------
    % Binning by number of particles encountered:
    Edges_Particles = -.5:1:max(X.NumParticlesFilt + .5);

    Bins_Particles  = (0:max(X.NumParticlesFilt))';
    [~,~, BinIx_Particles] = histcounts(X.NumParticlesFilt, Edges_Particles);

    % Pre-allocate space for a bunch of arrays:
    Group_Particles = cell(length(Bins_Particles), 3);
    Group_Particles_Out = Group_Particles;
    Avg_Particles = zeros(length(Bins_Particles), 2);
    Fluo_Particles = Avg_Particles;
    N_Particles = zeros(length(Bins_Particles),1);

    % Loop through particle bins and calculate/assign:
    for bb = 1:length(Bins_Particles)
        ix_Particles = find(BinIx_Particles == bb);
        N_Particles(bb) = length(ix_Particles);
        Group_Particles{bb}(:,1) = X.AggregateRadiusFilt(ix_Particles);
        Group_Particles{bb}(:,2) = X.NumParticlesFilt(ix_Particles);
        Group_Particles{bb}(:,3) = X.IntegratedFLFilt(ix_Particles);

        % Find outliers within this particular bin:
        NotOutliers = find(~isoutlier(Group_Particles{bb}(:,1),'mean'));
        Group_Particles_Out{bb} = Group_Particles{bb}(NotOutliers,:);

        % Measure averages on non-outliers:
        Avg_Particles(bb,:) = [mean(Group_Particles_Out{bb}(:,1)), std(Group_Particles_Out{bb}(:,1))];
        Fluo_Particles(bb,:) = [mean(Group_Particles_Out{bb}(:,3)), std(Group_Particles_Out{bb}(:,3))];
    end
    X.Group_Particles = Group_Particles;
    X.Group_Particles_Out = Group_Particles_Out;
    X.N_Particles = N_Particles;
    X.Avg_Particles = Avg_Particles;
    X.Bins_Particles = Bins_Particles;
    X.Fluo_Particles = Fluo_Particles;

    % Binning by aggregate radius: ----------------------------------------
    % nbins_guess = calcnbins(X.AggregateRadiusFilt,'all');
    % Nbins = nbins_guess.scott;
    % Edges_Rad = linspace(min(X.AggregateRadiusFilt), max(X.AggregateRadiusFilt), Nbins);
    % Edges_Rad = logspace(log10(min(X.AggregateRadiusFilt)), log10(max(X.AggregateRadiusFilt)), Nbins);
    Edges_Rad = logspace(log10(3), log10(100), 40);
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

    % Regression:
    [X.X_Particles, X.Y_Particles, X.R_Particles, X.B_Particles] = REGRESSION_AVGS(X.Avg_Particles(:,1), X.Bins_Particles, MicrobeadDiam);

end