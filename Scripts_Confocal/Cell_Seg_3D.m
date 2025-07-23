% Thomas C. Day
% Segmentation script 20240418

%% INPUTS
RUNSTEPS        = 2;
FIGVIZ          = 2;
SAVEPROG        = 0;
NET_THRESH      = 15;       % Initial binarization to get the whole group and nothing extra
INT_THRESH      = 0.25;     % Niblack threshold value (0.5 is pretty good generally)
VOL_THRESH      = 100;      % Size threshold (in voxels) of a cell

% filepath = 'F:\00_USC\01_Data_ProjectAggregateDevelopment\Webb_PointScanConfocal\20240925good\ISOs\6E02\';
filepath = 'E:\My Drive\01_Data\Project_Encounter\PackingMeasurements\';
filelist = dir([filepath,'00_Raw\','*.tif']);

for ff = 1:length(filelist)
    % Progress bar:
    disp(repmat('-',1,50));
    disp(['Segmenting file: ',num2str(ff)]);

    % Custom image stack segmentation:
    CELL_SEGMENTATION_TCD(filepath, filelist(ff).name, NET_THRESH, INT_THRESH, VOL_THRESH, SAVEPROG, RUNSTEPS, FIGVIZ);
end

%% Functions
function CELL_SEGMENTATION_TCD(filepath, filename, NET_THRESH, INT_THRESH, VOL_THRESH, SAVEPROG, RUNSTEPS, FIGVIZ)

    %% PART 0: Initial Preprocessing:
    
    disp('Preprocessing...');

    % Load data:
    savename = [filename(1:end-4),'_Seg3D_progress.mat'];
    orig     = tiffreadVolume([filepath,'00_Raw\',filename]);
    zslice   = 1:size(orig,3);
    
    % Big blur and binarize to get whole group:
    gb = fspecial('gaussian',11,5);
    BW = false(size(orig));
    for zz = 1:size(orig,3)
        temp = imfilter(orig(:,:,zz),gb);
        BW(:,:,zz) = temp > NET_THRESH;
    end
    
    % Dilate a bit to make sure everything is obtained:
    BW = imdilate(BW, strel('cube',3));
    BW = imfill(BW,'holes');
    if FIGVIZ == 1
        volshow(BW);
    end
    
    % Take the largest connected component:
    CC = bwconncomp(BW);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    filtered_BW = false(size(BW));
    filtered_BW(CC.PixelIdxList{idx}) = true;
    if FIGVIZ == 1
        volshow(filtered_BW);
    end
    
    % Mask the original:
    M1 = uint8(filtered_BW).*orig;
    
    % Write to file:
    write_tiff_stack(M1,[filepath,'01_Preprocess\',filename]);
    
    %% PART 1: Initial Binarization, slice by slice
    
    if RUNSTEPS > 0
    disp('Initial Binarization...');

    % Load pre-processed-image:
    orig = tiffreadVolume([filepath,'01_Preprocess\',filename]);
    
    % Binarize each z-slice separately:
    BW = false(size(orig));
    for zz = zslice
        % disp(['Binarizing z-slice: ',num2str(zz)]); % Progress bar
        tempPV = double(orig(:,:,zz));
    
        % Niblack thresholding:
        EX  = conv2(tempPV, fspecial('average',15),'same');
        EX2 = conv2(tempPV.^2, fspecial('average',15),'same');
        img_mu = EX;
        img_si = sqrt(EX2 - EX.^2);
        NiblackThresh = img_mu + INT_THRESH*img_si;
        BW_niblack = tempPV > NiblackThresh;
        BW(:,:,zz) = BW_niblack;
    end
    
    if FIGVIZ == 1
        figure;
        for zz = zslice
            B = bwboundaries(BW(:,:,zz));
            clf;
            imagesc(orig(:,:,zz)); hold on; colormap('gray');
            for bb = 1:length(B)
                if length(B{bb}) > 5
                    B{bb}(:,1) = sgolayfilt(B{bb}(:,1), 2, 5);
                    B{bb}(:,2) = sgolayfilt(B{bb}(:,2), 2, 5);
                    plot(B{bb}(:,2), B{bb}(:,1),'r-','linewidth',1);
                else
                    plot(B{bb}(:,2), B{bb}(:,1),'r-','linewidth',1);
                end
            end
            axis equal;
            drawnow;
        end
    end
    
    if SAVEPROG
        Part1.Orig   = orig;
        Part1.BW     = BW;
        Part1.Zslice = zslice;
        Part1.Net_Thresh = NET_THRESH;
        Part1.Int_Thresh = INT_THRESH;
        Part1.Vol_Thresh = VOL_THRESH;
        save([filepath,'02_Segs\',savename],'Part1');
    end
    
    end
    
    %% PART 2: A full 3D watershed
    
    if RUNSTEPS > 1
    disp('Watershed...');

    % Load last step:
    load([filepath,'02_Segs\',savename]);
    BW = Part1.BW;
    TI = Part1.Orig;
    
    % Distance transform:
    D = bwdist(~BW,'cityblock');
    D = round(medfilt3(D));
    D = imdilate(D, strel('cube',3));

    % Blur original by a bit:
    BI = medfilt3(TI);

    % Combine fluorescence image and distance transformed image for a
    % watershed-ready image:
    WI = double(BI) .* double(D);

    % Watershed:
    W = watershed(-WI, 6);
    L = W==0;
    
    % Cut along watershed lines:
    BW(L) = 0;
    TI(L) = 0;
    LI = bwlabeln(BW,6);
    
    % Remove items below a size threshold:
    BW2 = bwareaopen(BW, VOL_THRESH, 6);
    LI2 = bwlabeln(BW2,6);

    % Show:
    if FIGVIZ == 2
        CMAP = rand(max(LI,[],'all'), 3);
        RGB  = label2rgb3d(LI,CMAP,[0,0,0]);
        volshow(RGB);
        CMAP2 = rand(max(LI2,[],'all'), 3);
        RGB2  = label2rgb3d(LI2,CMAP2,[0,0,0]);
        volshow(RGB2);
    end

    % Save:
    if SAVEPROG
        Part2.Image          = TI;
        Part2.WatershedImage = WI;
        Part2.Watershed      = W;
        Part2.BW             = BW2;
        Part2.LI             = LI2;
        save([filepath,'02_Segs\',savename],'Part1','Part2');
    end
    
    end

    %% PART 3: Obtain isosurfaces within the basins:

    if RUNSTEPS > 2
    disp('Surfaces...');

    % Load last step:
    load([filepath,'02_Segs\',savename]);

    % Within each watershed basin, get an isosurface and weighted centroid fit:
    Surfaces = struct('faces',[],'vertices',[]);
    Centers  = zeros(max(Part2.LI, [],'all'),3);
    
    for ww = 1:max(Part2.LI,[],'all')
        
        if mod(ww,50) == 0
            disp(['Surfacing basin: ',num2str(ww),'/',num2str(max(Part2.LI,[],'all'))]); % Progress Bar
        end
    
        % Calculate an isovalue, isolate the basin from the image:
        ix           = find(Part2.LI==ww); % this particular basin's pixels
        V            = Part2.Image(ix); % pixel values
        ix_iso       = ix;
        ix_iso(V==0) = []; % remove all zero value pixels from the list
        IsoValue     = median(Part2.Image(ix)); % take the median of the rest
        [xi,yi,zi]   = ind2sub(size(Part2.Image), ix); % grab the basin
        temp = zeros(size(Part2.Image));
        for jj = 1:length(xi)
            temp(xi(jj),yi(jj),zi(jj)) = Part2.Image(xi(jj),yi(jj),zi(jj));
        end
        Surfaces(ww,1) = isosurface(temp,IsoValue);
    
        % Get the weighted centroid:
        ri = [xi,yi,zi]; % voxel location list
        wc_list = double(V)/sum(V) .* ri;
        wc      = sum(wc_list);
        Centers(ww,:) = wc;

    end
    
    % Show:
    if FIGVIZ == 3
        figure('units','centimeters','position',[0,0,20,20]);
        hold on; box on; set(gca,'linewidth',1);
        for ii = 1:length(Surfaces)
            Color = rand(1,3);
            patch('faces',Surfaces(ii).faces,'vertices',Surfaces(ii).vertices,'edgecolor','none','facealpha',0.5,'facecolor',Color);
            plot3(Centers(ii,2), Centers(ii,1), Centers(ii,3), 'x','linewidth',1,'markersize',10,'color',Color);
        end
        view(3); axis equal;
        camlight; material dull; lighting gouraud;
        % print([filepath,savename],'-dpng','-r300');
    end
    
    % Save:
    if SAVEPROG
        Part3.Surfaces  = Surfaces;
        Part3.Centers   = Centers;
        save([filepath,'02_Segs\',savename],'Part1','Part2','Part3');
    end
    
    end

    %% PART 4: Fit ellipsoids:
    
    if RUNSTEPS > 3
    disp('Ellipsoids...');

    % Load last step:
    load([filepath,'02_Segs\',savename]);
    
    if FIGVIZ == 4
        figure; hold on; box on; set(gca,'linewidth',1);
    end
    
    % Trim surfaces down to only ones that can be fit by an ellipsoid:
    TempSurfaces = [];
    for ii = 1:length(Part3.Surfaces)
        % Fit ellipsoid:
        try
            x = Part3.Surfaces(ii).vertices;
            warning('off','all');
            [t, s, r, ~, ~] = ellipsoid_fit_new(x);
            warning('on','all');
            TempSurfaces = [TempSurfaces, Part3.Surfaces(ii)];
        catch
            TempSurfaces = TempSurfaces;
        end
    end
    
    % For each surface, fit an ellipsoid:
    ELL = struct('T',[],'R',[],'S',[]);
    for ii = 1:length(TempSurfaces)
        
        if mod(ii,50) == 0
            disp(['Fitting ellipsoid: ',num2str(ii)]); % Progress bar
        end
    
        if ~isempty(TempSurfaces(ii).faces)
    
            % Fit ellipsoid:
            v = TempSurfaces(ii).vertices;
            c = mean(v);
            r = (v-c);
            rm = r./vecnorm(r,2,2);
            x = r + 0.5*rm;
            x = x + c;
            warning('off','all');
            [t, s, r, ~, ~] = ellipsoid_fit_new(x);
            S = [s(1); s(2); s(2)]; % duplicate the second largest radius, since the cells are spheroids
            S = eye(3).*S; S = [S(1,:),0; S(2,:),0; S(3,:),0; 0,0,0,1]; 
            R = [r(1,:),0; r(2,:),0; r(3,:),0; 0,0,0,1];
            T = [1,0,0,t(1); 0,1,0,t(2); 0,0,1,t(3); 0,0,0,1];
            warning('on','all');
        
            % Send to structure:
            ELL(ii,1).T = T;
            ELL(ii,1).R = R;
            ELL(ii,1).S = S;
        
            % Show:
            if FIGVIZ == 4
                M = T*R*S;
                [x,y,z] = sphere(20);
                x = x(:); y = y(:); z = z(:); r = [x,y,z,ones(length(x),1)];
                rp = M * r';
                K = convhull(rp(1:3,:)');    
            
                Color = rand(1,3);
                patch('faces',TempSurfaces(ii).faces,'vertices',Surfaces(ii).vertices,'edgecolor','none','facealpha',1,'facecolor',Color);
                trisurf(K, rp(1,:), rp(2,:), rp(3,:),'edgecolor','none','facealpha',0.5,'facecolor',Color);
            end
        end
    end
    if FIGVIZ == 4    
        view(3); axis equal;
        camlight; material dull; lighting gouraud;
    end
    
    if SAVEPROG
        Part4.Surfaces   = TempSurfaces;
        Part4.Ellipsoids = ELL;
        save([filepath,'02_Segs/',savename],'Part1','Part2','Part3','Part4');
    end
    
    end

    %% Part 5: Find overlapping ellipsoids:
    %{
    % Load last step:
    load([filepath,'01_Segs\',savename]);
    ELL = Part4.Ellipsoids;
    
    % Identify strongly overlapping ellipsoids:
    VolInt = zeros(length(ELL));
    for ii = 1:length(ELL)
        
        disp(['Checking for intersections with ellipsoid: ',num2str(ii)]);
        
        for jj = 1:length(ELL)
            if jj ~= ii
                % First determine if the two cells overlap:
                % Check all vertices of ellipsoid 1 to see if any reside inside
                % ellipsoid 2
                [x,y,z] = sphere(10);
                x = x(:); y = y(:); z = z(:); r_test = [x,y,z,ones(size(x))]';
                M1 = ELL(ii).T * ELL(ii).R * ELL(ii).S;
                M2 = ELL(jj).T * ELL(jj).R * ELL(jj).S;
                r_prime = M2 \ (M1*r_test);
                r_prime_norm = vecnorm(r_prime(1:3,:));
                if any(r_prime_norm <= 1)
                    % there is an intersection
                    int_yn = 1;
                    
                    % find a point common to both ellipsoids:
                    ix = find(r_prime_norm<=1);
                    common_point = r_test(1:3,ix(1));
                else
                    % there is no intersection
                    int_yn = 0;
                end
        
                % % If there is an intersection, find the intersection volume:
                % if int_yn
                %     v1 = M1*r_test;
                %     v2 = M2*r_test;
                %     [v_int] = INTERSECTION(v1(1:3,:), v2(1:3,:), common_point);
                %     [~,VolInt(ii,jj)] = convhull(v_int);
                % end
    
                % Record the intersection to a matrix:
                VolInt(ii,jj) = int_yn;
                
            end
        end
    end
    
    if FIGVIZ == 5
        figure; imagesc(VolInt);
    end
    %}


end