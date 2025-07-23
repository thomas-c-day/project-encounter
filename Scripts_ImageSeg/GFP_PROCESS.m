function [Xpk, Ipk, wiGFP1] = GFP_PROCESS(IMG_GFP, SEGTHRESH_GFP, EXP_BEAD_DIAM, FIGVIZ)
    % Thomas C. Day
    % Segment a GFP file for microbeads of a certain size.
    
    disp('Segmenting GFP image...');
    
    % Binarize:
    BW_GFP = IMG_GFP > SEGTHRESH_GFP;
    % BW_GFP = imdilate(BW_GFP, strel('disk',EXP_BEAD_DIAM,6));

    % Make a weighted image:
    wiGFP1 = double(IMG_GFP) .* double(BW_GFP);  % mask noise
    diGFP  = bwdist(~BW_GFP, 'euclidean');       % distance transform to try to isolate peaks
    wiGFP  = diGFP .* wiGFP1;                    % combine distance and brightness to determine "true peaks"

    % Watershed the BW image, get centroids:
    W = watershed(-wiGFP);
    BW_GFP(W==0) = 0;
    M = regionprops(BW_GFP,'Centroid');
    for mm = 1:length(M)
        Xpk(1,mm) = M(mm).Centroid(1);
        Xpk(2,mm) = M(mm).Centroid(2);
    end
    Ppk = round(Xpk'); 
    Ppk = sub2ind(size(IMG_GFP), Ppk(:,2), Ppk(:,1));
    Ipk = zeros(size(BW_GFP));
    Ipk(Ppk) = 1;

    % Old method that I decided I didn't like:
    % Get peaks, with filter for size:
    % sizefilter = fspecial('gaussian', 2*EXP_BEAD_DIAM, 1*EXP_BEAD_DIAM);
    % [Xpk, Ipk] = FastPeakFind(wiGFP, 0, sizefilter);

    % Show:    
    if FIGVIZ==2
        
        % Show pixel histogram:
        figure; histogram(IMG_GFP); set(gca,'yscale','log');

        % Show binarization mask:
        figure;
        imshowpair(imadjust(IMG_GFP), BW_GFP,'falsecolor');

        % Show raw and identified particles:
        figure;
        imagesc(imadjust(IMG_GFP,[0,0.2])); % increase the contrast
        colormap('hot');
        hold on;
        plot(Xpk(1,:), Xpk(2,:), 'b.','markersize',10);
        axis equal;
        title('(Adjusted) Raw image and ID beads');

        % Show overlay of Ipk and raw image:
        figure;
        imshowpair(imadjust(IMG_GFP,[0,0.2]), Ipk, 'falsecolor');
        
    end
end