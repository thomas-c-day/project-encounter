function [BW_TD, Bounds] = DIC_PROCESS(IMG_TD, IMG_SP, SEGTHRESH_AG, AREA_MIN, DILATION_SIZE, FIGVIZ)
    % Thomas C. Day
    % generate aggregate boundaries from a probabilities segmentation output
    % from Ilastik.
    
    disp('Segmenting BF image...');

    % Show label probabilities:
    if FIGVIZ==1
        figure; imagesc(IMG_SP(:,:,1)); axis equal; colorbar; title('Label 1 Probability');
        figure; imagesc(IMG_SP(:,:,2)); axis equal; colorbar; title('Label 2 Probability');
    end

    % Binarize:
    BW0 = IMG_SP(:,:,1) > SEGTHRESH_AG;

    % Remove small objects:
    BW_TD = bwareafilt(BW0, [AREA_MIN, Inf]);

    % Clear the boundary:
    BW_TD = imclearborder(BW_TD);

    % Dilate and fill holes:
    BW_TD = imdilate(BW_TD, strel('disk', DILATION_SIZE, 6));
    BW_TD = imfill(BW_TD,'holes');

    % Get and show boundaries:
    Bounds = bwboundaries(BW_TD, 4);

    if FIGVIZ==1
        figure; imshowpair(BW0, BW_TD); title('Orig vs. Processed');
        
        figure; 
        imshow(imadjust(IMG_TD,[0,0.3])); 
        hold on; box on; set(gca,'linewidth',1);
        for bb = 1:length(Bounds)
            b2 = sgolayfilt(Bounds{bb}, 3, 11);
            % plot(Bounds{bb}(:,2), Bounds{bb}(:,1),'r-','linewidth',1);
            plot(b2(:,2), b2(:,1), 'r-','linewidth',1);
        end
    end

end