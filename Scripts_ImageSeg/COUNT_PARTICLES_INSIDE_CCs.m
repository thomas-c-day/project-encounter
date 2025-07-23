function [CC, AggregateArea, NumParticles, IntegratedFL] = COUNT_PARTICLES_INSIDE_CCs(BW_TD, IPK, IMG_GFP)

    % Thomas C. Day.
    % Count the beads within each connected component that defines an
    % aggregate.
    % Do this in two different ways:
    % 1. Count the number of peaks inside the connected component boundary.
    % 2. Integrate the total fluorescence (after thresholding) within the
    % boundary.

    disp('Counting particles within connected regions...');

    % Generate the connected components of the TD image:
    CC   = bwconncomp(BW_TD, 4); % get connected components from TD image
    PkIx = find(IPK); % obtain peak indices

    % Preallocate measurements:
    AggregateArea = zeros(CC.NumObjects,1);
    NumParticles  = AggregateArea;
    IntegratedFL  = AggregateArea;

    % Loop through each connected component, find attached particles:
    for ii = 1:CC.NumObjects
        % Progress bar:
        disp([num2str(ii),' / ',num2str(CC.NumObjects)]);

        % Version 1: No extra dilation:
        CCix = CC.PixelIdxList{ii};

        % Version 2: Small extra dilation:
        %{
        % Examine only object of interest:
        temp = zeros(CC.ImageSize);
        temp(CC.PixelIdxList{ii}) = 1;

        % Dilate a bit:
        temp = imdilate(temp, strel('disk',3,6));
        CCix = find(temp); % obtain new pixel indices
        %}

        % -----------------------------------------------
        % 1. Count number of peaks inside the CC boundary:
        % Find any instances in peak list:
        MemberTF = ismember(PkIx, CCix);
        MemberIx = find(MemberTF);

        % Record object properties:
        AggregateArea(ii) = length(CCix);
        NumParticles(ii)  = length(MemberIx);

        % -----------------------------------------------
        % 2. Integrate total fluorescence inside the boundary:
        GFP_vals = IMG_GFP(CCix);
        IntegratedFL(ii) = sum(GFP_vals);

    end



end