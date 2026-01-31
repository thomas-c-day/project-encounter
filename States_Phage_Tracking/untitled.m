datafolder = 'C:\Users\daythoma\Documents\Testing_Masking_Richard\';
jpgs = dir([datafolder,'JPGs\Series_017\','*.jpg']);
mask = dir([datafolder,'Masks\','*.tif']);

for jj = 1:length(jpgs)
    img = imread([datafolder,'JPGS\Series_017\',jpgs(jj).name]);
    IMG(:,:,jj) = img;
end
MASK = tiffreadVolume([datafolder,'Masks\',mask(2).name]);

figure;
for tt = 1:size(MASK,3)
    BB = bwboundaries(MASK(:,:,tt));
    cc = regionprops(MASK(:,:,tt),'EquivDiameter');
    D(tt) = cc.EquivDiameter;
    clf;
    imshow(imadjust(IMG(:,:,tt),[0.05,.4]));
    hold on;
    plot(BB{1}(:,2), BB{1}(:,1), 'r-','linewidth',2);
    hold off;
    F(tt) = getframe(gcf);
    drawnow;
end

% create the video writer with 1 fps
writerObj = VideoWriter('series_017.avi');
writerObj.FrameRate = 15;
% set the seconds per image

% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

figure;
plot(D,'-','linewidth',2);
hold on; box on; set(gca,'linewidth',1);
xlabel('Time');
ylabel('Eq. Diameter');