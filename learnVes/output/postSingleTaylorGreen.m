clear; clc;
imovie = 0;


fileName2 = 'taylorGreen_IC4_exactRelax2_predictNear_diff625kNetJune8_dt1e-05_speed500.bin';
fileName = 'taylorGreen_IC4_exactRelax_predictNear_diff625kNetJune8_dt1e-05_speed500.bin';
[vesx, vesy, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);
[vesx2, vesy2, ten, time2, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName2);
nsteps = numel(time); 

numberOfFrames = nsteps;
if imovie

hFigure = figure;
allTheFrames = cell(numberOfFrames,1);
vidHeight = 344;
vidWidth = 446;
allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
% Next get a cell array with all the colormaps.
allTheColorMaps = cell(numberOfFrames,1);
allTheColorMaps(:) = {zeros(256, 3)};
% Now combine these to make the array of structures.
myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);
% Create a VideoWriter object to write the video out to a new, different file.
% writerObj = VideoWriter('problem_3.avi');
% open(writerObj);
% Need to change from the default renderer to zbuffer to get it to work right.
% openGL doesn't work and Painters is way too slow.
set(gcf, 'renderer', 'zbuffer');

end

    



for k = 1 : 1 : numberOfFrames
 if imovie; cla reset; end;

  
 xvecN = [vesx(:,:,k);vesx(1,:,k)];
 yvecN = [vesy(:,:,k);vesy(1,:,k)];

 xvec2 = [vesx2(:,:,k);vesx2(1,:,k)];
 yvec2 = [vesy2(:,:,k);vesy2(1,:,k)];

 figure(1); clf; 
 
 % h2 = plot(xvecN, yvecN, 'Color',[215/255 25/255 28/255 1],'linewidth',2);
 h2 = plot(xvecN, yvecN, '-o','Color',[26/255 150/255 65/255 1],'linewidth',2);

 hold on

 h2 = plot(xvec2, yvec2,'--d', 'Color',[215/255 25/255 28/255 1],'linewidth',2);

 plot(-0.15, -0.15, 'k.','markersize',0.001)
 plot(1.5, 1.5, 'k.','markersize',0.001)

 axis equal

 xlim([-0.15 1.5])
 ylim([-0.15 1.5])


 set(gca,'xtick',[]);
 set(gca,'ytick',[]);
 set(gca,'ztick',[]);
    
 set(gca,'xcolor','w');
 set(gca,'ycolor','w');
 set(gca,'zcolor','w');
 box on
 set(gca,'visible','off')

 pause

 if imovie
 drawnow;
 myMovie(k) = getframe(gca);
 else
 pause(0.1)
 end

end


if imovie


% Get the name of the file that the user wants to save.
% Note, if you're saving an image you can use imsave() instead of uiputfile().
startingFolder = pwd;
defaultFileName = {'*.avi';'*.mp4';'*.mj2'}; %fullfile(startingFolder, '*.avi');
[baseFileName, folder] = uiputfile(defaultFileName, 'Specify a file');
if baseFileName == 0
	% User clicked the Cancel button.
	return;
end
fullFileName = fullfile(folder, baseFileName);
% Create a video writer object with that file name.
% The VideoWriter object must have a profile input argument, otherwise you get jpg.
% Determine the format the user specified:
[folder, baseFileName, ext] = fileparts(fullFileName);
switch lower(ext)
	case '.jp2'
		profile = 'Archival';
	case '.mp4'
		profile = 'MPEG-4';
	otherwise
		% Either avi or some other invalid extension.
		profile = 'Uncompressed AVI';
end
writerObj = VideoWriter(fullFileName, profile);
open(writerObj);
% Write out all the frames.
numberOfFrames = length(myMovie);
for frameNumber = 1 : numberOfFrames 
   writeVideo(writerObj, myMovie(frameNumber));
end
close(writerObj);
% Display the current folder panel so they can see their newly created file.
message = sprintf('Finished creating movie file\n      %s.\n\nDone with demo!', fullFileName);
uiwait(helpdlg(message));



end