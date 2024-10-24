clear; clc;
imovie = 1;



% fileName = '128modes_taylorGreen_IC5_GT50ves_dt1e-05_speed200.bin';
% [vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);


% load 32modes_TaylorGreen_50Ves_NearNet % nearnet simulation

% load 128modes_TaylorGreen_50Ves_BIEM % Ground truth

load 32modes_TaylorGreen_50Ves_BIEM % Low-Res BIEM

nsteps = numel(timeT);


if imovie

numberOfFrames = numel(1:5:nsteps);
hFigure = figure;
allTheFrames = cell(numberOfFrames,1);
vidHeight = 576;
vidWidth = 1024;
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

Vsize = 2.5;
[xx,yy] = meshgrid(linspace(-1,3.5,50)',linspace(-1,3.5,50)');
uu = sin(xx/Vsize*pi).*cos(yy/Vsize*pi); 
vv = -cos(xx/Vsize*pi).*sin(yy/Vsize*pi);


frameCount = 1;
for k = 1 : 5 : nsteps
 if imovie; cla reset; end;

 % kT = 2*k-1;
 kT = k;
 xvecT = [vesxT(:,:,kT);vesxT(1,:,kT)] ;
 yvecT = [vesyT(:,:,kT);vesyT(1,:,kT)];
  

 figure(1); clf; 
 h = plot(xvecT, yvecT, 'Color',[26/255 150/255 65/255 1],'linewidth',2);
 % h = plot(xvecT, yvecT, 'Color',[202/255 0 32/255 1],'linewidth',2);
 
 % h = plot(xvecT, yvecT, 'Color',[0 0 0 1],'linewidth',2);
 for j = 1 : 9
 set(h(j),'Color',[h(j).Color, 1],'linewidth',2)
 end
 hold on
 l = streamslice(xx,yy,uu,vv);
 set(l,'Color',[12/255,44/255,132/255, 0.75])
 set(l,'linewidth',0.5)
 
 plot(-0.5, -0.5, 'k.','markersize',0.001)
 plot(3, 3, 'k.','markersize',0.001)

 axis equal

 xlim([-0.5 3])
 ylim([-0.5 3])


 set(gca,'xtick',[]);
 set(gca,'ytick',[]);
 set(gca,'ztick',[]);
    
 set(gca,'xcolor','w');
 set(gca,'ycolor','w');
 set(gca,'zcolor','w');
 box on
 set(gca,'visible','off')

 
 

 if imovie
 drawnow;
 myMovie(frameCount) = getframe(gca);
 frameCount = frameCount + 1;
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
% cd(folder);
% filebrowser;
message = sprintf('Finished creating movie file\n      %s.\n\nDone with demo!', fullFileName);
uiwait(helpdlg(message));



end