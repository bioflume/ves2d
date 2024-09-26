clear; clc;
imovie = 1;


% load taylorGreen_IC3_true_long_dt5E6_speed500

% vesxL = vesxT; vesyL = vesyT; timeL = timeT;

% load taylorGreen_IC3_nearNetLonger_diff625kNetJune8_long_dt1E5_speed500

% load taylorGreen_IC3_trueFiner_long_dt5E6_speed500

fileName = 'taylorGreen_IC4_true_diff625kNetJune8_dt1e-05_speed500.bin';
[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

fileName = '32modes32_taylorGreen_IC4_nearNet_diff625kNetJune8_dt1e-05_speed500.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

% load taylorGreen_IC4_nearNetLongest_dt1E5_speed500

nsteps = min(numel(timeT),numel(timeN));


if imovie

numberOfFrames = nsteps;
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

 % kT = 2*k-1;
 kT = k;
 xvecT = [vesxT(:,:,kT);vesxT(1,:,kT)] ;
 yvecT = [vesyT(:,:,kT);vesyT(1,:,kT)];
  
 xvecN = [vesxN(:,:,k);vesxN(1,:,k)];
 yvecN = [vesyN(:,:,k);vesyN(1,:,k)];

 % xvecL = [vesxL(:,:,k);vesxL(1,:,k)];
 % yvecL = [vesyL(:,:,k);vesyL(1,:,k)];

 figure(1); clf; 
 h = plot(xvecT, yvecT, 'Color',[0 0 0 0.75],'linewidth',2);
 for j = 1 : 9
 set(h(j),'Color',[h(j).Color, 0.75],'linewidth',2)
 end
 hold on

 h2 = plot(xvecN, yvecN, 'Color',[215/255 25/255 28/255 1],'linewidth',2);
 for j = 1 : 9
 set(h2(j),'Color',[h2(j).Color, 1],'linewidth',2)
 end
 hold on
  % h2 = plot(xvecL, yvecL, 'Color',[26/255 150/255 65/255 1],'linewidth',2);
  % for j = 1 : 9
  % set(h2(j),'Color',[h2(j).Color, 1],'linewidth',2)
  % end

 
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
% cd(folder);
% filebrowser;
message = sprintf('Finished creating movie file\n      %s.\n\nDone with demo!', fullFileName);
uiwait(helpdlg(message));



end