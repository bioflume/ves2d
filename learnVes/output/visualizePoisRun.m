clear all;
clc;

[vesxN, vesyN, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile('poisDNNnewSingVes_speed12000_newNet_trueAdv_noSplit.bin');

imovie = 1;
XN = [vesxN;vesyN];

numberOfFrames = numel(time);
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

count = 1;
for it = 1 : 10 : numel(time)
figure(1); clf;
if imovie; cla reset; end;
Xst = XN(:,it);
Xst(1:end/2) = Xst(1:end/2) - mean(Xst(1:end/2));
plot([Xst(1:end/2); Xst(1)], [Xst(end/2+1:end);Xst(end/2+1)], 'r', 'linewidth',2)
hold on
plot(Xst(1), Xst(end/2+1),'ko','markersize',10,'MarkerFaceColor','k')
plot(mean(Xst(1:end/2)),mean(Xst(end/2+1:end)),'ko','markerfacecolor','k')
axis equal


plot(linspace(-0.25,0.25,100)',zeros(100,1),'Color',[253 219 199]/255,'linewidth',2)
xlim([-0.25 0.25])
ylim([-0.3 0.3])

title('Dt = 1E-6')
% title(['Time = ' num2str(time(it))])
if imovie
drawnow;
myMovie(count) = getframe(gca);
count = count + 1;
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
cd(folder);
filebrowser;
message = sprintf('Finished creating movie file\n      %s.\n\nDone with demo!', fullFileName);
uiwait(helpdlg(message));



end