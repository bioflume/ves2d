clear; clc;
imovie = 1;



% runNew = './mixedNetsFourierParaData_poisRuns_speed750_width0.32275.bin';
runNew = './test_exAdv_diff625kNet_poisRuns_speed400_width0.17213.bin';
runTrue = './truePoisRuns/poisTrueRuns_speed400_width0.17213.bin';
% runNew = './secondPhaseRuns128modes/poisRuns_speed1500_width0.32275.bin';


[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runTrue);
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runNew);


cxT = []; cyT = [];
cxN = []; cyN = [];
nsteps = min(numel(timeT),numel(timeN));
for k = 1 : nsteps
cxT(k,1) = mean(vesxT(:,k));
cyT(k,1) = mean(vesyT(:,k));

cxN(k,1) = mean(vesxN(:,k));
cyN(k,1) = mean(vesyN(:,k));

% figure(1);clf;
% plot(cxT(1:k),cyT(1:k),'r','linewidth',2)
% hold on
% plot(cxN(1:k),cyN(1:k),'b','linewidth',2)
% axis square
% legend('True','New')
% xlabel('cx')
% ylabel('cy')
% grid on
% pause(0.1)
end

pause

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

xT = vesxT(:,k) - mean(vesxT(:,k)) + 0.0;
xN = vesxN(:,k) - mean(vesxN(:,k)) + 0.0;


yT = vesyT(:,k); yN = vesyN(:,k);

figure(1); clf;
plot([xT; xT(1)], [yT; yT(1)],'r','linewidth',2)
hold on
plot([xN; xN(1)], [yN; yN(1)],'b','linewidth',2)

plot(xT(1), yT(1), 'ro','markerfacecolor','r')
plot(xN(1), yN(1), 'bo','markerfacecolor','b')


plot(mean(xT), mean(yT), 'ro','markerfacecolor','r')
plot(mean(xN), mean(yN), 'bo','markerfacecolor','b')



axis equal

plot(linspace(-0.5,0.5,100)',zeros(100,1),'Color',[253 219 199]/255,'linewidth',2)

legend('True', 'NN','Orientation','horizontal','Location','north')
%legend('Supposed2be', 'NetEquil','Orientation','horizontal','Location','north')
legend boxoff

xlim([-0.5 0.5])
ylim([-0.5 0.5])

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
cd(folder);
filebrowser;
message = sprintf('Finished creating movie file\n      %s.\n\nDone with demo!', fullFileName);
uiwait(helpdlg(message));



end