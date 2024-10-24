clear; clc;
imovie = 1;
set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


% fileNameNN = 'moreModes32_test_shear_interpNear_diff625kNetJune8_dt1e-05_speed2000.bin';
% fileNameNN = 'test_shear_ignoreNear_diff625kNetJune8_dt1e-05_speed2000.bin';
fileNameNN = '128modes_shear_nearNet_relaxNet_dt1e-05_speed2000.bin';
% fileNameTR = 'test_shear_ignoreNear_diff625kNetJune8_dt1e-05_speed2000.bin';

fileNameTR = 'N128again_shearTrueRuns_dt1e-05_speed2000.bin';

[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameNN);
[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);

% load testShearNearInterpSpeed2000 

nsteps = min([numel(timeN); numel(timeT)]);

cxTrue = mean(vesxT,1); cyTrue = mean(vesyT,1);
cxNN = mean(vesxN,1); cyNN = mean(vesyN,1);

if ~imovie
for k = 1 : 1 :nsteps
  figure(1); clf;
  plot([vesxN(:,:,k);vesxN(1,:,k)], [vesyN(:,:,k);vesyN(1,:,k)], 'r', 'linewidth', 2)
  hold on
  plot([vesxT(:,:,k);vesxT(1,:,k)], [vesyT(:,:,k);vesyT(1,:,k)], 'b', 'linewidth', 2)

  plot(vesxN(1,:,k),vesyN(1,:,k),'ro','markersize',8,'markerfacecolor','r')
  plot(vesxT(1,:,k),vesyT(1,:,k),'bo','markersize',8,'markerfacecolor','b')
  axis equal
  
  xlim([-0.5 1])
  ylim([-0.5 0.5])
  
  % legend('Network','Ignoring near-field')
  % legend boxoff
  title(['Time: ' num2str(timeN(k))])
   
  pause(0.1)

end
end

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


for k = 1 :  numberOfFrames
cla reset;

figure(1); clf;
xvec = [vesxT(:,1,k);vesxT(1,1,k)];
yvec = [vesyT(:,1,k);vesyT(1,1,k)];
plot(xvec, yvec, 'Color',[244,165,130]/255,'linewidth',2)
hold on
hFill = fill(xvec, yvec, [244,165,130]/255);
set(hFill,'EdgeColor', [244,165,130]/255);
plot(vesxT(1,1,k),vesyT(1,1,k),'ko','markersize',8,'markerfacecolor','k')

xvec = [vesxT(:,2,k);vesxT(1,2,k)];
yvec = [vesyT(:,2,k);vesyT(1,2,k)];
plot(xvec, yvec, 'Color',[146,197,222]/255,'linewidth',2)
hold on
hFill = fill(xvec, yvec, [146,197,222]/255);
hFill.FaceAlpha = 0.5;
set(hFill,'EdgeColor', [146,197,222]/255);
plot(vesxT(1,2,k),vesyT(1,2,k),'ko','markersize',8,'markerfacecolor','k')



% figure(1); clf;
xvec = [vesxN(:,1,k);vesxN(1,1,k)];
yvec = [vesyN(:,1,k);vesyN(1,1,k)];
plot(xvec, yvec, 'Color',[202,0,32]/255,'linewidth',2)
hold on
hFill = fill(xvec, yvec, [202,0,32]/255);
set(hFill,'EdgeColor', [202,0,32]/255);
plot(vesxN(1,1,k),vesyN(1,1,k),'ko','markersize',8,'markerfacecolor','k')

xvec = [vesxN(:,2,k);vesxN(1,2,k)];
yvec = [vesyN(:,2,k);vesyN(1,2,k)];
plot(xvec, yvec, 'Color',[5,113,176]/255,'linewidth',2)
hold on
hFill = fill(xvec, yvec, [5,113,176]/255);
hFill.FaceAlpha = 0.5;
set(hFill,'EdgeColor', [5,113,176]/255);
plot(vesxN(1,2,k),vesyN(1,2,k),'ko','markersize',8,'markerfacecolor','k')


set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);
        
set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')

axis equal
xlim([-0.75 1])
ylim([-0.35 0.35])

titleStr = ['t = ' num2str(timeT(k),'%.2f')];
title(titleStr,'FontSize',28)

drawnow;
myMovie(k) = getframe(gca);
end

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


end %iMovie
%%
figure(2);clf;
distTrue = sqrt((cxTrue(1,1,:)-cxTrue(1,2,:)).^2 + (cyTrue(1,1,:)-cyTrue(1,2,:)).^2);
distNN = sqrt((cxNN(1,1,:)-cxNN(1,2,:)).^2 + (cyNN(1,1,:)-cyNN(1,2,:)).^2);
distTrue = reshape(distTrue,numel(distTrue),1);
distNN = reshape(distNN,numel(distNN),1);
plot(timeT(1:nsteps),distTrue(1:nsteps),'r','linewidth',2)
hold on
plot(timeN(1:nsteps),distNN(1:nsteps),'b','linewidth',2)
axis square
grid on
legend('True','Network')
legend boxoff
xlabel('Time')
ylabel('Distance between vesicles')