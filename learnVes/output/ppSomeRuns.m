clear; clc;
imovie = 1;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

runNew = './test_exAdv_diff625kNet_poisRuns_speed400_width0.17213.bin';
runTrue = './truePoisRuns/poisTrueRuns_speed400_width0.17213.bin';

%runOld = './poisDNNnewSingVesInterp5_oldNN_higher.bin';
runNew1 = './diffNetOldDataNewStandPhaseRuns/testNewStandRun3_noInitExact_diff160kNet_init100exact_poisRuns_speed400_width0.17213.bin';


[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runTrue);
%[vesxO, vesyO, ten, timeO, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runOld);
[vesxN1, vesyN1, ten, timeN1, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runNew1);
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runNew);


% cxT = mean(vesxT,1);
% size(cxT)
% cyT = mean(vesyT,1);
% cxF = mean(vesxN1,1);
% cyF = mean(vesyN1,1);
% cxR = mean(vesxN,1);
% cyR = mean(vesyN,1);
% cxT = reshape(cxT,numel(cxT),1);
% cyT = reshape(cyT,numel(cxT),1);
% cxF = reshape(cxF,numel(cxF),1);
% cyF = reshape(cyF,numel(cyF),1);
% cxR = reshape(cxR,numel(cxR),1);
% cyR = reshape(cyR,numel(cyR),1);
% plot(cxT, cyT,'Color',[94 60 153]/255,'linewidth',2)
% hold on
% plot(cxF,cyF,'Color',[253 184 99]/255,'linewidth',2)
% plot(cxR, cyR,'Color',[230 97 1]/255,'linewidth',2)
% xlim([0 100])
% xlabel('center-x')
% ylabel('center-y')
% grid
% axis square
% legend('True','FourierIn','FourierIn-DriftNet')
% legend boxoff

nsteps = min([numel(timeT); numel(timeN1); numel(timeN)]);
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

    



for k = 1 : 1 : nsteps
kN = k;
% if k == 1
%   kN = k;
% else
%  kN = (k-1)*10;
% end

if imovie; cla reset; end;

xT = vesxT(:,k) - mean(vesxT(:,k)) + 0.0;
% xO = vesxO(:,kN) - mean(vesxO(:,kN)) + 0.5;
xN1 = vesxN1(:,kN) - mean(vesxN1(:,kN)) + 0.5;
xN = vesxN(:,kN) - mean(vesxN(:,kN)) + 1.0;


yT = vesyT(:,k); 
% yO = vesyO(:,kN); 
yN = vesyN(:,kN);
yN1 = vesyN1(:,kN);

figure(1); clf;
plot([xT; xT(1)], [yT; yT(1)],'Color',[94 60 153]/255,'linewidth',2)
hold on
% plot([xO; xO(1)], [yO; yO(1)],'Color',[178 171 210]/255,'linewidth',2)
plot([xN1; xN1(1)], [yN1; yN1(1)],'Color',[253 184 99]/255,'linewidth',2)
plot([xN; xN(1)], [yN; yN(1)],'Color',[230 97 1]/255,'linewidth',2)




hFill = fill([xT; xT(1)], [yT; yT(1)],[94 60 153]/255);
set(hFill,'EdgeColor',[94 60 153]/255)
% hFill = fill([xO; xO(1)], [yO; yO(1)],[178 171 210]/255);
% set(hFill,'EdgeColor',[178 171 210]/255)
hFill = fill([xN1; xN1(1)], [yN1; yN1(1)],[253 184 99]/255);
set(hFill,'EdgeColor',[253 184 99]/255)
hFill = fill([xN; xN(1)], [yN; yN(1)],[230 97 1]/255);
set(hFill,'EdgeColor',[230 97 1]/255)

plot(xT(1), yT(1), 'ko','markerfacecolor','k')
% plot(xO(1), yO(1), 'ko','markerfacecolor','k')
plot(xN1(1), yN1(1), 'ko','markerfacecolor','k')
plot(xN(1), yN(1), 'ko','markerfacecolor','k')


plot(mean(xT), mean(yT), 'ko','markerfacecolor','k')
% plot(mean(xO), mean(yO), 'ko','markerfacecolor','k')
plot(mean(xN), mean(yN), 'ko','markerfacecolor','k')
plot(mean(xN1), mean(yN1), 'ko','markerfacecolor','k')



axis equal

plot(linspace(-0.25,1.75,100)',zeros(100,1),'Color',[0.5 0.5 0.5],'linewidth',2)

% legend('Ground Truth','PRE (2019)', 'New-biased','New','Orientation','horizontal','Location','north')
legend('True','160kOld','625kNew','Orientation','horizontal','Location','north')
legend boxoff

xlim([-0.25 1.25])
ylim([-0.4 0.5])

titleStr = ['t = ' num2str(timeT(k),'%.2f')];
title(titleStr,'FontSize',28)

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
cd(folder);
filebrowser;
message = sprintf('Finished creating movie file\n      %s.\n\nDone with demo!', fullFileName);
uiwait(helpdlg(message));



end