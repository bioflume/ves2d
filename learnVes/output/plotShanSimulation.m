clear; clc;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

% load shanSim_matlab_version.mat
addpath ../ShanSims/

% filename = 'TG_N32_nv48.bin';
filename = 'TG_dilute_auglag_25Feb.bin';
[vesx, vesy, time, N, nv, xinit, yinit] = loadShanVesFile(filename);

nskip = 5;
nsteps = numel(time(10000:nskip:20000));

% numberOfFrames = nsteps;
% hFigure = figure;
% allTheFrames = cell(numberOfFrames,1);
% vidHeight = 344;
% vidWidth = 446;
% allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
% allTheColorMaps = cell(numberOfFrames,1);
% allTheColorMaps(:) = {zeros(256, 3)};
% myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);
% set(gcf, 'renderer', 'zbuffer');

cnt = 1;
for k = 10000 : nskip : 20000

figure(1); clf;


xvec = [vesx(:,:,k); vesx(1,:,k)];
yvec = [vesy(:,:,k); vesy(1,:,k)];

plot(xvec, yvec, 'Color',[202,0,32]/255,'linewidth',2)
hold on
% hFill = fill(xvec, yvec, [202,0,32]/255);
% set(hFill,'EdgeColor', [202,0,32]/255);



% plot(xvec, yvec, 'Color',[5,113,176]/255,'linewidth',2)
% hold on
% hFill = fill(xvec, yvec, [5,113,176]/255);
% hFill.FaceAlpha = 0.5;
% set(hFill,'EdgeColor', [5,113,176]/255);

box on

axis equal
xlim([-0.5 3])
ylim([-0.5 3])

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

title(k)

drawnow;
myMovie(cnt) = getframe(gca);
cnt = cnt + 1;

pause(0.1)
end

startingFolder = pwd;
fullFileName = '~/Desktop/newRun.avi'; 

[folder, baseFileName, ext] = fileparts(fullFileName);
switch lower(ext)
	case '.jp2'
		profile = 'Archival';
	case '.mp4'
		profile = 'MPEG-4';
    otherwise
		profile = 'Uncompressed AVI';
end
writerObj = VideoWriter(fullFileName, profile);
open(writerObj);

numberOfFrames = length(myMovie);
for frameNumber = 1 : numberOfFrames 
   writeVideo(writerObj, myMovie(frameNumber));
end
close(writerObj);
message = sprintf('Finished creating movie file\n      %s.\n\nDone with demo!', fullFileName);


%%
