clear; clc;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

% load shanSim_matlab_version.mat
addpath ../ShanSims/

% filename = 'TG_N32_nv48.bin';

filename = 'BIEM_N32_nv128_VF25_TG_4layers.bin';
[vesxB, vesyB, time, N, nv, xinit, yinit] = loadShanVesFile(filename);

filename = 'ML_N32_nv128_VF25_TG_job232268.bin';
[vesxM, vesyM, time, N, nv, xinit, yinit] = loadShanVesFile(filename);

filename = 'ML_N32_nv128_VF25_TG_repul_job222818.bin';
[vesxMR, vesyMR, time, N, nv, xinit, yinit] = loadShanVesFile(filename);

nskip = 5;
nsteps = numel(time(1:nskip:32000));

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
%%
cnt = 1;
for k = 1 : nskip : 32000



xvecB = [vesxB(:,:,k); vesxB(1,:,k)];
yvecB = [vesyB(:,:,k); vesyB(1,:,k)];

xvecMR = [vesxMR(:,:,k); vesxMR(1,:,k)];
yvecMR = [vesyMR(:,:,k); vesyMR(1,:,k)];

xvecM = [vesxM(:,:,k); vesxM(1,:,k)];
yvecM = [vesyM(:,:,k); vesyM(1,:,k)];

figure(1); clf;
subplot(1,3,1)

plot(xvecB, yvecB, 'k','linewidth',2)
hold on
box on

axis equal
xlim([-0.5 5.5])
ylim([-0.5 5.5])

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

subplot(1,3,2)
plot(xvecM, yvecM, 'Color',[26/255 150/255 65/255],'linewidth',2)
hold on
box on

axis equal
xlim([-0.5 5.5])
ylim([-0.5 5.5])

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

subplot(1,3,3)
plot(xvecMR, yvecMR, 'Color',[202/255 0 32/255],'linewidth',2)
hold on
box on

axis equal
xlim([-0.5 5.5])
ylim([-0.5 5.5])

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);


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
