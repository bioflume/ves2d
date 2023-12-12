clear; clc;
imovie = 0;
iplot = 1;
fileName = './newDataGenRuns/dataGenSpeed10K_IDId17.bin';
fid = fopen(fileName,'r');
val = fread(fid,'double');
fclose(fid);

% Some important numbers
dt = val(1);
nv = val(2);
N = val(3);
nvbdExt = val(4);
NbdExt = val(5);
nvbdInt = val(6);
NbdInt = val(7);

% Remove them and move on with rigid boundaries
val = val(8:end);
XwallsExt = val(1:2*NbdExt);
xxwallsExt = XwallsExt(1:end/2);
yywallsExt = XwallsExt(end/2+1:end);
val = val(2*NbdExt+1:end);
XwallsInt = val(1:2*nvbdInt*NbdInt);
XwallsInt = reshape(XwallsInt,[2*NbdInt nvbdInt]);
val = val(2*nvbdInt*NbdInt+1:end);

% Now read in time steps
X = []; sigma = []; time = [];
etaExt = []; etaInt = []; RS = [];
nEntry2Read = 1 + 3 * N * nv + 2 * NbdExt + 2 * nvbdInt * NbdInt + 3*nvbdInt;
ist = 1;
while(numel(val)>=nEntry2Read)
time(ist) = val(1);
Xst = reshape(val(2:2*N*nv+1),[2*N nv]);
X(:,:,ist) = Xst;

val = val(2*N*nv+2:end);
sigma(:,:,ist) = reshape(val(1:N*nv), [N nv]);

val = val(N*nv+1:end);
etaExt(:,ist) = val(1:2*NbdExt);
val = val(2*NbdExt+1:end);
etaInt(:,:,ist) = reshape(val(1:2*nvbdInt*NbdInt), [2*NbdInt nvbdInt]);
val = val(2*NbdInt*nvbdInt+1:end);
RS(:,:,ist) = reshape(val(1:3*nvbdInt+3), [3 nvbdInt+1]);
val = val(3*nvbdInt+4:end);

ist = ist + 1;
end

% plot data
if iplot

if imovie
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
end

for it = 1 : numel(time)
figure(1); clf;
if imovie; cla reset; end;
plot([XwallsExt(1:end/2); XwallsExt(1)], [XwallsExt(end/2+1:end); XwallsExt(end/2+1)],'k','linewidth',2)
hold on
axis equal
plot([XwallsInt(1:end/2,:); XwallsInt(1,:)], [XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k','linewidth',2)
hWalls = fill([XwallsInt(1:end/2,:); XwallsInt(1,:)],[XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k');
set(hWalls,'edgecolor','k')
Xst = X(:,:,it);
plot([Xst(1:end/2,:); Xst(1,:)], [Xst(end/2+1:end,:);Xst(end/2+1,:)], 'r', 'linewidth',2)
hVes = fill([Xst(1:end/2,:); Xst(1,:)],[Xst(end/2+1:end,:);Xst(end/2+1,:)],'r');
set(hVes,'edgecolor','r')
title(['Time = ' num2str(time(it))])
if imovie
drawnow;
myMovie(it) = getframe(gca);
else
pause(0.1)
end

end

if imovie
message = sprintf('Done creating movie\nDo you want to play it?');
button = questdlg(message, 'Continue?', 'Yes', 'No', 'Yes');
drawnow;	% Refresh screen to get rid of dialog box remnants.
close(hFigure);
if strcmpi(button, 'Yes')
	hFigure = figure;
	% Enlarge figure to full screen.
	% set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
	title('Playing the movie we created', 'FontSize', 15);
	% Get rid of extra set of axes that it makes for some reason.
	axis off;
	% Play the movie.
	movie(myMovie);
	close(hFigure);
end

promptMessage = sprintf('Do you want to save this movie to disk?');
titleBarCaption = 'Continue?';
button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
if strcmpi(button, 'yes')
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
else
	uiwait(helpdlg('Done with demo!'));
end
end

% % create the video writer with 1 fps
% writerObj = VideoWriter('DatGen1K_ID34.avi');
% writerObj.FrameRate = 10;
% 
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
% % convert the image to a frame
% frame = F(i) ;    
% writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);


end