clear; clc;
imovie = 0;
set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


R = 0.1291;
chanWidths = R ./ [0.2; 0.4; 0.6; 0.75];
speeds = [1500 3000 4500;
    750 1500 2250;
    500 1000 1500;
    400 800 1200];
% speeds = speeds/1.5;

if 1
count = 1;
for iw = 1 : numel(chanWidths)
    for is = 1 : 3
      w = chanWidths(iw);
      vmax = speeds(iw,is);
      Cks(count,1) = 2*vmax*R^3/w;
      Cns(count,1) = R/w;
      count = count + 1;
    end
end


figure(1);clf;hold on;

for iw = 1 : numel(chanWidths)
    for is = 1 : 3
    speed = speeds(iw,is);
    chanWidth = chanWidths(iw);

    runNew = ['./poisDNN_dt1e-05_oldNN_speed' num2str(speed) '_width' num2str(chanWidth) '.bin'];
    [vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runNew);
    vesxN = vesxN(:,1:end-100);
    vesyN = vesyN(:,1:end-100);
    timeN = timeN(:,1:end-100);

    cxNew = (iw-1)*0.5;
    cyNew = (is-1)*0.5;

    if 1%iw > 1
    figure(1);
    plot([vesxN(:,end);vesxN(1,end)]-mean(vesxN(:,end))+cxNew,[vesyN(:,end);vesyN(1,end)]-mean(vesyN(:,end))+cyNew,'r','linewidth',2)
    xT = [vesxN(:,end);vesxN(1,end)]-mean(vesxN(:,end))+cxNew;
    yT = [vesyN(:,end);vesyN(1,end)]-mean(vesyN(:,end))+cyNew;
    hFill = fill(xT, yT,'r');
    set(hFill,'EdgeColor','r')
    disp(['iw = ' num2str(iw) ' is = ' num2str(is)])
    pause
    end
    end
end
axis equal
xticks(([0:numel(chanWidths)-1]*0.5))
xticklabels({'0.2','0.4','0.6','0.75'})
yticks(([0;0.5;1]))
yticklabels({'5','10','15'})
box on
pause

end
for iw = 1 : numel(chanWidths)
  for is = 1 : 3
      speed = speeds(iw,is);
      chanWidth = chanWidths(iw);

      disp(['speed = ' num2str(speed), ' width = ' num2str(chanWidth)])

      runNew = ['./poisDNN_dt1e-05_oldNN_speed' num2str(speed) '_width' num2str(chanWidth) '.bin'];

      
      [vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runNew);
    
      numFrames = numel(timeN);

    
      cxN = []; cyN = [];
      for k = 1 : numFrames
        cxN(k,1) = mean(vesxN(:,k));
        cyN(k,1) = mean(vesyN(:,k));
      end

      figure(2);clf;
      plot(timeN,cyN,'b','linewidth',2)
      axis square
      xlabel('Time')
      ylabel('cy')
      title(['R/W = ' num2str(R/chanWidth) ' Speed = ' num2str(speed)])
      grid on

      pause

      if 0 %imovie
      numberOfFrames = numFrames;
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
      
      if 0
      for k = 1 : 10 : numFrames


      if imovie; cla reset; end;

      xN = vesxN(:,k) - mean(vesxN(:,k)) + 0.0;
      yN = vesyN(:,k);

      figure(1); clf;
      plot([xN; xN(1)], [yN; yN(1)],'b','linewidth',2)
      hold on
      plot(xN(1), yN(1), 'bo','markerfacecolor','b')
    
      plot(mean(xN), mean(yN), 'bo','markerfacecolor','b')

      axis equal

      plot(linspace(-0.25,0.25,100)',zeros(100,1),'Color',[253 219 199]/255,'linewidth',2)

      

      % xlim([-0.25 0.25])
      % ylim([-chanWidth chanWidth])
      title(timeN(k))
      if imovie
      drawnow;
      myMovie(k) = getframe(gca);
      else
      pause(0.1)
      end

      end
      ax = gca;
      exportgraphics(ax,['~/Desktop/poisTrueRuns_speed' num2str(speed) '_width' num2str(chanWidth) '.png'],'Resolution',300)
      end

      if 0 %imovie
      startingFolder = pwd;
      defaultFileName = '*.avi';{'*.avi';'*.mp4';'*.mj2'}; %fullfile(startingFolder, '*.avi');
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
      end
  end
end