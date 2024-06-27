clear; clc;
imovie = 0;
set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


R = 0.1291;
chanWidths = R./[0.2; 0.4; 0.6; 0.75];
speeds = [6000 7500;
          3000 3750;
          2000 2500;
          1600 2000];
% speeds = speeds/1.5;

if 1
count = 1;
for iw = 1 : numel(chanWidths)
    for is = 1 : 2
      w = chanWidths(iw);
      vmax = speeds(iw,is);
      Cks(count,1) = 2*vmax*R^3/w;
      Cns(count,1) = R/w;
      count = count + 1;
    end
end


figure(1);clf;hold on;

for iw = 1 : numel(chanWidths)
    for is = 1 : 2
    speed = speeds(iw,is);
    chanWidth = chanWidths(iw);

    runNew = ['./test_exAdv_diff625kNetJune8_poisRuns_speed' num2str(speed) '_width' num2str(chanWidth) '.bin'];
    % runNew = ['./poisTrueRuns_dt5e-06_speed' num2str(speed) '_width' num2str(chanWidth) '.bin'];

    [vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runNew);
    vesxN = vesxN(:,1:end);
    vesyN = vesyN(:,1:end);
    timeN = timeN(:,1:end);

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
yticks(([0;0.5]))
yticklabels({'20','25'})
box on
pause

end
for iw = 1 : numel(chanWidths)
  for is = 1 : 2
      speed = speeds(iw,is);
      chanWidth = chanWidths(iw);

      disp(['speed = ' num2str(speed), ' width = ' num2str(chanWidth)])

      runNew = ['./test_exAdv_diff625kNetJune8_poisRuns_speed' num2str(speed) '_width' num2str(chanWidth) '.bin'];
      % runNew = ['./poisTrueRuns_dt5e-06_speed' num2str(speed) '_width' num2str(chanWidth) '.bin'];
      
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

  end
end