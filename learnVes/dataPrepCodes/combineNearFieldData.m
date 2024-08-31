clear; clc;

nInstances = 156225;

for imode = 2 : 128
  tTime = tic;
  velOnLayersReal = zeros(256,3,nInstances);
  velOnLayersImag = zeros(256,3,nInstances);
  for ives = 1 : nInstances
    load(['./vesicleByvesicle/vesicleID_' num2str(ives) '.mat'])
    velOnLayersReal(:,1,ives) = selfVelModesReal(:,imode);
    velOnLayersReal(:,2:3,ives) = VelOnGridModesReal(:,:,imode);
    velOnLayersImag(:,1,ives) = selfVelModesImag(:,imode);
    velOnLayersImag(:,2:3,ives) = VelOnGridModesImag(:,:,imode);
  end
  
  fileName = ['./modeByMode/nearField_mode' num2str(imode) '.mat'];
  save(fileName, "velOnLayersReal", "velOnLayersImag","nInstances")
  disp(['Mode ' num2str(imode) ' is done'])
  tEnd = toc(tTime);
  disp(['It takes ' num2str(tEnd) ' seconds'])
end
    