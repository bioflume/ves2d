R = 0.1291;
chanWidths = R ./ [0.2; 0.4; 0.6; 0.75];
speeds = [1500 3000 4500;
    750 1500 2250;
    500 1000 1500;
    400 800 1200];
timeHorizons = [0.4 0.2 0.15;
    0.1 0.08 0.05;
    0.05 0.04 0.02;
    0.04 0.025 0.015];

dt = 1e-5;
for iw = 1 : numel(chanWidths)
  for is = 1 : 3
      %driver_mixedRuns(chanWidths(iw),speeds(iw,is),timeHorizons(iw,is));
      %driver_torchRuns(chanWidths(iw),speeds(iw,is),timeHorizons(iw,is));
      % driver_trueRuns(chanWidths(iw),speeds(iw,is),timeHorizons(iw,is),dt);
      DNNsolve_oldNetSimple(chanWidths(iw), speeds(iw,is), timeHorizons(iw,is))
  end
end



%%
R = 0.1291;
speeds = [800 4500 2250 1500 1200 1500 3000];
chanWidths = R./[0.75; 0.2; 0.4; 0.6; 0.75; 0.2; 0.2];
timeHorizons = [0.025; 0.15; 0.05; 0.02; 0.015; 0.4; 0.2];

for irun = 1 : 7
% driver_torchRuns(chanWidths(irun),speeds(irun),timeHorizons(irun));
DNNsolve_oldNetSimple(chanWidths(irun), speeds(irun), timeHorizons(irun));
end

%% Higher Cas
dt = 5E-6;
R = 0.1291;
dtTorch = [5E-6 5E-6;
           5E-6 5E-6;
           2.5E-6 2.5E-6;
           2.5E-6 2.5E-6];

chanWidths = R./[0.2; 0.4; 0.6; 0.75];
speeds = [6000 7500;
          3000 3750;
          2000 2500;
          1600 2000];

timeHorizons = [0.04 0.05;
    0.020 0.015; % first two is additional
    0.015 0.015;
    0.010 0.008];

for iw = 1 : 2
  for is = 1 : 2
    driver_resume_torchRuns(chanWidths(iw),speeds(iw,is),timeHorizons(iw,is),dtTorch(iw,is));
    driver_resume_trueRuns(chanWidths(iw),speeds(iw,is),timeHorizons(iw,is),dtTorch(iw,is));
  end
end

for iw = 3 : 4 %numel(chanWidths)
  for is = 1 : 2
    driver_torchRuns(chanWidths(iw),speeds(iw,is),timeHorizons(iw,is),dtTorch(iw,is));
  end
end


%% 
R = 0.1291;
speeds = [7500 3000 3750];
chanWidths = R./[0.2; 0.4; 0.4];
timeHorizons = [0.05; 0.02; 0.015];
for irun = 2 : 3
driver_resume_torchRuns(chanWidths(irun),speeds(irun),timeHorizons(irun),5E-6);
end
