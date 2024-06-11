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
% chanWidths = [chanWidths(1);chanWidths(4)];

dt = 5e-5;
for iw = 1 : numel(chanWidths)
  for is = 1 : 3
      %driver_mixedRuns(chanWidths(iw),speeds(iw,is),timeHorizons(iw,is));
      % driver_torchRuns(chanWidths(iw),speeds(iw,is),timeHorizons(iw,is));
      driver_trueRuns(chanWidths(iw),speeds(iw,is),timeHorizons(iw,is),dt);
  end
end



%%
R = 0.1291;
speeds = [800 4500 2250 1500 1200];
chanWidths = R./[0.75; 0.2; 0.4; 0.6; 0.75];
timeHorizons = [0.025; 0.15; 0.05; 0.02; 0.015];

for irun = 1 : 5
driver_torchRuns(chanWidths(irun),speeds(irun),timeHorizons(irun));
end


