clear; 
clc;

VC = 1; % 1; 5; 10; 50
Vinf = 5;
RA = 0.65;
kappa = 1e-1;
tube_length = 90;
tube_height = 6;
min_heights = (0.5:0.5:5)'; % 10 runs in parallel
contraction_width = 3; % 2; 3; 6; 12
iRunID = 0;

parfor i = 1 : numel(min_heights)
minH = min_heights(i);
disp(i)
runID = iRunID+i;
runName = ['runID' num2str(runID)];
function_stenosis_runs(runName,VC,Vinf,RA,kappa,tube_length,tube_height,minH,contraction_width)
end
