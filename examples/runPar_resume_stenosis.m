clear; 
clc;

VC = 1; % 1; 5; 10; 50
Vinf = 5;
RA = 0.65;
kappa = 1e-1;
tube_length = 90;
tube_height = 6;
min_heights = (0.5:0.5:4)'; % 8 runs in parallel
contraction_width = 3; % 2; 3; 6; 12
iRunID = 40;

VCs = [1;5;10];
minH = min_heights(3);
runIDs = [3; 11; 19];
parfor i = 1 : numel(VCs)
VC = VCs(i);
disp(i)
runID = runIDs(i);
runName = ['runID' num2str(runID)];
icFileName = ['runID' num2str(runID) '_resFine.mat'];
function_stenosis_resume_runs(runName,VC,Vinf,RA,kappa,tube_length,tube_height,minH,contraction_width,icFileName);
end
