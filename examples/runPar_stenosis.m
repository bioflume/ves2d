clear; 
clc;

VC = 1; % 1; 5; 10; 50
Vinf = 5;
RA = 0.65;
kappa = 1e-1;
tube_length = 90;
tube_height = 6;
min_heights = (0.5:0.5:4)'; % 8 runs in parallel
contraction_width = 2; % 2; 3; 6; 12
iRunID = 88;

min_heights = min_heights(1:3);
VCs = [1;5;10];
%VCs = [10];
%[mm,vv] = meshgrid(min_heights,VCs);

[vv, mm] = meshgrid(VCs,min_heights);
runIDs = [1;2;3;9;10;11;17;18;19];

parfor i = 1 : numel(mm(:))
minH = mm(i);
VC = vv(i);
disp(i)
%runID = iRunID+i;
runID = runIDs(i);
runName = ['runID' num2str(runID)];
function_stenosis_runs(runName,VC,Vinf,RA,kappa,tube_length,tube_height,minH,contraction_width)
end
