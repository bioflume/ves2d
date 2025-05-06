clear; clc;
dt = 1E-4;
Th = 2*dt;

iExactTension = 1;
iExactNear = 1;
iExact = 1; % exact relaxation
iIgnoreNear = 0;

addpath ../src/
addpath ../examples/
addpath ./denseTayGreenICs/

% FLAGS
%-------------------------------------------------------------------------
prams.bgFlow = 'tayGreen'; % 'shear','tayGreen','relax','parabolic'
prams.speed = 200; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
iplot = 0;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
errTol = 1e-2;
maxDt = dt; % dt = 1.28e-3,1e-3, 1.6e-4, 1e-5, 1e-6
prams.Th = Th;

% prams.Th = 0.05; % time horizon
prams.N = 32; % num. points for true solve in DNN scheme
prams.Nfmm = 32;
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = maxDt; % time step size
prams.dtRelax = prams.dt;
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;
Th = prams.Th; N = prams.N; dt = prams.dt; 

oc = curve;

% net parameters
Nnet = 32; % num. points


fileNames = {'nv48IC.mat','nv102IC.mat','nv126IC.mat','nv260IC.mat',...
    'nv504IC.mat','nv846IC.mat','nv1020IC.mat','nv2250IC.mat'};
cpuTimes = zeros(7,1);

for irun = 8 : 8

load(['./denseTayGreenICs/' fileNames{irun}])
X = [interpft(X(1:end/2,:),prams.N);interpft(X(end/2+1:end,:),prams.N)];

prams.chanWidth = chanWidth;
[~,area0,len0] = oc.geomProp(X);
prams.nv = numel(X(1,:));
prams.repStrength = 0;
% -------------------------------------------------------------------------

nv = prams.nv; 
bgFlow = prams.bgFlow; speed = prams.speed;
%-------------------------------------------------------------------------
disp(['Flow: ' prams.bgFlow ', N = ' num2str(N) ', nv = ' num2str(nv) ...
    ', Th = ' num2str(Th)])
%-------------------------------------------------------------------------

solveType = 'DNN';

fileName = ['./output/' fileNames{irun} '_timingSim.bin'];


fid = fopen(fileName,'w');
output = [N;nv];
fwrite(fid,output,'double');
x = X(1:end/2,:); y = X(end/2+1:end,:);
output = [x(:); y(:)];
fwrite(fid,output,'double');
fclose(fid);

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnn = dnnToolsManyVesFree(X,prams);

tt = dnn.tt; dnn.oc = oc; 
% -------------------------------------------------------------------------

% INITIALLY TAKE SMALL TIME STEPS WITH IMPLICIT SOLVER
% ------------------------------------------------------------------------
% necessary to find correct initial tension, density, rotlet and stokeslet
tt.dt = maxDt; sig = zeros(N,nv); eta = []; RS = [];
% ------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
time = [0];
Xhist = X; sigStore = sig; 
errALPred = 0;
ncountCNN = 0;
ncountExct = 0;
driftyNet = [];
driftyAdv = [];
% ------------------------------------------------------------------------
writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);

% TIME STEPPING
it = 1;
tTimeSteps = tic;
while time(end) < prams.Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(time(it))])
  
  
  disp('Taking an exact time step...');  
  [Xhist,sigStore] = dnn.DNNsolveTorchMany(Xhist,sigStore,area0,len0,iExactTension,iExactNear,iExact,iIgnoreNear,1);
  
  it = it + 1;
  time(it) = time(it-1) + prams.dt;  

  disp('********************************************') 
  disp(' ')
end
cpuTimes(irun) = toc(tTimeSteps)/2;
% save cpuTimesForOneStepN32 cpuTimes
% Save data to a mat-file:
% writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(filename,X,sigma,time,ncountNN,ncountExact)
x = X(1:end/2,:);
y = X(end/2+1:end,:);
output = [time;ncountNN;ncountExact;x(:);y(:);sigma(:)];

fid = fopen(filename,'a');
fwrite(fid,output,'double');
fclose(fid);


end
