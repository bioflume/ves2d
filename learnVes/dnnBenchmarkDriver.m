clear; clc;
addpath ../src/
addpath ../examples/
addpath ./shannets/
addpath ./shannets/ves_fft_models/

pathofDocument = fileparts(which('Net_ves_relax_midfat.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pathofDocument = fileparts(which('Net_ves_adv_fft.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pathofDocument = fileparts(which('ves_fft_mode2.pth'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pe = pyenv('Version', '/Users/gokberk/opt/anaconda3/envs/mattorch/bin/python');

% FLAGS
%-------------------------------------------------------------------------
prams.bgFlow = 'relax'; % 'shear','tayGreen','relax','parabolic'
prams.speed = 0; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
prams.chanWidth = 0;
iplot = 0;
exactFreq = 0;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
errTol = 1e-2;
maxDt = 1E-5; % dt = 1.28e-3,1e-3, 1.6e-4, 1e-5, 1e-6


prams.Th = 10*maxDt;

% prams.Th = 0.05; % time horizon
prams.N = 128*1; % num. points for true solve in DNN scheme
prams.nv = 1; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = maxDt; % time step size
prams.dtRelax = prams.dt;
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;
oc = curve;
Th = prams.Th; N = prams.N; nv = prams.nv; dt = prams.dt; 
bgFlow = prams.bgFlow; speed = prams.speed;

% net parameters
Nnet = 128; % num. points
%-------------------------------------------------------------------------
disp(['Flow: ' prams.bgFlow ', N = ' num2str(N) ', nv = ' num2str(nv) ...
    ', Th = ' num2str(Th)])
%-------------------------------------------------------------------------

% VESICLES and WALLS:
% -------------------------------------------------------------------------


X0 = oc.initConfig(N,'ellipse');

[~,~,len] = oc.geomProp(X0);
X0 = X0./len;
IA = 0;
cent = [0; 0];
X = zeros(size(X0));
X(1:N) = cos(IA) * X0(1:N) - ...
      sin(IA) * X0(N+1:2*N) + cent(1);
X(N+1:2*N) = sin(IA) * X0(1:N) +  ...
      cos(IA) * X0(N+1:2*N) + cent(2);
[~,area0,len0] = oc.geomProp(X);

% load('./trueEquilX.mat')
% X = Xinit;
% [~,area0,len0] = oc.geomProp(X);


% figure(1); clf;
% plot(X(1:end/2),X(end/2+1:end))
% axis equal
% pause
% -------------------------------------------------------------------------

solveType = 'DNN';

fileName = ['./output/benchmark_relax.bin'];
fid = fopen(fileName,'w');
output = [N;nv];
fwrite(fid,output,'double');
x = X(1:end/2,:); y = X(end/2+1:end,:);
output = [x(:); y(:)];
fwrite(fid,output,'double');
fclose(fid);

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnn = dnnToolsSingleVes(X,prams);

% LOAD NORMALIZATION PARAMETERS
load ./shannets/ves_fft_in_param.mat
load ./shannets/ves_fft_out_param.mat
dnn.torchAdvInNorm = in_param;
dnn.torchAdvOutNorm = out_param;

tt = dnn.tt; dnn.oc = oc; 
% -------------------------------------------------------------------------

% INITIALLY TAKE SMALL TIME STEPS WITH IMPLICIT SOLVER
% ------------------------------------------------------------------------
% necessary to find correct initial tension, density, rotlet and stokeslet
tt.dt = maxDt; sig = zeros(N,nv); eta = []; RS = [];
% for iter = 1 : 2
%   vesicle = capsules(X,[],[],prams.kappa,ones(nv,1),1); vesicle.setUpRate();
%   [X,sig,eta,RS] = tt.timeStepSimple(X,sig,[],[],ones(nv,1),[],vesicle);
% end
tt.dt = dt;
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

timePred_list = [];
timeStand_list = [];
timeDestand_list = [];
timeCorr_list = [];

exactCnter = 0; iExact = false;
op = dnn.tt.op;
while time(end) < prams.Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(time(it))])
  
  
  disp('Taking a step with DNNs...');  tStart = tic;    
  
  % [Xnew, timeStand, timePred, timeDestand] = dnn.relaxWTorchBenchmark(Xhist);
  
  tEx = tic;
  vesicle = capsules(Xhist, [], [], 1, 1, 0);
  % vesicle.setUpRate();
  G = op.stokesSLmatrix(vesicle);
  [Ben,Ten,Div] = vesicle.computeDerivs;
  M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
  rhs = Xhist;
  LHS = (eye(2*vesicle.N)-vesicle.kappa*dnn.dt*(-G*Ben+M*G*Ben));
  Xnew = LHS\rhs;
  timeExact = toc(tEx);

  tReparamI = tic;
  [Xnew,~] = oc.reparametrize(Xnew,[],6,20);
  tReparamO = toc(tReparamI);

  % AREA-LENGTH CORRECTION
  tCorrI = tic;
  [Xnew,ifail] = oc.correctAreaAndLength2(Xnew,area0,len0);
  tCorrO = toc(tCorrI);
  
  timeCorr = tCorrO + tReparamO;

  [xIntersect,~,~] = oc.selfintersect(Xnew);
  if ~isempty(xIntersect); disp('New vesicle shape is self-intersecting!!!'); end;
  
  [~,area,len] = oc.geomProp(Xnew);
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  
  timePred_list(it,1) = timeExact;
  % timeStand_list(it,1) = timeStand;
  % timeDestand_list(it,1) = timeDestand;
  timeCorr_list(it,1) = timeCorr;

  it = it + 1;
  Xhist = Xnew;
  errALPred(it) = max(errArea,errLen);
  time(it) = time(it-1) + prams.dt;  
  
  disp(['Error in area and length: ' num2str(max(errArea, errLen))])   
  disp('********************************************') 
  disp(' ')
  
  if rem(it,1) == 0
    writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);  
  end

end

% Save data to a mat-file:
writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(filename,X,sigma,time,ncountNN,ncountExact)
x = X(1:end/2,:);
y = X(end/2+1:end,:);
output = [time;ncountNN;ncountExact;x(:);y(:);sigma(:)];

fid = fopen(filename,'a');
fwrite(fid,output,'double');
fclose(fid);


end
