clear; clc;
dt = 1E-5;
Th = 0.03;

iExactTension = 1;
iExactNear = 0;
iExact = 1; % exact relaxation
iIgnoreNear = 0;

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
prams.bgFlow = 'tayGreen'; % 'shear','tayGreen','relax','parabolic'
prams.speed = 500; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
iplot = 0;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
errTol = 1e-2;
maxDt = dt; % dt = 1.28e-3,1e-3, 1.6e-4, 1e-5, 1e-6
prams.Th = Th;

% prams.Th = 0.05; % time horizon
prams.N = 128; % num. points for true solve in DNN scheme
prams.Nfmm = 128;
prams.nv = 9; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
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
[~,area0,len0] = oc.geomProp(X0);
scale = 1/len0;

sx = [.17 -.00 -.07 .12 .04 -.12 .15 .02 -.17]*scale*2; 
sy = [.15 -.010 +.08 .04 .02 -.01 -.01 .06 -.12]*scale*2; 
cenx = kron(ones(1,3),[pi/4 pi/2 3*pi/4]*scale*1.85);
ceny = kron([pi/4 pi/2 3*pi/4]*scale*1.85,ones(1,3));
% scale = 0.225;
cenx = cenx + sx;
ceny = ceny + sy;

angle = -ones(prams.nv,1);

X = oc.initConfig(prams.N,'nv',prams.nv,...
  'reducedArea',0.65,...
  'angle',angle,...
  'center',[cenx;ceny], 'scale',scale);

XOrig = X;
for it = 1 : 5
  X = oc.redistributeArcLength(X);
end
X = oc.alignCenterAngle(XOrig,X);

% load tayGreenStep140ic
% X = Xic;
% load ./output/taylorGreenFinalIC3_nearNet
% load ./output/taylorGreenFinalIC3_trueFiner
% load taylorGreenFinalIC3_trueLowRes
% X = Xic;

prams.chanWidth = pi*scale*2;
[~,area0,len0] = oc.geomProp(X);
% 

% figure(1); clf;
% plot(X(1:end/2,:),X(end/2+1:end,:),'k')
% hold on
% % plot(Xnew(1:end/2,:),Xnew(end/2+1:end,:),'r')
% axis equal
% pause
% -------------------------------------------------------------------------

solveType = 'DNN';
% fileName = ['./output/taylorGreen_IC4_ignoreNear_diff625kNetJune8_dt' num2str(dt) '_speed' num2str(prams.speed) '.bin'];
% fileName = ['./output/taylorGreen_IC4_true_diff625kNetJune8_dt' num2str(dt) '_speed' num2str(prams.speed) '.bin'];
% fileName = ['./output/taylorGreen_IC4_nearNet_diff625kNetJune8_dt' num2str(dt) '_speed' num2str(prams.speed) '.bin'];
% fileName = ['./output/taylorGreen_IC4_exactNear_diff625kNetJune8_dt' num2str(dt) '_speed' num2str(prams.speed) '.bin'];
fileName = ['./output/taylorGreen_IC4_exactRelax2_predictNear_diff625kNetJune8_dt' num2str(dt) '_speed' num2str(prams.speed) '.bin'];

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

% LOAD NORMALIZATION PARAMETERS
load ./shannets/ves_fft_in_param.mat
load ./shannets/ves_fft_out_param.mat
dnn.torchAdvInNorm = in_param;
dnn.torchAdvOutNorm = out_param;

% % LOAD NEAR-SINGULAR NORMALIZATION PARAMS
load ./shannets/nearInterp_fft_in_param.mat
load ./shannets/nearInterp_fft_out_param.mat
dnn.torchNearInNorm = in_param;
dnn.torchNearOutNorm = out_param;

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
while time(end) < prams.Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(time(it))])
  
  
  disp('Taking a step with DNNs...');  tStart = tic;    
  [Xhist,sigStore] = dnn.DNNsolveTorchMany(Xhist,sigStore,area0,len0,iExactTension,iExactNear,iExact,iIgnoreNear);


  [xIntersect,~,~] = oc.selfintersect(Xhist);
  if ~isempty(xIntersect); disp('New vesicle shape is self-intersecting!!!'); break; end;
  
  [~,area,len] = oc.geomProp(Xhist);
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  
  it = it + 1;
  time(it) = time(it-1) + prams.dt;  
  
  disp(['Error in area and length: ' num2str(max(errArea, errLen))])   
  disp('********************************************') 
  disp(' ')
  
  if rem(it,1) == 0
    writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);  
  end
  toc(tStart)

  if iplot
  figure(1);clf;
  hold on;
  x = [Xhist(1:end/2,:); Xhist(1,:)];
  y = [Xhist(1+end/2:end,:); Xhist(end/2+1,:)];
  plot(x,y,'r','linewidth',2)
  hold on
  plot(Xhist(1,:), Xhist(end/2+1,:),'o','markerfacecolor','r','markersize',8)
  xlim([0 pi*2*0.224])
  ylim([0 pi*2*0.224])
  axis equal
  pause(0.1)
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
