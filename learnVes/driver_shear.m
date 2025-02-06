clear; clc;
dt = 1E-5;
Th = dt;
% cx = [-0.4; 0];
% cy = [0.05; 0];

% First test IC
cx = [-0.3; 0]; % true
% cx = [-0.32; 0]; % true
cy = [0.040; 0];

IA = [0; pi/2];

prams.repStrength = 5E+4;

% Second test IC
% cx = [-0.28; 0];
% cy = [0.05;0];
% IA = [0; pi/2];

iExactTension = 0;
iExactNear = 0;
iExact = 0; % exact relaxation
iIgnoreNear = 0;
iAdv = 3; % 1: exact, 3: network

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

pathofDocument = fileparts(which('2024Nov_downsample32_ves_advten_mode1.pth'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pathofDocument = fileparts(which('ves_adv_downsample_fft_2024Oct_mode2.pth'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end


pe = pyenv('Version', '/Users/gokberk/opt/anaconda3/envs/mattorch/bin/python');

% FLAGS
%-------------------------------------------------------------------------
prams.bgFlow = 'shear'; % 'shear','tayGreen','relax','parabolic'
prams.speed = 2000; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
iplot = 1;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
errTol = 1e-2;
maxDt = dt; % dt = 1.28e-3,1e-3, 1.6e-4, 1e-5, 1e-6
prams.Th = Th;

% prams.Th = 0.05; % time horizon
prams.N = 32; % num. points for true solve in DNN scheme
prams.Nfmm = 32;
prams.nv = 2; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = maxDt; % time step size
prams.dtRelax = prams.dt;
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;
prams.chanWidth = 0;
oc = curve;
Th = prams.Th; N = prams.N; nv = prams.nv; dt = prams.dt; 
bgFlow = prams.bgFlow; speed = prams.speed;

% net parameters
Nnet = 32; % num. points
%-------------------------------------------------------------------------
disp(['Flow: ' prams.bgFlow ', N = ' num2str(N) ', nv = ' num2str(nv) ...
    ', Th = ' num2str(Th)])
%-------------------------------------------------------------------------

% VESICLES and WALLS:
% -------------------------------------------------------------------------
% X0 = oc.initConfig(N,'ellipse');
% 
% [~,~,len] = oc.geomProp(X0);
% X0 = X0./len;
% X = zeros(2*N,2);
% for k = 1 : 2
% X(1:N,k) = cos(IA(k)) * X0(1:N) - ...
%       sin(IA(k)) * X0(N+1:2*N) + cx(k);
% X(N+1:2*N,k) = sin(IA(k)) * X0(1:N)  + ...
%       cos(IA(k)) * X0(N+1:2*N) + cy(k);
% end
% 
% XOrig = X;
% for it = 1 : 5
%   X = oc.redistributeArcLength(X);
% end
% X = oc.alignCenterAngle(XOrig,X);

load('finalShearXclose.mat')
X = Xf; %[Xf(1:end/2,1)-mean(Xf(1:end/2,1)); Xf(end/2+1:end,1)-mean(Xf(end/2+1:end,1))];
for it = 1 : 5
  X = oc.redistributeArcLength(X);
end
X = [interpft(X(1:end/2,:),32);interpft(X(end/2+1:end,:),32)];
[~,area0,len0] = oc.geomProp(X);

 
figure(1); clf;
plot(X(1:end/2,:),X(end/2+1:end,:),'k-o')
hold on
% plot(Xnew(1:end/2,:),Xnew(end/2+1:end,:),'r')
axis equal
pause(0.1)
% -------------------------------------------------------------------------

solveType = 'DNN';
% fileName = ['./output/test_shear_ignoreNearN64_diff625kNetJune8_dt' num2str(dt) '_speed' num2str(prams.speed) '.bin'];
% fileName = ['./output/128modes_shear_nearNetrelaxNetTenNetAdvNet_noFiltering_dt' num2str(dt) '_speed' num2str(prams.speed) '.bin'];
fileName = ['./output/CheckingShansNetN32.bin'];
% fileName = ['./output/entangled_shear_biem_diff625kNetJune8_dt' num2str(dt) '_speed' num2str(prams.speed) '.bin'];
% fileName = ['./output/N64_shearTrueRuns_dt' num2str(dt) '_speed' num2str(speed) '.bin'];
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

if prams.N == 128
load ./shannets/mergedAdv_NormParams.mat
dnn.torchAdvInNorm = in_param;
dnn.torchAdvOutNorm = out_param;
elseif prams.N == 32
load ./shannets/32modes_Adv_NormParams_2024Nov.mat
dnn.torchAdvInNorm = in_param;
dnn.torchAdvOutNorm = out_param;    
end

% % LOAD NEAR-SINGULAR NORMALIZATION PARAMS
if prams.N == 128
% load ./shannets/nearInterp_fft_in_param.mat
% load ./shannets/nearInterp_fft_out_param.mat
load ./shannets/nearInterp_128modes_disth_params.mat
elseif prams.N == 32
load ./shannets/nearInterp_32modes_in_param.mat
load ./shannets/nearInterp_32modes_out_param.mat
end
dnn.torchNearInNorm = in_param;
dnn.torchNearOutNorm = out_param;

if prams.N == 128
load ./shannets/tensionAdv_NormParams_2024Oct.mat
dnn.torchTenAdvInNorm = in_param;
dnn.torchTenAdvOutNorm = out_param;
elseif prams.N == 32
load ./shannets/32modes_tensionAdv_NormParams_2024Nov.mat
dnn.torchTenAdvInNorm = in_param;
dnn.torchTenAdvOutNorm = out_param;
end

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

save('./nets32modeCheckData/input2sim.mat','Xhist');

% TIME STEPPING
it = 1;
while time(end) < prams.Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(time(it))])
  
  
  tStart = tic;    
  [Xhist,sigStore] = dnn.DNNsolveTorchMany(Xhist,sigStore,area0,len0,iExactTension,iExactNear,iExact,iIgnoreNear,iAdv);


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


  if iplot
  figure(1);clf;
  hold on;
  x = [Xhist(1:end/2,:); Xhist(1,:)];
  y = [Xhist(1+end/2:end,:); Xhist(end/2+1,:)];
  plot(x,y,'r','linewidth',2)
  hold on
  plot(Xhist(1,:), Xhist(end/2+1,:),'o','markerfacecolor','r','markersize',8)
  xlim([-0.5 0.5])
  ylim([-0.5 0.5])
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
