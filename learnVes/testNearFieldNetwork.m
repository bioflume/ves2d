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

% load finalShearX.mat
% X = Xhist;

load('~/Desktop/taylorClose4VesIC.mat')
X = Xic;

prams.bgFlow = 'tayGreen'; % 'shear','tayGreen','relax','parabolic'
prams.speed = 200; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
prams.Th = 0.75;

% prams.Th = 0.05; % time horizon
prams.N = 128; % num. points for true solve in DNN scheme
prams.nv = numel(X(1,:));
prams.Nfmm = 128;
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = 1E-5; % time step size
prams.dtRelax = 1E-5;
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;
Nnet = 128; % num. points
prams.chanWidth = 2.5;

dnn = dnnToolsManyVesFree(X,prams);
load ./shannets/nearInterp_128modes_disth_params.mat
dnn.torchNearInNorm = in_param;
dnn.torchNearOutNorm = out_param;

oc = curve;
nv = numel(X(1,:));
N = numel(X(:,1))/2;

vesicle = capsules(X,[],[],prams.kappa,1,0);
tracJump = vesicle.tracJump(X,zeros(N,nv));

[velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
    rotCentNear, scalingNear, sortIdxNear] = dnn.predictNearLayersOnceAllModes(vesicle.X);    

stokesInteract = dnn.computeStokesInteractionsNet_Alternative(vesicle, tracJump, dnn.opNfmm, oc, ...
    velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
    rotCentNear, scalingNear, sortIdxNear);


op = dnn.tt.op;
SLPnoCorr = []; % it would SLPdiag if fmm is on
G = op.stokesSLmatrix(vesicle);

% Get the near structure (this will be done using NN in Python)
NearV2V = vesicle.getZone([],1);    

kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;

SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
stokesInteractTrue = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,NearV2V,...
    kernel,kernelDirect,vesicle,true,false);

% figure(1);clf;
% plot(Xstand(1:end/2,k),Xstand(end/2+1:end,k),'linewidth',2)
% hold on
% plot(tracersX(1:end/2,:,k),tracersX(end/2+1:end,:,k),'k.','markersize',8)
% quiver(tracersX(1:end/2,:,k),tracersX(end/2+1:end,:,k),velx_stand,vely_stand)
% axis equal
% 
% figure(2); clf;
% plot(X(1:end/2,k),X(end/2+1:end,k),'linewidth',2)
% hold on
% plot(xlayers(:,:,k),ylayers(:,:,k),'k.','markersize',8)
% quiver(xlayers(:,:,k),ylayers(:,:,k),velx(:,:,k),vely(:,:,k))
% axis equal
% 
% title(k)
% pause




% From the tracJump build the velocity
% normalize the tracJump, it is a vector, so normalize as vinf
% then find the velocity on the layers, then denormalize both the velocity
% and the layers


% Rotate, translate tracersX to xlayers, ylayers -- ready for input form
% Also rotate velocity to velx and vely


