clear; clc;

addpath ../../src/
addpath ../
addpath ../../examples/
addpath ../shannets/
addpath ../shannets/ves_fft_models/

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
prams.bgFlow = 'shear'; % 'shear','tayGreen','relax','parabolic'
prams.speed = 2000; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
prams.Th = 1;

% prams.Th = 0.05; % time horizon
prams.N = 128; % num. points for true solve in DNN scheme
prams.Nfmm = 128;
prams.nv = 2; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = 1E-5; % time step size
prams.dtRelax = prams.dt;
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;
prams.chanWidth = 0;


load shanSim_matlab_version.mat
tsteps = [1;2;105;106;150;151];
istep = 3;

% These are the initial X and tension with which the rest is computed
% X = reshape(X(tsteps(istep),:,:),256,2);
shan_tenNew = reshape(ten(tsteps(istep)+1,:,:),128,2);
ten = reshape(ten(tsteps(istep),:,:),128,2);

load ~/Desktop/tensionCheck.mat

% Build capsules and curve
oc = curve;
vesicle = capsules(X,[],[],1,1,0);

% tension term -- applyes curve.arcDeriv
tenTerm = vesicle.tensionTerm(ten);

% bending term
benTerm = vesicle.bendingTerm(X);

tracJump = benTerm + tenTerm;


%%
opNfmm = poten(128);
dnn = dnnToolsManyVesFree(X,prams);
load ../shannets/nearInterp_128modes_disth_params.mat
dnn.torchNearInNorm = in_param;
dnn.torchNearOutNorm = out_param;

load ../shannets/tensionAdv_NormParams_2024Oct.mat
dnn.torchTenAdvInNorm = in_param;
dnn.torchTenAdvOutNorm = out_param;

[velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
    rotCentNear, scalingNear, sortIdxNear] = dnn.predictNearLayersOnceAllModes(vesicle.X);    

farFieldtracJump = dnn.computeStokesInteractionsNet_Alternative(vesicle, tracJump, opNfmm, oc, ...
    velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
    rotCentNear, scalingNear, sortIdxNear);

%%
% farFieldtracJump = oc.upsThenFilterShape(farFieldtracJump,4*128,16);
% farFieldtracJump1(tsteps(istep)+1,:,1) = oc.upsThenFilterShape(farFieldtracJump1(tsteps(istep)+1,:,1)',4*128,16);
% farFieldtracJump1(tsteps(istep)+1,:,2) = oc.upsThenFilterShape(farFieldtracJump1(tsteps(istep)+1,:,2)',4*128,16);


%%
vback = 2000*[X(end/2+1:end,:);zeros(size(X(1:end/2,:)))];
vBackSolve = dnn.invTenMatOnVback(X, vback + farFieldtracJump);
selfBendSolve = dnn.invTenMatOnSelfBend(X);
tenNew = -(vBackSolve + selfBendSolve);
tracJumpNew = benTerm + vesicle.tensionTerm(tenNew);

%%

tenNewExact = zeros(128,2);
selfBendExact = zeros(128,2);
vBackSolveExact = zeros(128,2);
G = opNfmm.stokesSLmatrix(vesicle);
[~,Ten,Div] = vesicle.computeDerivs;
for k = 1 : 2
LHS = (Div(:,:,k)*G(:,:,k)*Ten(:,:,k));
selfBend = G(:,:,k)*benTerm(:,k);
RHS = -Div(:,:,k)*(vback(:,k)+farFieldtracJump(:,k)+selfBend);
selfBendExact(:,k) = LHS\(-Div(:,:,k)*selfBend);
vBackSolveExact(:,k) = LHS\(-Div(:,:,k)*(vback(:,k)+farFieldtracJump(:,k)));
tenNewExact(:,k) = LHS\RHS;
end % k = 1 : nv

  %%

figure(1);clf;
plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
hold on
quiver(X(1:end/2,1),X(end/2+1:end,1),farFieldtracJump(1:end/2,1),farFieldtracJump(end/2+1:end,1),0.7)
quiver(X(1:end/2,2),X(end/2+1:end,2),farFieldtracJump(1:end/2,2),farFieldtracJump(end/2+1:end,2),0.7)
quiver(X(1:end/2,1),X(end/2+1:end,1),farFieldtracJump1(tsteps(istep)+1,1:end/2,1)',farFieldtracJump1(tsteps(istep)+1,end/2+1:end,1)',0.7)
quiver(X(1:end/2,2),X(end/2+1:end,2),farFieldtracJump1(tsteps(istep)+1,1:end/2,2)',farFieldtracJump1(tsteps(istep)+1,end/2+1:end,2)',0.7)
axis equal
legend('vesicle','vesicle','matlab','matlatb','python','python')


figure(2);clf;
plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
hold on
quiver(X(1:end/2,1),X(end/2+1:end,1),tracJumpNew(1:end/2,1),tracJumpNew(end/2+1:end,1),0.7)
quiver(X(1:end/2,2),X(end/2+1:end,2),tracJumpNew(1:end/2,2),tracJumpNew(end/2+1:end,2),0.7)
quiver(X(1:end/2,1),X(end/2+1:end,1),tracJump2(tsteps(istep)+1,1:end/2,1)',tracJump2(tsteps(istep)+1,end/2+1:end,1)',0.7)
quiver(X(1:end/2,2),X(end/2+1:end,2),tracJump2(tsteps(istep)+1,1:end/2,2)',tracJump2(tsteps(istep)+1,end/2+1:end,2)',0.7)
axis equal
legend('vesicle','vesicle','matlab','matlatb','python','python')

%%
XnewRelax = dnn.relaxWTorchNet(X);