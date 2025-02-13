clear; clc;

% load checkCollIC.mat
load ./ShanSims/crashingTGdata.mat

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


%% 
oc = curve;
dt = 1E-4;
Th = 250*dt;
N = 32; nv = 48;
prams.nv = 48;
prams.repStrength = 0;
prams.bgFlow = 'tayGreen'; % 'shear','tayGreen','relax','parabolic'
prams.speed = 200; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
prams.chanWidth = 2.5;
prams.Th = Th;
prams.N = 32; % num. points for true solve in DNN scheme
prams.Nfmm = 32;
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = dt; % time step size
prams.dtRelax = dt;
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;

X = [vesx(:,:,1);vesy(:,:,1)];
dnn = dnnToolsManyVesFree(X,prams);

% LOAD NORMALIZATION PARAMETERS
load ./shannets/32modes_Adv_NormParams_2024Nov.mat
dnn.torchAdvInNorm = in_param;
dnn.torchAdvOutNorm = out_param;    
% % LOAD NEAR-SINGULAR NORMALIZATION PARAMS
load ./shannets/nearInterp_32modes_in_param.mat
load ./shannets/nearInterp_32modes_out_param.mat
dnn.torchNearInNorm = in_param;
dnn.torchNearOutNorm = out_param;

load ./shannets/32modes_tensionAdv_NormParams_2024Nov.mat
dnn.torchTenAdvInNorm = in_param;
dnn.torchTenAdvOutNorm = out_param;

opNfmm = dnn.opNfmm;
tt = dnn.tt;
op = tt.op;
%% Now take time steps
tenOld = zeros(N,nv);
X = [vesx(:,:,1);vesy(:,:,1)];
[~,area0,len0] = oc.geomProp(X);
for it = 1 : 50
  % X = [vesx_coll(:,:,it);vesy_coll(:,:,it)];  
  vback = dnn.vinf(X); %+ vInfs(:,:,it);

  % build vesicle class at the current step
  vesicle = capsules(X,[],[],1,1,0);
  nv = vesicle.nv;
  N = vesicle.N;

  % Compute bending forces + old tension forces
  fBend = vesicle.tracJump(X,zeros(N,nv));
  fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
  tracJump = fBend+fTen;
  
  % Near-field velocity
  G = op.stokesSLmatrix(vesicle);
  
  [velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
    rotCentNear, scalingNear, sortIdxNear] = dnn.predictNearLayersOnce32modes(vesicle.X);

  farFieldtracJump = dnn.computeStokesInteractionsNet_Alternative(vesicle, tracJump, opNfmm, oc, ...
    velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
    rotCentNear, scalingNear, sortIdxNear);
  
  figure(3);clf;
  plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
  hold on
  quiver(X(1:end/2,:),X(end/2+1:end,:),farFieldtracJump(1:end/2,:),farFieldtracJump(end/2+1:end,:),0.5,'r')
  axis equal
  title('FarField for Tension Solve')
  pause 

  % tension solves
  vBackSolve = dnn.invTenMatOnVback32modes(X, vback + farFieldtracJump);

  selfBendSolve = dnn.invTenMatOnSelfBend(X);
  tenNew = -(vBackSolve + selfBendSolve);

  % tenNewFilt = oc.filterTension(tenNew,32,4);
  % disp('Tension: ')
  % disp(tenNew)
  % figure(5);clf;
  % plot(tenNew,'r','linewidth',2)
  % hold on
  % plot(tenNewFilt,'b','linewidth',2)
  % axis square
  % title('tension')
  % pause
  % tenNew = tenNewFilt;
  tenOld = tenNew;

  % update the traction jump calculation
  fTen = vesicle.tracJump(zeros(2*N,nv), tenNew); 
  tracJump = fBend + fTen;

  % compute velocity again
  farFieldtracJump = dnn.computeStokesInteractionsNet_Alternative(vesicle, tracJump, opNfmm, oc, ...
    velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
    rotCentNear, scalingNear, sortIdxNear);
  
  figure(3);
  quiver(X(1:end/2,:),X(end/2+1:end,:),farFieldtracJump(1:end/2,:),farFieldtracJump(end/2+1:end,:),0.5,'b')
  axis equal
  legend('Vesicle','Vesicle','FarField for tension','FarField for Xsolve')
  title('FarField for X Solve')
  pause 


  % Total background velocity
  vbackTotal = vback + farFieldtracJump;

  % advection solve
  Xadv = dnn.translateVinfwTorch32modes(X, vbackTotal);
  
  figure(4);clf;
  plot(Xadv(1:end/2,:),Xadv(end/2+1:end,:),'k','linewidth',2)
  axis equal
  title('Advection solve')
  hold on
  pause
  % relaxation
  Xnew = dnn.relaxWTorchNet(Xadv);    
  
  plot(Xnew(1:end/2,:),Xnew(end/2+1:end,:),'r','linewidth',2)
  legend('Advection solve','Advection solve','Relaxation solve','Relaxation solve')
  axis equal
  pause
  % 
  XnewO = Xnew;
  for iter = 1 : 5
    Xnew = oc.redistributeArcLength(Xnew);
  end
  Xnew = oc.alignCenterAngle(XnewO,Xnew);

  % AREA-LENGTH CORRECTION
  disp('Area-Length correction after relaxation step')
  [Xnew,ifail] = oc.correctAreaAndLength2(Xnew,area0,len0);
  if ifail; disp('Error in AL cannot be corrected!!!'); end;
  

  figure(1);clf;
  plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
  hold on
  plot(Xnew(1:end/2,:),Xnew(end/2+1:end,:),'r','linewidth',2)
  axis equal
  title(it)

  figure(2);clf;
  plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
  hold on
  quiver(X(1:end/2,:),X(end/2+1:end,:),vback(1:end/2,:),vback(end/2+1:end,:),0.5,'b')
  quiver(X(1:end/2,:),X(end/2+1:end,:),farFieldtracJump(1:end/2,:),farFieldtracJump(end/2+1:end,:),0.5,'r')
  axis equal
  title(it)
    
  X = Xnew;

  % figure(3);clf;
  % plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
  % hold on
  % quiver(X(1:end/2,:),X(end/2+1:end,:),farFieldtracJump(1:end/2,:),farFieldtracJump(end/2+1:end,:),'r')
  % axis equal
  % title(it)

  pause()

end


