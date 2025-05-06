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



% load checkCollIC.mat
load crashTGtestICs.mat
% load randomSamplesTG

for iter = 1 : 21
X = Xics(:,:,iter);
% figure(1); clf;
% plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
% axis equal
% for k = 1 : 32
% figure(1);clf;
% plot(X(1:end/2,k),X(end/2+1:end,k),'k','linewidth',2)
% title(k)
% axis equal
% pause
% end


%% 
oc = curve;
dt = 1E-5;
Th = 5*dt;
N = 32; nv = 32;
prams.nv = nv;
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

% X = [vesx(:,:,1);vesy(:,:,1)];
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
uprate = 1; %ceil(sqrt(N));

%% Now take time steps
tenOld = zeros(N,nv);
% X = [vesx(:,:,1);vesy(:,:,1)];
[~,area0,len0] = oc.geomProp(X);
cmap = colorcube(nv);

% for it = 1 : 1

  xnear = zeros(N,2,nv); ynear = zeros(N,2,nv);
  for k = 1 : nv
    [~,tang] = oc.diffProp(X(:,k));
    nx = tang(N+1:2*N);
    ny = -tang(1:N);
    xnear(:,1,k) = X(1:end/2,k) - nx*len0(k)/N;
    xnear(:,2,k) = X(1:end/2,k) + nx*len0(k)/N;
    ynear(:,1,k) = X(end/2+1:end,k) - ny*len0(k)/N;
    ynear(:,2,k) = X(end/2+1:end,k) + ny*len0(k)/N;
  end

  % X = [vesx_coll(:,:,it);vesy_coll(:,:,it)];  
  vback = dnn.vinf(X); %+ vInfs(:,:,it);

  % build vesicle class at the current step
  vesicle = capsules(X,[],[],1,1,0);
  nv = vesicle.nv;
  N = vesicle.N;

  % Compute bending forces + old tension forces
  fBend = vesicle.tracJump_upsample(X,zeros(N,nv),uprate);
  fTen = vesicle.tracJump_upsample(zeros(2*N,nv),tenOld,uprate);
  
  tracJump = fBend+fTen;
  
  tracJumpBeforeTenSolve = tracJump;

  % Near-field velocity
  % [velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
  %   rotCentNear, scalingNear, sortIdxNear] = dnn.predictNear5LayersOnce32modes(vesicle.X);
  % 
  % farFieldtracJumpNet = dnn.computeStokesInteractionsNet_5Layers(vesicle, tracJump, opNfmm, oc, ...
  %   velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
  %   rotCentNear, scalingNear, sortIdxNear);

  
  SLPnoCorr = []; % it would SLPdiag if fmm is on
  G = op.stokesSLmatrix(vesicle);

  % Get the near structure (this will be done using NN in Python)
  NearV2V = vesicle.getZone([],1);    

  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;

  SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
  farFieldtracJumpTrue = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,NearV2V,...
      kernel,kernelDirect,vesicle,true,false);

  farFieldtracJumpForTenSolve = farFieldtracJumpTrue;

  % figure(4);clf;
  % for k = 1 : nv
  % plot(X(1:end/2,k),X(end/2+1:end,k),'Color',cmap(k,:),'linewidth',2)
  % hold on
  % plot(xnear(:,:,k),ynear(:,:,k),'--','Color',cmap(k,:),'linewidth',2)
  % quiver(X(1:end/2,k),X(end/2+1:end,k),farFieldtracJumpNet(1:end/2,k),farFieldtracJumpNet(end/2+1:end,k),0.5,'Color',cmap(k,:));
  % quiver(X(1:end/2,k),X(end/2+1:end,k),farFieldtracJumpTrue(1:end/2,k),farFieldtracJumpTrue(end/2+1:end,k),0.5,'c');
  % end
  % hold on
  % axis equal
  % title('FarField for Tension Solve')

  % figure(5);clf;
  % for k = 1 : nv
  % plot(X(1:end/2,k),X(end/2+1:end,k),'Color',cmap(k,:),'linewidth',2)
  % hold on
  % plot(xnear(:,:,k),ynear(:,:,k),'--','Color',cmap(k,:),'linewidth',2)
  % quiver(X(1:end/2,k),X(end/2+1:end,k),fTen(1:end/2,k),fTen(end/2+1:end,k),0.5,'Color',cmap(k,:))
  % end
  % hold on
  % axis equal
  % title('tension force for Tension Solve')

  % figure(6);clf;
  % for k = 1 : nv
  % plot(X(1:end/2,k),X(end/2+1:end,k),'Color',cmap(k,:),'linewidth',2)
  % hold on
  % plot(xnear(:,:,k),ynear(:,:,k),'--','Color',cmap(k,:),'linewidth',2)
  % quiver(X(1:end/2,k),X(end/2+1:end,k),fBend(1:end/2,k),fBend(end/2+1:end,k),0.5,'Color',cmap(k,:))
  % end
  % hold on
  % axis equal
  % title('bending force for Tension Solve')

  % pause 
  %%
  % tension solves
  tenTrue = zeros(N,nv);
  trueVbackSolve = zeros(N,nv);
  trueSelfBendSolve = zeros(N,nv);

  G = op.stokesSLmatrix(vesicle);
  [~,Ten,Div] = vesicle.computeDerivs;
  for k = 1 : nv
    LHS = (Div(:,:,k)*G(:,:,k)*Ten(:,:,k));
    selfBend = G(:,:,k)*fBend(:,k);
    RHS = -Div(:,:,k)*(vback(:,k)+farFieldtracJumpTrue(:,k)+selfBend);
    trueSelfBendSolve(:,k) = LHS\(Div(:,:,k)*selfBend);
    trueVbackSolve(:,k) = LHS\(Div(:,:,k)*(vback(:,k)+farFieldtracJumpTrue(:,k)));
    tenTrue(:,k) = LHS\RHS;
  end % k = 1 : nv
  
  tension = tenTrue;
  % 
  % vBackSolve = dnn.invTenMatOnVback32modes(X, vback + farFieldtracJumpNet);
  % % 
  % selfBendSolve = dnn.invTenMatOnSelfBend(X);
  % tenNet = -(vBackSolve + selfBendSolve);

  % pause
  % tenNet = tenTrue;
  
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

  %%
  % tenOld = tenNew;

  % update the traction jump calculation
  fTenTrue = vesicle.tracJump_upsample(zeros(2*N,nv),tenTrue,uprate);
  tracJumpTrue = fBend + fTenTrue;
  % 
  % fTenNet = vesicle.tracJump_upsample(zeros(2*N,nv),tenNet,uprate);
  % tracJumpNet = fBend + fTenNet;
  % 
  % figure(5);clf;
  % for k = 1 : nv
  % plot(X(1:end/2,k),X(end/2+1:end,k),'Color',cmap(k,:),'linewidth',2)
  % hold on
  % quiver(X(1:end/2,k),X(end/2+1:end,k),tracJumpNet(1:end/2,k),tracJumpNet(end/2+1:end,k),0.5,'Color',cmap(k,:))
  % end
  % hold on
  % axis equal
  % title('Bending + Tension Force (Network)')

  % figure(6);clf;
  % for k = 1 : nv
  % plot(X(1:end/2,k),X(end/2+1:end,k),'Color',cmap(k,:),'linewidth',2)
  % hold on
  % quiver(X(1:end/2,k),X(end/2+1:end,k),tracJumpTrue(1:end/2,k),tracJumpTrue(end/2+1:end,k),0.5,'Color',cmap(k,:))
  % end
  % hold on
  % axis equal
  % title('Bending + Tension Force (True)')
  % pause
  %%
  % % compute velocity again
  % farFieldtracJumpNet = dnn.computeStokesInteractionsNet_5Layers(vesicle, tracJumpNet, opNfmm, oc, ...
  %   velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, transNear, rotateNear, ...
  %   rotCentNear, scalingNear, sortIdxNear);

  SLPnoCorr = []; % it would SLPdiag if fmm is on
  G = op.stokesSLmatrix(vesicle);

  % Get the near structure (this will be done using NN in Python)
  NearV2V = vesicle.getZone([],1);    

  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;

  SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
  farFieldtracJumpTrue = op.nearSingInt(vesicle,tracJumpTrue,SLP,SLPnoCorr,NearV2V,...
      kernel,kernelDirect,vesicle,true,false);
  

  farFieldtracJumpForXsolve = farFieldtracJumpTrue;

  vbackTotalTrue = vback + farFieldtracJumpTrue;
  % vbackTotalNet = vback + farFieldtracJumpNet;

  % figure(4);clf;
  % for k = 1 : nv
  % plot(X(1:end/2,k),X(end/2+1:end,k),'Color',cmap(k,:),'linewidth',2)
  % hold on
  % plot(xnear(:,:,k),ynear(:,:,k),'--','Color',cmap(k,:),'linewidth',2)
  % quiver(X(1:end/2,k),X(end/2+1:end,k),vbackTotalNet(1:end/2,k),vbackTotalNet(end/2+1:end,k),0.5,'Color',cmap(k,:));
  % quiver(X(1:end/2,k),X(end/2+1:end,k),vbackTotalTrue(1:end/2,k),vbackTotalTrue(end/2+1:end,k),0.5,'c');
  % end
  % hold on
  % axis equal
  % title('FarField+vback for X Solve')

  % figure(5);clf;
  % for k = 1 : nv
  % plot(X(1:end/2,k),X(end/2+1:end,k),'Color',cmap(k,:),'linewidth',2)
  % hold on
  % plot(xnear(:,:,k),ynear(:,:,k),'--','Color',cmap(k,:),'linewidth',2)
  % quiver(X(1:end/2,k),X(end/2+1:end,k),fTen(1:end/2,k),fTen(end/2+1:end,k),0.5,'Color',cmap(k,:))
  % end
  % hold on
  % axis equal
  % title('tension force for X Solve')
  % 
  % figure(6);clf;
  % for k = 1 : nv
  % plot(X(1:end/2,k),X(end/2+1:end,k),'Color',cmap(k,:),'linewidth',2)
  % hold on
  % plot(xnear(:,:,k),ynear(:,:,k),'--','Color',cmap(k,:),'linewidth',2)
  % quiver(X(1:end/2,k),X(end/2+1:end,k),fBend(1:end/2,k),fBend(end/2+1:end,k),0.5,'Color',cmap(k,:))
  % end
  % hold on
  % axis equal
  % title('bending force for X Solve')

  % pause 

  %%
  % Total background velocity
  vbackTotalTrue = vback + farFieldtracJumpTrue;
  % vbackTotalNet = vback + farFieldtracJumpNet;

  % % advection solve
  % XadvNet = dnn.translateVinfwTorch32modes(X, vbackTotalNet);
  % 
  XadvTrue = zeros(2*N,nv);
  G = op.stokesSLmatrix(vesicle);
  [Ben,Ten,Div] = vesicle.computeDerivs;
  for k = 1 : nv
    M = G(:,:,k)*Ten(:,:,k)*((Div(:,:,k)*G(:,:,k)*Ten(:,:,k))\eye(vesicle.N))*Div(:,:,k);
    XadvTrue(:,k) = X(:,k) + dt*(eye(2*vesicle.N)-M)*vbackTotalTrue(:,k);
  end
  % Xadv = oc.upsThenFilterShape(Xadv,N,8);

  % figure(5);clf;
  % plot(XadvTrue(1:end/2,:),XadvTrue(end/2+1:end,:),'k','linewidth',2)
  % axis equal
  % title('True Advection solve')
  % hold on
  % 
  % figure(6);clf;
  % plot(XadvNet(1:end/2,:),XadvNet(end/2+1:end,:),'k','linewidth',2)
  % axis equal
  % title('Ves-Net Advection solve')
  % hold on
  %%
  % pause
  % relaxation
  % XnewNet = dnn.relaxWTorchNet(XadvNet);    

  XnewTrue = zeros(2*N,nv);
  % vesicle = capsules(XadvTrue,[],[],1,1,0); 
  G = op.stokesSLmatrix(vesicle);
  % Bending, tension and surface divergence
  [Ben,Ten,Div] = vesicle.computeDerivs;
  for k = 1 : nv
  M = G(:,:,k)*Ten(:,:,k)*((Div(:,:,k)*G(:,:,k)*Ten(:,:,k))\eye(vesicle.N))*Div(:,:,k);
  rhs = XadvTrue(:,k);
  LHS = (eye(2*vesicle.N)-vesicle.kappa*dt*(-G(:,:,k)*Ben(:,:,k)+M*G(:,:,k)*Ben(:,:,k)));
  XnewTrue(:,k) = LHS\rhs;
  end
  % 
  % figure(5);
  % plot(XnewTrue(1:end/2,:),XnewTrue(end/2+1:end,:),'r','linewidth',2)
  % axis equal
  % 
  % figure(6);
  % plot(XnewNet(1:end/2,:),XnewNet(end/2+1:end,:),'r','linewidth',2)
  % axis equal
  % %%
  % % 
  % 
  % figure(7); clf;
  % plot(XnewTrue(1:end/2,:),XnewTrue(end/2+1:end,:),'k','linewidth',2)
  % hold on
  % plot(XnewNet(1:end/2,:),XnewNet(end/2+1:end,:),'r','linewidth',2)
  % axis equal
  % title('After a time step')
  % % pause
  
  fname = ['~/Desktop/comparisonX' num2str(iter) '.mat'];
  save(fname,'X', 'tracJumpBeforeTenSolve' , 'farFieldtracJumpForTenSolve', 'tension', 'trueVbackSolve','trueSelfBendSolve','farFieldtracJumpForXsolve', 'XadvTrue', 'XnewTrue')

  % XnewO = Xnew;
  % for iter = 1 : 5
  %   Xnew = oc.redistributeArcLength(Xnew);
  % end
  % Xnew = oc.alignCenterAngle(XnewO,Xnew);
  % 
  % % AREA-LENGTH CORRECTION
  % disp('Area-Length correction after relaxation step')
  % [Xnew,ifail] = oc.correctAreaAndLength2(Xnew,area0,len0);
  % if ifail; disp('Error in AL cannot be corrected!!!'); end;
  % 

  % figure(1);clf;
  % plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
  % hold on
  % plot(Xnew(1:end/2,:),Xnew(end/2+1:end,:),'r','linewidth',2)
  % axis equal
  % title(it)
  % 
  % figure(2);clf;
  % plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
  % hold on
  % quiver(X(1:end/2,:),X(end/2+1:end,:),vback(1:end/2,:),vback(end/2+1:end,:),0.5,'b')
  % quiver(X(1:end/2,:),X(end/2+1:end,:),farFieldtracJump(1:end/2,:),farFieldtracJump(end/2+1:end,:),0.5,'r')
  % axis equal
  % title(it)
    
  % X = Xnew;

  % figure(3);clf;
  % plot(X(1:end/2,:),X(end/2+1:end,:),'k','linewidth',2)
  % hold on
  % quiver(X(1:end/2,:),X(end/2+1:end,:),farFieldtracJump(1:end/2,:),farFieldtracJump(end/2+1:end,:),'r')
  % axis equal
  % title(it)

  % pause()

% end
end

