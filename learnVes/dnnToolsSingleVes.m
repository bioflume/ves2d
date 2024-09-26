classdef dnnToolsSingleVes

properties
KbDts    
variableKbDt
tt
dt
dtRelax
vinf
bendNets
muChan_bend
sdevChan_bend
scale_bend
offset_bend
MVnets
muChan_MV
sdevChan_MV
scale_MV
offset_MV
MVoutSize
tenBendNets
muChan_tenBend
sdevChan_tenBend
scale_tenBend
offset_tenBend
nTenModes
tenPredModes
tenVnets
muChan_tenV
sdevChan_tenV
scale_tenV
offset_tenV
tenVoutSize
nTenVmodes
TenVactiveModes
nVelModes
velActiveModes
nComp
nCompRelax
colMeans
evects
oc
kappa
confined
interpOrder
torchAdvInNorm
torchAdvOutNorm
driftNet
driftNet_muChan
driftNet_sdevChan
driftNet_muOutput 
driftNet_sdevOutput
end

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = dnnToolsSingleVes(X,prams)
o.confined = false;
o.kappa = 1;
o.interpOrder = prams.interpOrder;
o.dt = prams.dt;    
o.tt = o.buildTstep(X,prams);  
o.vinf = o.setBgFlow(prams.bgFlow,prams.speed,prams.chanWidth);  
o.dtRelax = prams.dtRelax;
o.oc = curve;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,dyNet,dyAdv] = DNNsolveTorchNoSplit(o,Xold,area0,len0,iExact)
oc = o.oc;
vback = o.vinf(Xold);
advType = 1; % 1: exact, 2: with old net, 3: with Torch net
% 1) COMPUTE THE ACTION OF dt*(1-M) ON Xold  
if advType == 1
  vback = o.vinf(Xold);

  tt = o.tt;
  op = tt.op;
  vesicle = capsules(Xold,[],[],1,1,0);
 
  G = op.stokesSLmatrix(vesicle);
  [~,Ten,Div] = vesicle.computeDerivs;

  M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
  Xadv = Xold + o.dt*(eye(2*vesicle.N)-M)*vback;
elseif advType == 2
  % Take a step due to advection
  Xadv = o.translateVinfwNN(Xold,vback);
elseif advType == 3
  Xadv = o.translateVinfwTorch(Xold,vback);
  % filter shape
  Xadv = oc.upsThenFilterShape(Xadv,512,32);
end

dyAdv = mean(Xadv(end/2+1:end))-mean(Xold(end/2+1:end));

% AREA-LENGTH CORRECTION
% disp('Area-Length correction after advection step')
% [XadvC,ifail] = oc.correctAreaAndLength2(Xadv,area0,len0);
% if ifail; disp('Error in AL cannot be corrected!!!'); end;
% Xadv = oc.alignCenterAngle(Xadv,XadvC);

% 2) COMPUTE THE ACTION OF RELAX OP. ON Xold + Xadv
if iExact
disp('ExactSolve')
vesicle = capsules(Xold,[],[],1,1,0); 
G = op.stokesSLmatrix(vesicle);
% Bending, tension and surface divergence
[Ben,Ten,Div] = vesicle.computeDerivs;
M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
rhs = Xadv;
LHS = (eye(2*vesicle.N)-vesicle.kappa*o.dt*(-G*Ben+M*G*Ben));
Xnew = LHS\rhs;
disp('Taking exact relaxation step')
else
Xnew = o.relaxWTorchNet(Xadv);    
end

dyNet = mean(Xnew(end/2+1:end))-mean(Xadv(end/2+1:end));

% First reparameterize
XnewO = Xnew;
for it = 1 : 5
  Xnew = oc.redistributeArcLength(Xnew);
end
Xnew = oc.alignCenterAngle(XnewO,Xnew);

% AREA-LENGTH CORRECTION
disp('Area-Length correction after relaxation step')
[Xnew,ifail] = oc.correctAreaAndLength2(Xnew,area0,len0);
if ifail; disp('Error in AL cannot be corrected!!!'); end;
% Xnew = oc.alignCenterAngle(Xnew,XnewC);


% filter shape
Xnew = oc.upsThenFilterShape(Xnew,512,32);

  
end % DNNsolveTorchNoSplit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,dyNet,dyAdv] = DNNsolve(o,Xold,area0,len0,iExact)
oc = o.oc;
vback = o.vinf(Xold);
advType = 1; % 1: exact, 2: with old net, 3: with Torch net
% 1) COMPUTE THE ACTION OF dt*(1-M) ON Xold  
if advType == 1
  [XoldC,~] = oc.reparametrize(Xold,[],6,20);
  Xold = oc.alignCenterAngle(Xold,XoldC);

  vback = o.vinf(Xold);

  tt = o.tt;
  op = tt.op;
  vesicle = capsules(Xold,[],[],1,1,0);
 
  G = op.stokesSLmatrix(vesicle);
  [~,Ten,Div] = vesicle.computeDerivs;

  M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
  Xadv = Xold + o.dt*(eye(2*vesicle.N)-M)*vback;
elseif advType == 2
  % Take a step due to advection
  Xadv = o.translateVinfwNN(Xold,vback);
elseif advType == 3
  Xadv = o.translateVinfwTorch(Xold,vback);
end

dyAdv = mean(Xadv(end/2+1:end))-mean(Xold(end/2+1:end));

% AREA-LENGTH CORRECTION
% disp('Area-Length correction after advection step')
% [XadvC,ifail] = oc.correctAreaAndLength2(Xadv,area0,len0);
% if ifail; disp('Error in AL cannot be corrected!!!'); end;
% Xadv = oc.alignCenterAngle(Xadv,XadvC);

% 2) COMPUTE THE ACTION OF RELAX OP. ON Xold + Xadv
if iExact
disp('ExactSolve')
vesicle = capsules(Xold,[],[],1,1,0); 
G = op.stokesSLmatrix(vesicle);
% Bending, tension and surface divergence
[Ben,Ten,Div] = vesicle.computeDerivs;
M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
rhs = Xadv;
LHS = (eye(2*vesicle.N)-vesicle.kappa*o.dt*(-G*Ben+M*G*Ben));
Xnew = LHS\rhs;
disp('Taking exact relaxation step')
else
N = numel(Xadv)/2;
Xnew = o.relaxWNNvariableKbDt(Xadv,N,256);
end

dyNet = mean(Xnew(end/2+1:end))-mean(Xadv(end/2+1:end));

% First reparameterize
[XnewC,~] = oc.reparametrize(Xnew,[],6,20);
Xnew = oc.alignCenterAngle(Xnew,XnewC);

% AREA-LENGTH CORRECTION
disp('Area-Length correction after relaxation step')
[Xnew,ifail] = oc.correctAreaAndLength2(Xnew,area0,len0);
if ifail; disp('Error in AL cannot be corrected!!!'); end;
% Xnew = oc.alignCenterAngle(Xnew,XnewC);


  
end % DNNsolve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = DNNsolveTorchSplitTime(o,Xold,area0,len0,iBoth,Xlast,dXdt,timeLastRelax)
oc = o.oc;
vback = o.vinf(Xold);
advType = 3; % 1: exact, 2: with old net, 3: with Torch net
% 1) COMPUTE THE ACTION OF dt*(1-M) ON Xold  
if advType == 1
  [XoldC,~] = oc.reparametrize(Xold,[],6,20);
  Xold = oc.alignCenterAngle(Xold,XoldC);

  tt = o.tt;
  op = tt.op;
  vesicle = capsules(Xlast,[],[],1,1,0);
 
  G = op.stokesSLmatrix(vesicle);
  [~,Ten,Div] = vesicle.computeDerivs;

  M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
  Xadv = Xold + o.dt*(eye(2*vesicle.N)-M)*vback;
elseif advType == 2
  % Take a step due to advection
  Xadv = o.translateVinfwNN(Xold,vback);
elseif advType == 3
  Xadv = o.translateVinfwTorch(Xold,vback);
end

% AREA-LENGTH CORRECTION
disp('Area-Length correction after advection step')
[XadvC,ifail] = oc.correctAreaAndLength2(Xadv,area0,len0);
if ifail; disp('Error in AL cannot be corrected!!!'); end;
Xadv = oc.alignCenterAngle(Xadv,XadvC);

if iBoth
% 2) COMPUTE THE ACTION OF RELAX OP. ON Xold + Xadv
Xnew = o.relaxWTorchNet(Xadv);  
disp('Relaxing with Net')
else
Xnew = Xadv;
% Xnew = Xadv + dXdt * timeLastRelax;
% disp('Relaxing with extrapolation')
end

% AREA-LENGTH CORRECTION
disp('Area-Length correction after relaxation step')
[XnewC,ifail] = oc.correctAreaAndLength2(Xnew,area0,len0);
if ifail; disp('Error in AL cannot be corrected!!!'); end;
Xnew = oc.alignCenterAngle(Xnew,XnewC);


  
end % DNNsolveTorchSplitTime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateVinfwNN(o,Xinput,vinf,Xold,rotate,sortIdx)
% Xinput is equally distributed in arc-length
% Xold as well. So, we add up coordinates of the same points.

N = numel(Xold(:,1))/2;
nv = numel(Xold(1,:));
Nnet = numel(sortIdx(:,1));

% load network files
FCnets = o.MVnets; activeModes = o.velActiveModes; outputSize = o.MVoutSize;

% Approximate the multiplication M*(FFTBasis)     
Z11r = zeros(Nnet,numel(activeModes),nv); Z12r = Z11r;
Z21r = Z11r; Z22r = Z11r;

for k = 1 : numel(activeModes)
  pred = predict(FCnets{k},Xinput)';
  Z11r(:,k,:) = interpft(pred(1:outputSize/4,:),Nnet);
  Z21r(:,k,:) = interpft(pred(outputSize/4+1:outputSize/2,:),Nnet);
  Z12r(:,k,:) = interpft(pred(outputSize/2+1:3*outputSize/4,:),Nnet);
  Z22r(:,k,:) = interpft(pred(3*outputSize/4+1:outputSize,:),Nnet);
end
% upsample vinf
vinfUp = [interpft(vinf(1:end/2,:),Nnet);interpft(vinf(end/2+1:end,:),Nnet)];
% Take fft of the velocity (should be standardized velocity)

% Here the output of MVinf is destandardized
MVinfMat = zeros(2*N,nv);
for k = 1 : nv
  % only sort points and rotate to pi/2 (no translation, no scaling)
  vinfStand = o.standardize(vinfUp(:,k),[0;0],rotate(k),[0;0],1,sortIdx(:,k));
  z = vinfStand(1:end/2)+1i*vinfStand(end/2+1:end);

  zh = fft(z);
  V1 = real(zh(activeModes)); V2 = imag(zh(activeModes));
  % Compute the approximate value of the term M*vinf
  MVinfFull = [Z11r(:,:,k)*V1+Z12r(:,:,k)*V2; Z21r(:,:,k)*V1+Z22r(:,:,k)*V2];
  % Need to destandardize MVinf (take sorting and rotation back)
  MVinf = zeros(size(MVinfFull));
  MVinf([sortIdx(:,k);sortIdx(:,k)+Nnet]) = MVinfFull;
  MVinf = o.rotationOperator(MVinf,-rotate(k),[0; 0]);
  % downsample MVinf
  MVinfMat(:,k) = [interpft(MVinf(1:end/2),N);interpft(MVinf(end/2+1:end),N)];
end
% Update the position
Xnew = Xold + o.dt*vinf-o.dt*MVinfMat;
end % translateVinfwNN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateVinfwTorch(o,Xold,vinf)
% Xinput is equally distributed in arc-length
% Xold as well. So, we add up coordinates of the same points.
N = numel(Xold(:,1))/2;
nv = numel(Xold(1,:));
Nnet = N;
oc = o.oc;
disp('taking advection step with nets')
modes = [(0:Nnet/2-1) (-Nnet/2:-1)];
modesInUse = 48;
modeList = find(abs(modes)<=modesInUse);

% Standardize input
[Xstand,scaling,rotate,rotCent,trans,sortIdx] = o.standardizationStep(Xold,Nnet);

in_param = o.torchAdvInNorm;
out_param = o.torchAdvOutNorm;

% Normalize input
input_list = []; 
cnt = 1;
for imode = modeList
  if imode ~= 1
  input_net = zeros(nv,2,Nnet);  
  x_mean = in_param(imode-1,1);
  x_std = in_param(imode-1,2);
  y_mean = in_param(imode-1,3);
  y_std = in_param(imode-1,4);
  input_net(:,1,:) = (Xstand(1:end/2)-x_mean)/x_std;
  input_net(:,2,:) = (Xstand(end/2+1:end)-y_mean)/y_std;
  input_list{cnt} = py.numpy.array(input_net);
  cnt = cnt + 1;
  end
end % imode

tS = tic;
[Xpredict] = pyrunfile("advect_predict.py","output_list",input_shape=input_list,num_ves=py.int(nv),modesInUse=py.int(modesInUse));
tPyCall = toc(tS);

disp(['Calling python to predict MV takes ' num2str(tPyCall) ' seconds'])
% we have 128 modes
% Approximate the multiplication M*(FFTBasis)     
Z11r = zeros(Nnet,Nnet,nv); Z12r = Z11r;
Z21r = Z11r; Z22r = Z11r;

tS = tic;
for ij = 1 : numel(modeList)-1
  
  imode = modeList(ij+1); % mode index # skipping the first mode
  pred = double(Xpredict{ij}); % size(pred) = [1 2 256]


  % denormalize output
  real_mean = out_param(imode-1,1);
  real_std = out_param(imode-1,2);
  imag_mean = out_param(imode-1,3);
  imag_std = out_param(imode-1,4);
  
  % first channel is real
  pred(1,1,:) = (pred(1,1,:)*real_std) + real_mean;
  % second channel is imaginary
  pred(1,2,:) = (pred(1,2,:)*imag_std) + imag_mean;

  Z11r(:,imode,1) = pred(1,1,1:end/2);
  Z21r(:,imode,1) = pred(1,1,end/2+1:end);
  Z12r(:,imode,1) = pred(1,2,1:end/2);
  Z22r(:,imode,1) = pred(1,2,end/2+1:end);
end
tOrganize = toc(tS);
disp(['Organizing MV output takes ' num2str(tOrganize) ' seconds'])

% Take fft of the velocity (should be standardized velocity)
% only sort points and rotate to pi/2 (no translation, no scaling)
vinfStand = o.standardize(vinf,[0;0],rotate,[0;0],1,sortIdx);
z = vinfStand(1:end/2)+1i*vinfStand(end/2+1:end);

zh = fft(z);
V1 = real(zh); V2 = imag(zh);
% Compute the approximate value of the term M*vinf
MVinfStand = [Z11r*V1+Z12r*V2; Z21r*V1+Z22r*V2];
% Need to destandardize MVinf (take sorting and rotation back)
MVinf = zeros(size(MVinfStand));
MVinf([sortIdx;sortIdx+Nnet]) = MVinfStand;
MVinf = o.rotationOperator(MVinf,-rotate,[0;0]);

Xnew = Xold + o.dt * vinf - o.dt*MVinf;

% XnewStand = Xstand + o.dt*vinfStand - o.dt*MVinf;
% Update the position
% Xnew = o.destandardize(XnewStand,trans,rotate,scaling,sortIdx);
end % translateVinfwTorch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = relaxWNNvariableKbDt(o,Xmid,N,Nnet)
% load network
nets = o.bendNets; muChan1 = o.muChan_bend; sdevChan1 = o.sdevChan_bend;
scale = o.scale_bend; offset = o.offset_bend;
KbDts = o.KbDts; flowKbDt = o.dt * o.kappa;

% number of nets used
nnets = numel(KbDts);

% Get 5 approximations, then interpolate then


% 2) RELAXATION w/ NETWORK
% Standardize vesicle Xmid
[Xstand,scaling,rotate,rotCent, trans,sortIdx] = ...
  o.standardizationStep(Xmid,Nnet);
Xin = (Xstand'-o.colMeans)*o.evects(:,1:o.nCompRelax);

for inet = 1 : nnets
Xinput{inet}(1:16,1,1) = scale(inet,1)*(Xin(1:16)-muChan1(inet,1))/...
  sdevChan1(inet,1)+offset(inet,1);
if o.nCompRelax > 16
Xinput{inet}(17:32,1,1) = scale(inet,2)*(Xin(17:32)-muChan1(inet,2))/...
  sdevChan1(inet,2)+offset(inet,2);
end
end % inet = 1 : nnets


Ypred = zeros(o.nCompRelax,nnets);
for inet = 1 : nnets
% predict PCA coefficients for Xnew
Ypred(1:16,inet) = predict(nets{inet,1},Xinput{inet}(1:16,:,:));
% we use normalized output so take that back
Ypred(1:16,inet) = (Ypred(1:16,inet)-offset(inet,1))*sdevChan1(inet,1)/...
    scale(inet,1)+muChan1(inet,1); 
if o.nCompRelax > 16
  % predict PCA coefficients for Xnew
  Ypred(17:32,inet) = predict(nets{inet,2},Xinput{inet}(17:32,:,:));
  % we use normalized output so take that back
  Ypred(17:32,inet) = (Ypred(17:32,inet)-offset(inet,2))*sdevChan1(inet,2)/...
      scale(inet,2)+muChan1(inet,2); 
end
end

% Build Lagrange Interpolation Function
KbDts = log(KbDts); flowKbDt = log(flowKbDt);
funcl = zeros(nnets,1);
for inet = 1 : nnets
  xm = KbDts(inet); xj = [(1:inet-1)';(inet+1:nnets)'];
  funcl(inet) = prod(flowKbDt-KbDts(xj))/prod(xm-KbDts(xj));
end

YpredInt = Ypred(:,1)*funcl(1);
for inet = 2 : nnets
YpredInt = YpredInt + Ypred(:,inet)*funcl(inet);
end
% reconstruct Xnew using PCA basis
Xpred = (YpredInt'*o.evects(:,1:o.nCompRelax)'+o.colMeans)';

% destandardize
Xpred = o.destandardize(Xpred,trans,rotate,rotCent,scaling,sortIdx);

% downsample to N
Xnew = [interpft(Xpred(1:end/2),N);interpft(Xpred(end/2+1:end),N)];

end % relaxWNNvariableKbDt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = relaxWTorchNet(o,Xmid)  

% 1) RELAXATION w/ NETWORK
% Standardize vesicle Xmid
Xin = zeros(size(Xmid));
nv = numel(Xmid(1,:));
N = numel(Xmid(:,1))/2;

scaling = zeros(nv,1); rotate = zeros(nv,1); 
rotCent = zeros(2,nv); trans = zeros(2,nv);
sortIdx = zeros(N,nv);

for k = 1 : nv
  [Xin(:,k),scaling(k),rotate(k),rotCent(:,k),trans(:,k),sortIdx(:,k)] = ...
    o.standardizationStep(Xmid(:,k),N);
end

% INPUT NORMALIZATION INFO

if N == 128
% For the 625k - June8 - Dt = 1E-5 data
x_mean = -8.430413700466488e-09; 
x_std = 0.06278684735298157;
y_mean = 6.290720477863943e-08; 
y_std = 0.13339413702487946;
elseif N == 32
x_mean = -1.5200416214611323e-07; 
x_std = 0.06278670579195023;
y_mean = -2.5547041104800883e-07; 
y_std = 0.13339416682720184;
end



% INPUT NORMALIZING

% REAL SPACE
Xstand = Xin; % before normalization
Xin(1:end/2,:) = (Xin(1:end/2,:)-x_mean)/x_std;
Xin(end/2+1:end,:) = (Xin(end/2+1:end,:)-y_mean)/y_std;
XinitShape = zeros(nv,2,N);
for k = 1 : nv
XinitShape(k,1,:) = Xin(1:end/2,k)'; 
XinitShape(k,2,:) = Xin(end/2+1:end,k)';
end
XinitConv = py.numpy.array(XinitShape);


% OUTPUT

% [DXpredictStand] = pyrunfile("relax_predict_DIFF_IT3_dt1E5.py", "predicted_shape", input_shape=XinitConv);

% June8 - Dt1E5
if N == 128
[DXpredictStand] = pyrunfile("relax_predict_DIFF_June8_dt1E5.py", "predicted_shape", input_shape=XinitConv);
elseif N == 32
[DXpredictStand] = pyrunfile("32modes_relax_predict_DIFF_dt1E5.py", "predicted_shape", input_shape=XinitConv);
end

% June8 - Dt5E5
% [DXpredictStand] = pyrunfile("relax_predict_DIFF_June8_dt5E5.py", "predicted_shape", input_shape=XinitConv);

% June8 - Dt1E5
% [DXpredictStand] = pyrunfile("relax_predict_DIFF_June8_dt1E4.py", "predicted_shape", input_shape=XinitConv);


% 625k IT3 DIFF net
% x_mean = -1.1276009015404043e-09; 
% x_std = 0.00020668248180299997;
% y_mean = 1.3305034157751194e-11; 
% y_std = 0.00017724868666846305;

if N == 128
% For the 625k - June8 - Dt = 1E-5 data
x_mean = -2.884585348361668e-10; 
x_std = 0.00020574081281665713;
y_mean = -5.137390512999218e-10; 
y_std = 0.0001763451291481033;
elseif N == 32
x_mean = -2.329148207635967e-09; 
x_std = 0.00020403489179443568;
y_mean = -1.5361016902915026e-09; 
y_std = 0.00017457373905926943;
end

% For the 625k - June8 - Dt = 5E-5 data
% x_mean = 2.9649513400009653e-10; 
% x_std = 0.0008968145120888948;
% y_mean = 1.698907792224702e-09; 
% y_std = 0.0007752174278721213;

% For the 625k - June8 - Dt = 1E-4 data
% x_mean = 4.258172037197028e-09; 
% x_std = 0.001633652369491756;
% y_mean = 7.698989001880818e-09; 
% y_std = 0.0014213572721928358;



DXpred = zeros(size(Xin));
DXpredictStand = double(DXpredictStand);
Xnew = zeros(size(Xmid));

for k = 1 : nv
% normalize output
DXpred(1:end/2,k) = DXpredictStand(k,1,:)*x_std + x_mean;
DXpred(end/2+1:end,k) = DXpredictStand(k,2,:)*y_std + y_mean;


DXpred(:,k) = DXpred(:,k)/1E-5 * o.dt;
Xpred = Xstand(:,k) + DXpred(:,k);

Xnew(:,k) = o.destandardize(Xpred,trans(:,k),rotate(k),rotCent(:,k),...
    scaling(k),sortIdx(:,k));
end



end % relaxWTorchNet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew, timeStand, timePred, timeDestand] = relaxWTorchBenchmark(o,Xmid)  

% 1) RELAXATION w/ NETWORK
tStandI = tic;
[Xin,scaling,rotate,rotCent,trans,sortIdx] = ...
  o.standardizationStep(Xmid,128);

% INPUT NORMALIZATION INFO

% For the 625k (IT3) data
% x_mean = -3.775884049872502e-09; 
% x_std = 0.06278640776872635;
% y_mean = -5.037749133407488e-07; 
% y_std = 0.1333947628736496;

% For the 625k - June8 - Dt = 1E-5 data
x_mean = -8.430413700466488e-09; 
x_std = 0.06278684735298157;
y_mean = 6.290720477863943e-08; 
y_std = 0.13339413702487946;

% For the 625k - June8 - Dt = 5E-5 data
% x_mean = 2.4658157826706883e-07; 
% x_std = 0.06278616935014725;
% y_mean = -4.5408405924263207e-08; 
% y_std = 0.13339488208293915;

% For the 625k - June8 - Dt = 1E-4 data
% x_mean = 7.684710645605719e-09; 
% x_std = 0.06278636306524277;
% y_mean = 7.071167829053593e-08; 
% y_std = 0.13339479267597198;



% INPUT NORMALIZING

% REAL SPACE
Xstand = Xin; % before normalization
Xin(1:end/2) = (Xin(1:end/2)-x_mean)/x_std;
Xin(end/2+1:end) = (Xin(end/2+1:end)-y_mean)/y_std;
XinitShape = zeros(1,2,128);
XinitShape(1,1,:) = Xin(1:end/2)'; 
XinitShape(1,2,:) = Xin(end/2+1:end)';

tStandO = toc(tStandI);
timeStand = tStandO;

XinitConv = py.numpy.array(XinitShape);

tPredI = tic;
% OUTPUT

% [DXpredictStand] = pyrunfile("relax_predict_DIFF_IT3_dt1E5.py", "predicted_shape", input_shape=XinitConv);

% June8 - Dt1E5
[DXpredictStand] = pyrunfile("relax_predict_DIFF_June8_dt1E5.py", "predicted_shape", input_shape=XinitConv);

% June8 - Dt5E5
% [DXpredictStand] = pyrunfile("relax_predict_DIFF_June8_dt5E5.py", "predicted_shape", input_shape=XinitConv);

% June8 - Dt1E5
% [DXpredictStand] = pyrunfile("relax_predict_DIFF_June8_dt1E4.py", "predicted_shape", input_shape=XinitConv);

tPredO = toc(tPredI);
timePred = tPredO;

% 625k IT3 DIFF net
% x_mean = -1.1276009015404043e-09; 
% x_std = 0.00020668248180299997;
% y_mean = 1.3305034157751194e-11; 
% y_std = 0.00017724868666846305;

% For the 625k - June8 - Dt = 1E-5 data
x_mean = -2.884585348361668e-10; 
x_std = 0.00020574081281665713;
y_mean = -5.137390512999218e-10; 
y_std = 0.0001763451291481033;

% For the 625k - June8 - Dt = 5E-5 data
% x_mean = 2.9649513400009653e-10; 
% x_std = 0.0008968145120888948;
% y_mean = 1.698907792224702e-09; 
% y_std = 0.0007752174278721213;

% For the 625k - June8 - Dt = 1E-4 data
% x_mean = 4.258172037197028e-09; 
% x_std = 0.001633652369491756;
% y_mean = 7.698989001880818e-09; 
% y_std = 0.0014213572721928358;


DXpred = zeros(size(Xin));
DXpredictStand = double(DXpredictStand);

tDestandI = tic;
% normalize output
DXpred(1:end/2) = DXpredictStand(1,1,:)*x_std + x_mean;
DXpred(end/2+1:end) = DXpredictStand(1,2,:)*y_std + y_mean;

DXpredScaled = DXpred; %/1E-5 * o.dt;
Xpred = Xstand + DXpredScaled;
Xnew = o.destandardize(Xpred,trans,rotate,rotCent,scaling,sortIdx);
tDestandO = toc(tDestandI);
timeDestand = tDestandO ;

end % relaxBenchmark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xinput = prepareInputForNet(o,X,netType)
% Normalize input
if strcmp(netType,'tensionOnFourier')
scale = o.scale_tenV;
muChan = o.muChan_tenV;
sdevChan = o.sdevChan_tenV;
offset = o.offset_tenV;
elseif strcmp(netType,'tensionSelfBend')
scale = o.scale_tenBend;
muChan = o.muChan_tenBend;
sdevChan = o.sdevChan_tenBend;
offset = o.offset_tenBend;    
elseif strcmp(netType,'advection')
scale = o.scale_MV;
muChan = o.muChan_MV;
sdevChan = o.sdevChan_MV;
offset = o.offset_MV;    
elseif strcmp(netType,'relaxation')
scale = o.scale_bend;
muChan = o.muChan_bend;
sdevChan = o.sdevChan_bend;
offset = o.offset_bend;     
end

if strcmp(netType,'relaxation')
  % find PCA coefficients  
  Xinput(:,1,1,1) = (X'-o.colMeans)*o.evects(:,1:o.nCompRelax);  
  Xinput(1:16,1,1,1) = scale(1)*(Xinput(1:16,1,1,1)-muChan(1))/...
    sdevChan(1)+offset(1);
  if o.nCompRelax > 16
    Xinput(17:32,1,1,1) = scale(2)*(Xinput(17:32,1,1,1)-muChan(2))/...
        sdevChan(2)+offset(2);
  end
  if o.nCompRelax > 32
    Xinput(33:48,1,1,1) = scale(3)*(Xinput(33:48,1,1,1)-muChan(3))/...
        sdevChan(3)+offset(3);
  end
  if o.nCompRelax > 48
    Xinput(49:64,1,1,1) = scale(4)*(Xinput(49:64,1,1,1)-muChan(4))/...
        sdevChan(4)+offset(4);
  end
else
  % find PCA coefficients  
  Xinput(:,1,1,1) = (X'-o.colMeans)*o.evects(:,1:o.nComp);  
  Xinput(:,1,1,1) = scale*(Xinput(:,1,1,1)-muChan)/sdevChan+offset;        
end     
end % prepareInputForNet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(o,Xin,Nnet)
oc = o.oc;
N = numel(Xin)/2;
if Nnet ~= N
  Xin = [interpft(Xin(1:end/2),Nnet);interpft(Xin(end/2+1:end),Nnet)];    
end

% Equally distribute points in arc-length
for iter = 1 : 10
  [Xin,~,~] = oc.redistributeArcLength(Xin);
end


X = Xin;
[trans,rotate,rotCent,scaling,sortIdx] = o.referenceValues(X);

% Fix misalignment in center and angle due to reparametrization
% X = oc.alignCenterAngle(Xin,X);

% standardize angle, center, scaling and point order

X = o.standardize(X,trans,rotate,rotCent,scaling,sortIdx);
end % standardizationStep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XrotSort = standardize(o,X,translation,rotation,rotCent,scaling,sortIdx)
N = numel(sortIdx);

% translate, rotate and scale configuration

Xrotated = o.rotationOperator(X,rotation,rotCent);   
Xrotated = o.translateOp(Xrotated,translation);

% now order the points
XrotSort = [Xrotated(sortIdx);Xrotated(sortIdx+N)];

XrotSort = scaling*XrotSort;

end % standardize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = destandardize(o,XrotSort,translation,rotation,rotCent,scaling,sortIdx)

N = numel(sortIdx);    
    
% scaling back
XrotSort = XrotSort/scaling;

% change ordering back 
X = zeros(size(XrotSort));
X([sortIdx;sortIdx+N]) = XrotSort;

% take translation back
X = o.translateOp(X,-translation);

% take rotation back
X = o.rotationOperator(X,-rotation,rotCent);


end % destandardize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [translation,rotation,rotCent,scaling,sortIdx] = referenceValues(o,Xref)
oc = o.oc;
N = numel(Xref)/2;

% find translation, rotation and scaling
center = oc.getPhysicalCenterShan(Xref);
V = oc.getPrincAxesGivenCentroid(Xref,center);
% % find rotation angle
w = [0;1]; % y-axis
rotation = atan2(w(2)*V(1)-w(1)*V(2), w(1)*V(1)+w(2)*V(2));


% translation = [-mean(Xref(1:end/2));-mean(Xref(end/2+1:end))];
% rotation = pi/2-oc.getIncAngle2(Xref);
       
% find the ordering of the points
rotCent = center;
Xref = o.rotationOperator(Xref, rotation, center);
center = oc.getPhysicalCenterShan(Xref);
translation = -center;

Xref = o.translateOp(Xref, translation);

firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

% amount of scaling
[~,~,length] = oc.geomProp(Xref);
scaling = 1/length;
end % referenceValues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xrot = rotationOperator(o,X,theta, rotCent)
% Get x-y coordinates
Xrot = zeros(size(X));
x = X(1:end/2)-rotCent(1); y = X(end/2+1:end)-rotCent(2);

% Rotated shape
xrot = (x)*cos(theta) - (y)*sin(theta);
yrot = (x)*sin(theta) + (y)*cos(theta);

Xrot(1:end/2) = xrot+rotCent(1);
Xrot(end/2+1:end) = yrot+rotCent(2);
end % rotationOperator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateOp(o,X,transXY)
Xnew = zeros(size(X));
Xnew(1:end/2) = X(1:end/2)+transXY(1);
Xnew(end/2+1:end) = X(end/2+1:end)+transXY(2);
end  % translateOp  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaNew,RSnew] = DNNlikeExactSolve(o,Xold,tenOld,etaOld,RSold)
 
tt = o.tt;
op = tt.op;
if o.confined
  opWall = tt.opWall;
  invWallsMat = tt.bdiagWall;
  walls = o.walls;
  uWalls = walls.u;
  nvbd = walls.nv;
  Nbd = walls.N;
  etaNew = zeros(2*Nbd,nvbd);
  RSnew = zeros(3,nvbd);
else
  vback = o.vinf(Xold);
  etaNew = []; RSnew = [];
end

% Semi-implicit time stepping w/splitting
vesicle = capsules(Xold,[],[],o.kappa,1,0);
nv = vesicle.nv;
N = vesicle.N;
% Build the operators
G = op.stokesSLmatrix(vesicle);
[~,Ten,Div] = vesicle.computeDerivs;
if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xold(1:end/2,:),Nup);interpft(Xold(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  SLPnoCorr = [];
end

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

% 1) EXPLICIT TENSION AT THE CURRENT STEP
% neglect near interactions

% Get the near structure, so that we know what to neglect
if o.confined
  [NearV2V,NearV2W] = vesicle.getZone(walls,3);
  [~,NearW2V] = walls.getZone(vesicle,2);
  if nvbd > 1
    NearW2W = o.NearW2W;  
  end % if nvbd > 1
else
  NearV2V = vesicle.getZone([],1);    
end % if o.confined

if nv > 1
% Far-field due to traction jump
kernelDirect = @op.exactStokesSL;
if tt.fmm
  kernel = @op.exactStokesSL;
else
  kernel = @op.exactStokesSL;
end

if o.useTrueNear
  farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLPnoCorr,NearV2V,...
      kernel,kernelDirect,vesicle,true,false);
else
  [~,nearField,farFieldtracJump] = op.divideNearFarSLP(vesicle,tracJump,SLPnoCorr,...
      NearV2V,kernel,kernelDirect,vesicle,true);      
  if o.useNear
    farFieldtracJump = farFieldtracJump + nearField;
  end
end

else
farFieldtracJump = zeros(2*N,1);
end

% Far-field due to walls or background flow
if o.confined
  kernelDirect = @opWall.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWall.exactStokesDL;
  else
    kernel = @opWall.exactStokesDL;
  end
  if o.useTrueNear
    vback = opWall.nearSingInt(walls,etaOld,[],NearW2V,kernel,...
        kernelDirect,vesicle,false,false);
  else
    [~,nearField,vback] = opWall.divideNearFarSLP(walls,etaOld,[],NearW2V,...
        kernel,kernelDirect,vesicle,false);
    if o.useNear
      vback = vback + nearField;
    end
  end
  
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
  end
end

tenNew = zeros(N,nv);
for k = 1 : nv
  LHS = (Div(:,:,k)*G(:,:,k)*Ten(:,:,k));
  % DNN approximates inv(Div*G*Ten)*G*(-Ben)*x and inv(Div*G*Ten)*Div*Fourier
  
  % SLP due to bending force on itself 
  selfBend = G(:,:,k)*fBend(:,k);
  
  RHS = -Div(:,:,k)*(vback(:,k)+farFieldtracJump(:,k)+selfBend);
  %[inv(LHS)*(Div(:,:,k)*(vback(:,k)+farFieldtracJump(:,k))) inv(LHS)*(Div(:,:,k)*selfBend)]
  %pause
  tenNew(:,k) = LHS\RHS;
end
% compute force due to tension
fTen = vesicle.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;

if o.confined
  % 2) SOLVE FOR DENSITY and RS ON WALLS
  
  % vesicle2wall interactions
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  [~,nearField,ves2wallInt] = op.divideNearFarSLP(vesicle,tracJump,[],NearV2W,...
      kernel,kernelDirect,walls,0);
  if o.useNear
    ves2wallInt = ves2wallInt + nearField;
  end
  RHS = [uWalls(:)-ves2wallInt(:);zeros(3*(nvbd-1),1)];
  etaRS = invWallsMat*RHS;
  for iw = 1 : nvbd
    etaNew(:,iw) = etaRS((iw-1)*2*Nbd+1:iw*2*Nbd);  
    if iw <= nvbd-1
      RSnew(:,iw+1) = etaRS(2*Nbd*nvbd+(iw-1)*3+1:2*Nbd*nvbd+iw*3);
    end
  end
  % Update vback due to walls
  kernelDirect = @opWall.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWall.exactStokesDL;
  else
    kernel = @opWall.exactStokesDL;
  end
  [~,nearField,vback] = opWall.divideNearFarSLP(walls,etaNew,[],NearW2V,...
      kernel,kernelDirect,vesicle,false);
  if o.useNear
    vback = vback + nearField;
  end
  for k = 2:nvbd
    stokeslet = RSnew(1:2,k); rotlet = RSnew(3,k);  
    vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
  end
end % if iconfined

% 3) SOLVE FOR NEXT TIME STEP USING OPERATOR SPLITTING
% EXPLICIT-TREATMENT OF INTER-VESICLE BENDING AND TENSION
% TENSION AND X are coupled (semi-implicit).
Xnew = zeros(size(Xold));
if nv > 1
% compute traction jump due to explicit tension and old bending
kernelDirect = @op.exactStokesSL;
if tt.fmm
  kernel = @op.exactStokesSL;
else
  kernel = @op.exactStokesSL;
end
[~,nearField,farFieldtracJump] = op.divideNearFarSLP(...
        vesicle,tracJump,[],NearV2V,kernel,kernelDirect,vesicle,1);
if o.useNear 
  farFieldtracJump = farFieldtracJump + nearField;
end
else
farFieldtracJump = zeros(2*N,1);
end

% Walls2Ves is already computed if confined (vback)    
for k = 1 : nv
  % First, translate with vext
  vext = vback(:,k) + farFieldtracJump(:,k);
  %Xmid = o.translateVinf(vext,o.dt,Xold(:,k),op);
  
  % Second, solve bending problem
  vesicle = capsules(Xold(:,k),[],[],1,1,0); 
  Xnew(:,k) = o.relaxExactSolve(vesicle,vext,o.dt,Xold(:,k),op);
end


end % DNNlikeExactSolve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = relaxExactSolve(o,vesicle,vinf,dt,Xold,op)
% HERE, BENDING IN TENSION SOLVE IS IMPLICIT

% SLP
G = op.stokesSLmatrix(vesicle);
% Bending, tension and surface divergence
[Ben,Ten,Div] = vesicle.computeDerivs;
M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
rhs = Xold + dt*(eye(2*vesicle.N)-M)*vinf;
LHS = (eye(2*vesicle.N)-vesicle.kappa*dt*(-G*Ben+M*G*Ben));
Xnew = LHS\rhs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateVinf(o,vinf,dt,Xold,op)
% ADVECTION PART FOR OPERATOR SPLITTING
vesicle = capsules(Xold,[],[],1,1,0);
 
G = op.stokesSLmatrix(vesicle);
[~,Ten,Div] = vesicle.computeDerivs;

M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
Xnew = Xold + dt*(eye(2*vesicle.N)-M)*vinf;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vinf = setBgFlow(o,bgFlow,speed,chanWidth)

vinf = @(X) zeros(size(X));      
if strcmp(bgFlow,'relax')
  vinf = @(X) zeros(size(X));  % relaxation
elseif strcmp(bgFlow,'shear') 
  vinf = @(X) speed*[X(end/2+1:end,:);zeros(size(X(1:end/2,:)))]; 
elseif strcmp(bgFlow,'tayGreen')
  vinf = @(X) speed*[sin(X(1:end/2,:)).*cos(X(end/2+1:end,:));-...
    cos(X(1:end/2,:)).*sin(X(end/2+1:end,:))]; % Taylor-Green
elseif strcmp(bgFlow,'parabolic')
  % in our earlier simulations curvature was 1.3, now 0.65  
  vinf = @(X) [speed*(1-(X(end/2+1:end,:)/chanWidth).^2);...
      zeros(size(X(1:end/2,:)))];
elseif strcmp(bgFlow,'rotation')
  vinf = @(X) [-sin(atan2(X(end/2+1:end,:),X(1:end/2,:)))./sqrt(X(1:end/2,:).^2+X(end/2+1:end,:).^2);...
    cos(atan2(X(end/2+1:end,:),X(1:end/2,:)))./sqrt(X(1:end/2,:).^2+X(end/2+1:end,:).^2)]*speed;
end
    
end % setBgFlow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tt = buildTstep(o,X,pramsIn)
N = pramsIn.N; nv = pramsIn.nv;

prams.N = N;
prams.nv = nv;
prams.Nbd = pramsIn.Nbd;
prams.nvbd = pramsIn.nvbd;
prams.kappa = pramsIn.kappa;
options.diffDiscWalls = 0;
prams.NbdInt = 0; prams.NbdExt = 0; prams.nvbdInt = 0; prams.nvbdExt = 0;
options.verbose = 0;
options.saveData = 0;
options.usePlot = 0;
options.track = 0;
options.quiver = 0;
options.axis = 0;
prams.T = pramsIn.Th;
prams.gmresTol = 1e-8;
prams.m = pramsIn.Th/pramsIn.dt;
prams.errTol = 1e-1;
options.tracers = 0;
options.timeAdap = 0;
options.order = 1;
prams.areaLenTol = 1e-1;
options.fmm = pramsIn.fmm;
options.fmmDLP = pramsIn.fmmDLP;

options.correctShape = 1;
options.adhesion = 0;
options.repulsion = 0;
if strcmp(pramsIn.bgFlow,'couette')
options.confined = true;
options.farField = 'couette';
options.farFieldSpeed = pramsIn.speed;
else
options.confined = false;
end

[options,prams] = initVes2D(options,prams);
    
om = monitor(X,options,prams);
tt = tstep(options,prams,om);
if ~options.confined % if not confined, then assign vinf
tt.farField = o.setBgFlow(pramsIn.bgFlow,pramsIn.speed); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % methods

end % dnnTools
