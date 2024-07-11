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

load ./shannets/nearInterp_fft_in_param.mat
load ./shannets/nearInterp_fft_out_param.mat

load finalShear2.mat
X = Xhist;
sig = sigStore;

Nnet = 128;
N = 128;
nv = 2;
oc = curve;

vesicle = capsules(X,sig,[],1,1,0);
tracJump = vesicle.tracJump(X,sig);

%% wrong mode cutting in tracJump
wrongTrac = zeros(2*N,nv);
for k = 1 : nv
z = tracJump(1:end/2,k) + 1i*tracJump(end/2+1:end,k);
z = fft(z); z = z(1:16);
z2 = real(z) + 1i*imag(z);
z2 = ifft([z2;zeros(112,1)]);
wrongTrac(:,k) = [real(z2);imag(z2)];
end

%% 
maxLayerDist = @(h) sqrt(h);
maxLayerDistNet = sqrt(1/Nnet); % length = 1, h = 1/Nnet;
nlayers = 3;

nmodes = 128;
theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N * exp(1i*theta*ks');

op = poten(N);

%% Calculate near field with exact calculation
[~,tang] = oc.diffProp(X);
% get x and y components of normal vector at each point
nx = tang(N+1:2*N,:);
ny = -tang(1:N,:);

% Points where velocity is calculated involve the points on vesicle
tracersX = zeros(2*N, nlayers-1,nv);
tracersX(:,1,1) = X(:,1);
tracersX(:,1,2) = X(:,2);

h = vesicle.length/vesicle.N;  % arc-length spacing
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist(h);

for k = 1 : nv
  for il = 2 : nlayers
    tracersX(:,il,k) = [X(1:end/2,k)+nx(:,k)*dlayer(il);X(end/2+1:end,k)+ny(:,k)*dlayer(il)];
  end
end

% do it for both vesicles
G = op.stokesSLmatrix(vesicle);
velx_true = zeros(N,3,nv);
vely_true = zeros(N,3,nv);

kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;


for k = 1 : nv
tracers.N = N;
tracers.nv = 2;
tracers.X = tracersX(:,2:nlayers,k);

ves2 = capsules(X(:,k),sig(:,k),[],1,1,0);

[~,NearV2T] = ves2.getZone(tracers,2);
SLP = @(X) op.exactStokesSLdiag(ves2,G(:,:,k),X);

V = op.nearSingInt(ves2,wrongTrac(:,k),SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);
velx_true(:,2:end,k) = V(1:end/2,:);
vely_true(:,2:end,k) = V(end/2+1:end,:);

Vs = G(:,:,k)*wrongTrac(:,k);
velx_true(:,1,k) = Vs(1:end/2);
vely_true(:,1,k) = Vs(end/2+1:end);
end

%% Calculate near field with network

% standardize input
Xstand = zeros(size(X));
scaling = zeros(nv,1);
rotate = zeros(nv,1);
rotCent = zeros(2,nv);
trans = zeros(2,nv);
sortIdx = zeros(N,nv);

tracersXnet = zeros(2*N,3,nv);

dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDistNet;
for k = 1 : nv
  [Xstand(:,k),scaling(k),rotate(k),rotCent(:,k),trans(:,k),sortIdx(:,k)] = standardizationStep(X(:,k),128);
  [~,tang] = oc.diffProp(Xstand(:,k));
  nx = tang(Nnet+1:2*Nnet);
  ny = -tang(1:Nnet);

  tracersXnet(:,1,k) = Xstand(:,k);
  for il = 2 : nlayers 
    tracersXnet(:,il,k) = [Xstand(1:end/2,k)+nx*dlayer(il); Xstand(end/2+1:end,k)+ny*dlayer(il)];
  end
end

% Normalize input
input_net = zeros(nv,2,128);

for k = 1 : nv
  input_net(k,1,:) = (Xstand(1:end/2,k)-in_param(1,1))/in_param(1,2);
  input_net(k,2,:) = (Xstand(end/2+1:end,k)-in_param(1,3))/in_param(1,4);
end

fstandRe = zeros(16, nv);
fstandIm = zeros(16, nv);
tracStand = zeros(256,nv);
tracStand2 = zeros(256,nv);
for k = 1 : nv
  % fstand = standardize(tracJump(:,k),[0;0], rotate(k), [0;0], 1, sortIdx(:,k));
  fRot = rotationOperator(tracJump(:,k),rotate(k),[0;0]);
  fstand = [fRot(sortIdx(:,k));fRot(sortIdx(:,k)+N)];
  tracStand(:,k) = fstand;
  z = fstand(1:end/2) + 1i*fstand(end/2+1:end);
  zh = fft(z);
  fstandRe(:,k) = real(zh(1:16)); 
  fstandIm(:,k) = imag(zh(1:16));
  z2 = fstandRe(:,k) + 1i*fstandIm(:,k);
  z2 = [ifft([z2;zeros(112,1)])];
  tracStand2(:,k) = [real(z2);imag(z2)];
end

input_conv = py.numpy.array(input_net);
[Xpredict] = pyrunfile("near_vel_predict.py","output_list",input_shape=input_conv,num_ves=py.int(nv));

for k = 1 : nv
velx_real{k} = zeros(128,16,3);
vely_real{k} = zeros(128,16,3);
velx_imag{k} = zeros(128,16,3);
vely_imag{k} = zeros(128,16,3);
end

% denormalize output
for imode = 1 : 16
  pred = double(Xpredict{imode});
  % its size is (nv x 12 x 128) 
  % channel 1-3: vx_real_layers 0, 1, 2
  % channel 4-6; vy_real_layers 0, 1, 2
  % channel 7-9: vx_imag_layers 0, 1, 2
  % channel 10-12: vy_imag_layers 0, 1, 2

  % denormalize output
  for k = 1 : nv
    velx_real{k}(:,imode,1) = (pred(k,1,:)*out_param(imode,2,1))  + out_param(imode,1,1);
    velx_real{k}(:,imode,2) = (pred(k,2,:)*out_param(imode,2,2))  + out_param(imode,1,2);
    velx_real{k}(:,imode,3) = (pred(k,3,:)*out_param(imode,2,3))  + out_param(imode,1,3);
    vely_real{k}(:,imode,1) = (pred(k,4,:)*out_param(imode,2,4))  + out_param(imode,1,4);
    vely_real{k}(:,imode,2) = (pred(k,5,:)*out_param(imode,2,5))  + out_param(imode,1,5);
    vely_real{k}(:,imode,3) = (pred(k,6,:)*out_param(imode,2,6))  + out_param(imode,1,6);

    velx_imag{k}(:,imode,1) = (pred(k,7,:)*out_param(imode,2,7))  + out_param(imode,1,7);
    velx_imag{k}(:,imode,2) = (pred(k,8,:)*out_param(imode,2,8))  + out_param(imode,1,8);
    velx_imag{k}(:,imode,3) = (pred(k,9,:)*out_param(imode,2,9))  + out_param(imode,1,9);
    vely_imag{k}(:,imode,1) = (pred(k,10,:)*out_param(imode,2,10))  + out_param(imode,1,10);
    vely_imag{k}(:,imode,2) = (pred(k,11,:)*out_param(imode,2,11))  + out_param(imode,1,11);
    vely_imag{k}(:,imode,3) = (pred(k,12,:)*out_param(imode,2,12))  + out_param(imode,1,12);
  end
end


velx = zeros(128,3,nv); 
vely = zeros(128,3,nv);
xlayers = zeros(128,3,nv);
ylayers = zeros(128,3,nv);
velx_stand = zeros(128,3,nv);
vely_stand = zeros(128,3,nv);
for k = 1 : nv
  
  for il = 1 : 3
     velx_stand(:,il,k) = velx_real{k}(:,:,il) * fstandRe(:,k) + velx_imag{k}(:,:,il)*fstandIm(:,k);
     vely_stand(:,il,k) = vely_real{k}(:,:,il) * fstandRe(:,k) + vely_imag{k}(:,:,il)*fstandIm(:,k);
     vx = zeros(128,1);
     vy = zeros(128,1);
     vx(sortIdx(:,k)) = velx_stand(:,il,k);
     vy(sortIdx(:,k)) = vely_stand(:,il,k);
     
     VelBefRot = [vx;vy];
     
     VelRot = rotationOperator(VelBefRot, -rotate(k), [0;0]);
     velx(:,il,k) = VelRot(1:end/2); vely(:,il,k) = VelRot(end/2+1:end);
        
     Xl = destandardize(tracersXnet(:,il,k),trans(:,k),rotate(k),rotCent(:,k),scaling(k),sortIdx(:,k));
     xlayers(:,il,k) = Xl(1:end/2);
     ylayers(:,il,k) = Xl(end/2+1:end);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(Xin,Nnet)
oc = curve;
N = numel(Xin)/2;
if Nnet ~= N
  Xin = [interpft(Xin(1:end/2),Nnet);interpft(Xin(end/2+1:end),Nnet)];    
end

% Equally distribute points in arc-length
for iter = 1 : 10
  [Xin,~,~] = oc.redistributeArcLength(Xin);
end


X = Xin;
[trans,rotate,rotCent,scaling,sortIdx] = referenceValues(X);

% Fix misalignment in center and angle due to reparametrization
% X = oc.alignCenterAngle(Xin,X);

% standardize angle, center, scaling and point order

X = standardize(X,trans,rotate,rotCent,scaling,sortIdx);
end % standardizationStep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XrotSort = standardize(X,translation,rotation,rotCent,scaling,sortIdx)
N = numel(sortIdx);

% translate, rotate and scale configuration

Xrotated = rotationOperator(X,rotation,rotCent);   
Xrotated = translateOp(Xrotated,translation);

% now order the points
XrotSort = [Xrotated(sortIdx);Xrotated(sortIdx+N)];

XrotSort = scaling*XrotSort;

end % standardize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = destandardize(XrotSort,translation,rotation,rotCent,scaling,sortIdx)

N = numel(sortIdx);    
    
% scaling back
XrotSort = XrotSort/scaling;

% change ordering back 
X = zeros(size(XrotSort));
X([sortIdx;sortIdx+N]) = XrotSort;

% take translation back
X = translateOp(X,-translation);

% take rotation back
X = rotationOperator(X,-rotation,rotCent);


end % destandardize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [translation,rotation,rotCent,scaling,sortIdx] = referenceValues(Xref)
oc = curve;
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
Xref = rotationOperator(Xref, rotation, center);
center = oc.getPhysicalCenterShan(Xref);
translation = -center;

Xref = translateOp(Xref, translation);

firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

% amount of scaling
[~,~,length] = oc.geomProp(Xref);
scaling = 1/length;
end % referenceValues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xrot = rotationOperator(X,theta, rotCent)
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
function Xnew = translateOp(X,transXY)
Xnew = zeros(size(X));
Xnew(1:end/2) = X(1:end/2)+transXY(1);
Xnew(end/2+1:end) = X(end/2+1:end)+transXY(2);
end  % translateOp  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%