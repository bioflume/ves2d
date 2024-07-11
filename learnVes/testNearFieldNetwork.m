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

load finalShearX.mat
X = Xhist;

Nnet = 128;
oc = curve;


tracJump = zeros(size(X));
tracJump(:,1) = [X(end/2+1:end,1);zeros(Nnet,1)];
tracJump(:,2) = [zeros(Nnet,1); X(1:end/2,2)];


load ./shannets/nearInterp_fft_in_param.mat
load ./shannets/nearInterp_fft_out_param.mat

maxLayerDist = sqrt(1/Nnet); % length = 1, h = 1/Nnet;
% Predictions on three layers
nlayers = 3;
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist;

% standardize input
Xstand = zeros(size(X));
nv = numel(X(1,:));
scaling = zeros(nv,1);
rotate = zeros(nv,1);
rotCent = zeros(2,nv);
trans = zeros(2,nv);
sortIdx = zeros(Nnet,nv);

tracersX = zeros(2*Nnet,3,nv);
for k = 1 : nv
  [Xstand(:,k),scaling(k),rotate(k),rotCent(:,k),trans(:,k),sortIdx(:,k)] = standardizationStep(X(:,k),128);
  [~,tang] = oc.diffProp(Xstand(:,k));
  nx = tang(Nnet+1:2*Nnet);
  ny = -tang(1:Nnet);

  tracersX(:,1,k) = Xstand(:,k);
  for il = 2 : nlayers 
    tracersX(:,il,k) = [Xstand(1:end/2,k)+nx*dlayer(il); Xstand(end/2+1:end,k)+ny*dlayer(il)];
  end
end

% figure(1); clf;
% subplot(2,1,1);
% plot(X(1:end/2,1),X(end/2+1:end,1),'linewidth',2)
% axis equal
% 
% subplot(2,1,2)
% plot(Xstand(1:end/2,1),Xstand(end/2+1:end,1),'linewidth',2)
% axis equal
% hold on
% plot(tracersX(1:end/2,:,1),tracersX(end/2+1:end,:,1),'k.','markersize',5)
% 
% figure(2); clf;
% subplot(2,1,1);
% plot(X(1:end/2,2),X(end/2+1:end,2),'linewidth',2)
% axis equal
% 
% subplot(2,1,2)
% plot(Xstand(1:end/2,2),Xstand(end/2+1:end,2),'linewidth',2)
% axis equal
% hold on
% plot(tracersX(1:end/2,:,2),tracersX(end/2+1:end,:,2),'k.','markersize',5)
% 
% pause


% Normalize input
input_net = zeros(nv,2,128);

for k = 1 : nv
  input_net(k,1,:) = (Xstand(1:end/2,k)-in_param(1,1))/in_param(1,2);
  input_net(k,2,:) = (Xstand(end/2+1:end,k)-in_param(1,3))/in_param(1,4);
end

%% standardize tracJump
fstandRe = zeros(16, nv);
fstandIm = zeros(16, nv);
for k = 1 : nv
  fstand = standardize(tracJump(:,k),[0;0], rotate(k), [0;0], 1, sortIdx(:,k));
  z = fstand(1:end/2) + 1i*fstand(end/2+1:end);
  zh = fft(z);
  fstandRe(:,k) = real(zh(1:16)); 
  fstandIm(:,k) = imag(zh(1:16));
end


input_conv = py.numpy.array(input_net);
[Xpredict] = pyrunfile("near_vel_predict.py","output_list",input_shape=input_conv,num_ves=py.int(nv));
%%

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
for k = 1 : nv
  velx_stand = zeros(128,3);
  vely_stand = zeros(128,3);
  for il = 1 : 3
     velx_stand(:,il) = velx_real{k}(:,:,il) * fstandRe(:,k) + velx_imag{k}(:,:,il)*fstandIm(:,k);
     vely_stand(:,il) = vely_real{k}(:,:,il) * fstandRe(:,k) + vely_imag{k}(:,:,il)*fstandIm(:,k);
     vx = zeros(128,1);
     vy = zeros(128,1);
     vx(sortIdx(:,k)) = velx_stand(:,il);
     vy(sortIdx(:,k)) = vely_stand(:,il);
     
     VelBefRot = [vx;vy];
     
     VelRot = rotationOperator(VelBefRot, -rotate(k), [0;0]);
     velx(:,il,k) = VelRot(1:end/2); vely(:,il,k) = VelRot(end/2+1:end);
        
     Xl = destandardize(tracersX(:,il,k),trans(:,k),rotate(k),rotCent(:,k),scaling(k),sortIdx(:,k));
     xlayers(:,il,k) = Xl(1:end/2);
     ylayers(:,il,k) = Xl(end/2+1:end);
  end

  figure(1);clf;
  plot(Xstand(1:end/2,k),Xstand(end/2+1:end,k),'linewidth',2)
  hold on
  plot(tracersX(1:end/2,:,k),tracersX(end/2+1:end,:,k),'k.','markersize',8)
  quiver(tracersX(1:end/2,:,k),tracersX(end/2+1:end,:,k),velx_stand,vely_stand)
  axis equal

  figure(2); clf;
  plot(X(1:end/2,k),X(end/2+1:end,k),'linewidth',2)
  hold on
  plot(xlayers(:,:,k),ylayers(:,:,k),'k.','markersize',8)
  quiver(xlayers(:,:,k),ylayers(:,:,k),velx(:,:,k),vely(:,:,k))
  axis equal

  title(k)
  pause


end

% From the tracJump build the velocity
% normalize the tracJump, it is a vector, so normalize as vinf
% then find the velocity on the layers, then denormalize both the velocity
% and the layers


% Rotate, translate tracersX to xlayers, ylayers -- ready for input form
% Also rotate velocity to velx and vely

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
