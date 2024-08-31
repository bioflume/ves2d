clear; 
clc;
set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


iCalcGT = ~false;
idIC = 5;

fname = ['compareNearMeths_VesID' num2str(idIC) '.mat'];
addpath ../src/
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
N = 128;
Nup = 128;
op = poten(N,4,0);

if idIC == 1
% Parabolic
load('True_speed6000_width0.6455_FinalIC.mat')
X = [Xic(1:end/2)-mean(Xic(1:end/2)); Xic(end/2+1:end)-mean(Xic(end/2+1:end))];
elseif idIC == 2 
% from shear
load('finalShearXclose.mat')
X = [Xf(1:end/2,1)-mean(Xf(1:end/2,1)); Xf(end/2+1:end,1)-mean(Xf(end/2+1:end,1))];
elseif idIC == 3
% parabolic - symmetric
load('True_speed3000_width0.32275_FinalIC.mat')
X = [Xic(1:end/2)-mean(Xic(1:end/2)); Xic(end/2+1:end)-mean(Xic(end/2+1:end))];
elseif idIC == 4
% From taylorgreen
load('tayGreenStep140ic.mat')
X = [Xic(1:end/2,1)-mean(Xic(1:end/2,1)); Xic(end/2+1:end,1)-mean(Xic(end/2+1:end,1))];
elseif idIC == 5
% Weird shape from set-1
load('weirdICfromSet.mat')
X = Xic(:,1);
elseif idIC == 6
% Weird shape from set-2
load('weirdICfromSet2.mat')
X = Xic(:,1);
elseif idIC == 7
% Weird shape from set-2
load('weirdICfromSet2.mat')
X = Xic(:,2);
elseif idIC == 8
% Weird shape from set-2
load('weirdICfromSet2.mat')
X = Xic(:,3);
elseif idIC == 9
% Weird shape from set-2
load('weirdICfromSet2.mat')
X = Xic(:,4);
end

X = oc.upsThenFilterShape(X,512,16);
for it = 1 : 5
X = oc.redistributeArcLength(X);
end

vesicle = capsules(X,[],[],1,1,0);
h = vesicle.length/vesicle.N;

[jac,tan,~] = oc.diffProp(X);
normx = tan(end/2+1:end);
normy = -tan(1:end/2);

idP = 12;
x0 = X(idP); y0 = X(idP+N);
Ngrid = 1000; Ninit = 2*Ngrid;
% cgrid = (-cos(pi*(2*(0:Ninit-1)+1)/(2*Ninit))')+1;
cgrid = (-cos((0:Ninit-1)/(Ninit-1)*pi)'+0.50);
d2ves = h*(cgrid(1:Ngrid));
xgrid = x0 + normx(idP)*d2ves;
ygrid = y0 + normy(idP)*d2ves;

Xtra = [xgrid(:);ygrid(:)];
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;
Ntra = tracers.N;


% bendF = vesicle.tracJump(X,zeros(N,1));
bendF = zeros(2*N,1);
theta = (0:N-1)'/N * 2 * pi;
bendF(1:end/2) = sin(theta); bendF(end/2+1:end) = cos(theta);

G = op.stokesSLmatrix(vesicle);
[~,NearV2T] = vesicle.getZone(tracers,2);



%% now calculate with near-singular
kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
velTraNear = op.nearSingInt(vesicle,bendF,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false); 

velxTraNear = velTraNear(1:end/2);
velyTraNear = velTraNear(end/2+1:end);

%% now calculate without near-singular
[~,velTraNoNear] = op.exactStokesSL(vesicle,bendF,[],tracers.X,1);

velxTraNoNear = velTraNoNear(1:end/2);
velyTraNoNear = velTraNoNear(end/2+1:end);

%% now calculate with neural network
load ./shannets/near_vel_allModes_normParams/nearInterp_allModes_in_param.mat
load ./shannets/near_vel_allModes_normParams/nearInterp_allModes_out_param.mat
Nnet = 128;

%maxLayerDist = sqrt(1/Nnet); % length = 1, h = 1/Nnet;
maxLayerDist = 1/Nnet; % h
% Predictions on three layers
nlayers = 3;
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist;

% standardize input
tracersX = zeros(2*Nnet,3);

[Xstand,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(X,128);
[~,tang] = oc.diffProp(Xstand);
nx = tang(Nnet+1:2*Nnet);
ny = -tang(1:Nnet);
% 
tracersX(:,1) = Xstand;
for il = 2 : nlayers 
tracersX(:,il) = [Xstand(1:end/2)+nx*dlayer(il); Xstand(end/2+1:end)+ny*dlayer(il)];
end

% Normalize input
input_net = zeros(1,2*128,128);

for ij = 1 : 128
input_net(1,(ij-1)*2+1,:) = (Xstand(1:end/2)-in_param(1,1))/in_param(1,2);
input_net(1,2*ij,:) = (Xstand(end/2+1:end)-in_param(1,3))/in_param(1,4);
end


modes = [(0:Nnet/2-1) (-Nnet/2:-1)];
modesInUse = 128; % DO THIS WITH MORE MODES
modeList = find(abs(modes)<=modesInUse);

nv = 1;
input_conv = py.numpy.array(input_net);
[Xpredict] = pyrunfile("near_vel_allModes_predict.py","output_list",input_shape=input_conv);


velx_real = zeros(128,128,3);
vely_real = zeros(128,128,3);
velx_imag = zeros(128,128,3);
vely_imag = zeros(128,128,3);

Xpredict = double(Xpredict);
% denormalize output
for ij = 1 : numel(modeList)
  imode = modeList(ij);
  pred = Xpredict(1,(ij-1)*12+1:ij*12,:);
 
  % denormalize output

  velx_real(:,imode,1) = (pred(1,1,:)*out_param(imode,2,1))  + out_param(imode,1,1);
  velx_real(:,imode,2) = (pred(1,2,:)*out_param(imode,2,2))  + out_param(imode,1,2);
  velx_real(:,imode,3) = (pred(1,3,:)*out_param(imode,2,3))  + out_param(imode,1,3);
  vely_real(:,imode,1) = (pred(1,4,:)*out_param(imode,2,4))  + out_param(imode,1,4);
  vely_real(:,imode,2) = (pred(1,5,:)*out_param(imode,2,5))  + out_param(imode,1,5);
  vely_real(:,imode,3) = (pred(1,6,:)*out_param(imode,2,6))  + out_param(imode,1,6);

  velx_imag(:,imode,1) = (pred(1,7,:)*out_param(imode,2,7))  + out_param(imode,1,7);
  velx_imag(:,imode,2) = (pred(1,8,:)*out_param(imode,2,8))  + out_param(imode,1,8);
  velx_imag(:,imode,3) = (pred(1,9,:)*out_param(imode,2,9))  + out_param(imode,1,9);
  vely_imag(:,imode,1) = (pred(1,10,:)*out_param(imode,2,10))  + out_param(imode,1,10);
  vely_imag(:,imode,2) = (pred(1,11,:)*out_param(imode,2,11))  + out_param(imode,1,11);
  vely_imag(:,imode,3) = (pred(1,12,:)*out_param(imode,2,12))  + out_param(imode,1,12);
  
end

% outputs
% velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, trans, rotate, rotCent, scaling, sortIdx

xlayers = zeros(128,3);
ylayers = zeros(128,3);
% 
for il = 1 : 3
 Xl = destandardize(tracersX(:,il),trans,rotate,rotCent,scaling,sortIdx);
 xlayers(:,il) = Xl(1:end/2);
 ylayers(:,il) = Xl(end/2+1:end);
end


fstand = standardize(bendF,[0;0], rotate, [0;0], 1, sortIdx);
z = fstand(1:end/2) + 1i*fstand(end/2+1:end);
zh = fft(z);
fstandRe = real(zh); 
fstandIm = imag(zh);


velx = zeros(128,3); 
vely = zeros(128,3);

velx_stand = zeros(128,3);
vely_stand = zeros(128,3);
for il = 1 : 3
     velx_stand(:,il) = velx_real(:,:,il) * fstandRe + velx_imag(:,:,il)*fstandIm;
     vely_stand(:,il) = vely_real(:,:,il) * fstandRe + vely_imag(:,:,il)*fstandIm;
     
     vx = zeros(128,1);
     vy = zeros(128,1);
     
     vx(sortIdx) = velx_stand(:,il);
     vy(sortIdx) = vely_stand(:,il);

     VelBefRot = [vx;vy];
     
     VelRot = rotationOperator(VelBefRot, -rotate, [0;0]);
     velx(:,il) = VelRot(1:end/2); vely(:,il) = VelRot(end/2+1:end);
end

% 
% if sqrt(h)/2 > 1.1*h % calculate layer velocities exactly
% XlayDirect = [xlayers(:,2);xlayers(:,3);ylayers(:,2);ylayers(:,3)];
% [~,velLayer] = op.exactStokesSL(vesicle,bendF,[],XlayDirect,1);
% velx(:,2:3) = reshape(velLayer(1:end/2),Nnet,2);
% vely(:,2:3) = reshape(velLayer(end/2+1:end),Nnet,2);
% end

kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);

Xlayers = [xlayers(:,2);xlayers(:,3);ylayers(:,2);ylayers(:,3)];
traLayers.N = numel(Xlayers)/2;
traLayers.nv = 1;
traLayers.X = Xlayers;
[~,NearV2Tlayers] = vesicle.getZone(traLayers,2);

velNearLayer = op.nearSingInt(vesicle,bendF,SLP,[],NearV2Tlayers,kernel,kernelDirect,traLayers,false,false); 
velx(:,2:3) = reshape(velNearLayer(1:end/2),Nnet,2);
vely(:,2:3) = reshape(velNearLayer(end/2+1:end),Nnet,2);

%% Interpolate
iRBF = 0;
iMat = 1;

if iRBF % rbf
tRBF = tic;
Xin = [reshape(xlayers,1,3*numel(xlayers(:,1))); reshape(ylayers,1,3*numel(xlayers(:,1)))];
velXInput = reshape(velx, 1, 3*numel(xlayers(:,1))); 
velYInput = reshape(vely, 1, 3*numel(xlayers(:,1)));  

opX = rbfcreate(Xin,velXInput,'RBFFunction','linear');
opY = rbfcreate(Xin,velYInput,'RBFFunction','linear');
tRBF = toc(tRBF);
disp(['RBF build takes ' num2str(tRBF)])
end



if iMat % matlabs own functions
tMat = tic;
opX = scatteredInterpolant(xlayers(:),ylayers(:),velx(:),'linear');
opY = scatteredInterpolant(xlayers(:),ylayers(:),vely(:),'linear');
tMat = toc(tMat);
disp(['MATLAB build takes ' num2str(tMat)])
end


% Now predict
[~,velTraDirect] = op.exactStokesSL(vesicle,bendF,[],tracers.X,1);
farField = velTraDirect;
nearField = zeros(size(farField));

zone = NearV2T.zone;

J = find(zone{1}(:,1) == 1);
if numel(J) ~= 0
  [~,potTar] = op.exactStokesSL(vesicle, bendF, [], [Xtra(J,1);Xtra(J+Ntra,1)],1);
  nearField(J) = nearField(J) - potTar(1:numel(J));
  nearField(J+Ntra) = nearField(J+Ntra) - potTar(numel(J)+1:end);
  
  if iMat
  tMat = tic;    
  rbfVelX = opX(Xtra(J,1),Xtra(J+Ntra,1));
  rbfVelY = opY(Xtra(J,1),Xtra(J+Ntra,1));
  nearField(J) = nearField(J) + rbfVelX;
  nearField(J+Ntra) = nearField(J+Ntra) + rbfVelY;
  tMat = toc(tMat);
  disp(['MATLAB prediction takes ' num2str(tMat)])
  end
  % now interpolate
  if iRBF
  tRBF = tic;          
  for i = 1 : numel(J)
    
    pointsIn = [Xtra(J(i),1);Xtra(J(i)+Ntra,1)];
    % interpolate for the k2th vesicle's points near to the k1th vesicle
    rbfVelX = rbfinterp(pointsIn, opX);
    rbfVelY = rbfinterp(pointsIn, opY);    

    nearField(J(i)) = nearField(J(i)) + rbfVelX;
    nearField(J(i)+Ntra) = nearField(J(i)+Ntra) + rbfVelY; 
  end
  tRBF = toc(tRBF);
  disp(['RBF prediction takes ' num2str(tRBF)])
  end
end


velTraMLARM = farField + nearField;
velxTraMLARM = velTraMLARM(1:end/2);
velyTraMLARM = velTraMLARM(end/2+1:end);


%%
if 1
figure(3);clf;
velTraNear = sqrt(velxTraNear.^2 + velyTraNear.^2);
velTraMLARM = sqrt(velxTraMLARM.^2 + velyTraMLARM.^2);
velTraNoNear = sqrt(velxTraNoNear.^2 + velyTraNoNear.^2);

errMag = max(abs(velTraNear-velTraMLARM))/max(abs(velTraNear))
errMag2 = max(abs(velTraNear-velTraNoNear))/max(abs(velTraNear))

plot(d2ves/h,velTraNear/max(abs(velTraNear)),'k','linewidth',2)
hold on
plot(d2ves/h,velTraNoNear/max(abs(velTraNear)),'b','linewidth',2)
plot(d2ves/h,velTraMLARM/max(abs(velTraNear)),'r','linewidth',2)
axis square
xlim([-0.5 0.5])
ylim([0.9 1.2])
% xlabel('Distance to vesicle (h)')
% ylabel('Velocity magnitude')
grid on
box on
% legend('Near-singular int.','Direct evalaution','Proposed scheme')
% legend boxoff
figName = ['~/Desktop/IC' num2str(idIC) '_AlongLine.png'];
ax = gca;
exportgraphics(ax,figName,'Resolution',300)


% xlim([0 0.01])
% ylim([0.95 1.05])
% figName = ['~/Desktop/IC' num2str(idIC) '_ZoomInAlongLine.png'];
% ax = gca;
% exportgraphics(ax,figName,'Resolution',300)
end



if 0
figure(1);clf;


errX = max(abs(velxTraNear-velxTraMLARM))/max(abs(velxTraNear))
errX2 = max(abs(velxTraNear-velxTraNoNear))/max(abs(velxTraNear))

plot(d2ves/h,velxTraNear/max(abs(velxTraNear)),'k','linewidth',2)
hold on
plot(d2ves/h,velxTraNoNear/max(abs(velxTraNear)),'b','linewidth',2)
plot(d2ves/h,velxTraMLARM/max(abs(velxTraNear)),'r','linewidth',2)
axis square
xlim([-0.5 0.5])
ylim([0.9 1.2])
% xlabel('Distance to vesicle (h)')
% ylabel('Velocity magnitude')
grid on
box on
% legend('Near-singular int.','Direct evalaution','Proposed scheme')
% legend boxoff
figName = ['~/Desktop/IC' num2str(idIC) '_velx_AlongLine.png'];
ax = gca;
exportgraphics(ax,figName,'Resolution',300)

figure(2); clf;
errY = max(abs(velyTraNear-velyTraMLARM))/max(abs(velyTraNear))
errY2 = max(abs(velyTraNear-velyTraNoNear))/max(abs(velyTraNear))

plot(d2ves/h,velyTraNear/max(abs(velyTraNear)),'k','linewidth',2)
hold on
plot(d2ves/h,velyTraNoNear/max(abs(velyTraNear)),'b','linewidth',2)
plot(d2ves/h,velyTraMLARM/max(abs(velyTraNear)),'r','linewidth',2)
axis square
xlim([0 0.5])
ylim([0.9 1.2])
% xlabel('Distance to vesicle (h)')
% ylabel('Velocity magnitude')
grid on
box on
% legend('Near-singular int.','Direct evalaution','Proposed scheme')
% legend boxoff
figName = ['~/Desktop/IC' num2str(idIC) '_vely_AlongLine.png'];
ax = gca;
exportgraphics(ax,figName,'Resolution',300)


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
function LP = lagrangeInterp
% interpMap = lagrangeInterp builds the Lagrange interpolation
% matrix that takes seven function values equally distributed
% in [0,1] and returns the seven polynomial coefficients

interpMat = zeros(7);
LP(1,1) = 6.48e1;
LP(1,2) = -3.888e2;
LP(1,3) = 9.72e2;
LP(1,4) = -1.296e3;
LP(1,5) = 9.72e2;
LP(1,6) = -3.888e2;
LP(1,7) = 6.48e1;

LP(2,1) = -2.268e2;
LP(2,2) = 1.296e3;
LP(2,3) = -3.078e3;
LP(2,4) = 3.888e3;
LP(2,5) = -2.754e3;
LP(2,6) = 1.0368e3;
LP(2,7) = -1.62e2;

LP(3,1) = 3.15e2;
LP(3,2) = -1.674e3;
LP(3,3) = 3.699e3;
LP(3,4) = -4.356e3;
LP(3,5) = 2.889e3;
LP(3,6) = -1.026e3;
LP(3,7) = 1.53e2;

LP(4,1) = -2.205e2;
LP(4,2) = 1.044e3;
LP(4,3) = -2.0745e3;
LP(4,4) = 2.232e3;
LP(4,5) = -1.3815e3;
LP(4,6) = 4.68e2;
LP(4,7) = -6.75e1;

LP(5,1) = 8.12e1;
LP(5,2) = -3.132e2;
LP(5,3) = 5.265e2;
LP(5,4) = -5.08e2;
LP(5,5) = 2.97e2;
LP(5,6) = -9.72e1;
LP(5,7) = 1.37e1;

LP(6,1) = -1.47e1;
LP(6,2) = 3.6e1;
LP(6,3) = -4.5e1;
LP(6,4) = 4.0e1;
LP(6,5) = -2.25e1;
LP(6,6) = 7.2e0;
LP(6,7) = -1e0;

LP(7,1) = 1e0;
% rest of the coefficients are zero

end % lagrangeInterp