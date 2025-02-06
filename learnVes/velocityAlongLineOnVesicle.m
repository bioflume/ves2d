clear; 
clc;
set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


iCalcGT = ~false;
idIC = 2;

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
Nup1 = ceil(sqrt(N))*N;
Nup2 = ceil(sqrt(32))*32;

op128 = poten(N,4,0);
opNup1 = poten(Nup1,4,0);
opNup2 = poten(Nup2,4,0);
op32 = poten(32,4,0);

% opUp = poten(Nup,4,0);

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
XOrig = X;
for it = 1 : 5
  X = oc.redistributeArcLength(X);
end
X = oc.alignCenterAngle(XOrig,X);

vesicle = capsules(X,[],[],1,1,0);
h = vesicle.length/vesicle.N;

% Xup = [interpft(X(1:end/2),Nup);interpft(X(end/2+1:end),Nup)];
[jac,tan,~] = oc.diffProp(X);
normx = tan(end/2+1:end);
normy = -tan(1:end/2);
% 
% idP = 104;
% x0 = X(idP); y0 = X(idP+N);
% Ngrid = 2000; 
% cgrid = linspace(-2,2,Ngrid);
% d2ves = h*(cgrid);
% xgrid104 = x0 + normx(idP)*d2ves;
% ygrid104 = y0 + normy(idP)*d2ves;
% 
% idP = 118;
% x0 = X(idP); y0 = X(idP+N);
% Ngrid = 2000; 
% cgrid = linspace(-2,2,Ngrid);
% d2ves = h*(cgrid);
% xgrid118 = x0 + normx(idP)*d2ves;
% ygrid118 = y0 + normy(idP)*d2ves;

idP = 118;
x0 = X(idP); y0 = X(idP+N);
Ngrid = 2000; 
cgrid = linspace(-2,2,Ngrid);
d2ves = h*(cgrid);
xgrid = x0 + normx(idP)*d2ves;
ygrid = y0 + normy(idP)*d2ves;

% 
% figure(3);clf;
% plot([X(1:end/2);X(1)],[X(end/2+1:end);X(end/2+1)],'k','linewidth',2)
% axis equal
% hold on
% plot(X([12;104;118]),X([12;104;118]+128),'ro','markersize',8,'markerfacecolor','r')
% plot(xgrid, ygrid, 'r','linewidth',2)
% plot(xgrid104, ygrid104, 'r','linewidth',2)
% plot(xgrid118, ygrid118, 'r','linewidth',2)
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'ztick',[]);
% 
% set(gca,'xcolor','w');
% set(gca,'ycolor','w');
% set(gca,'zcolor','w');
% box on
% set(gca,'visible','off')
% 
% ax = gca;
% exportgraphics(ax,['~/Desktop/IC' num2str(idIC) '.png'],'Resolution',300)
% pause

% [xgrid, ygrid] = meshgrid(linspace(-0.25,0.25,500)',linspace(-0.25,0.25,500)');
Xtra = [xgrid(:);ygrid(:)];
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;
Ntra = tracers.N;


% bendF = vesicle.tracJump(X,zeros(N,1));
bendF = zeros(2*N,1);
theta = (0:N-1)'/N * 2 * pi;
bendF(1:end/2) = sin(theta); bendF(end/2+1:end) = cos(theta);

G = op128.stokesSLmatrix(vesicle);
[~,NearV2T] = vesicle.getZone(tracers,2);

%% now calculate with near-singular
kernel = @op128.exactStokesSL;
kernelDirect = @op128.exactStokesSL;
SLP = @(X) op128.exactStokesSLdiag(vesicle,G,X);
velTraNear = op128.nearSingInt(vesicle,bendF,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false); 

velxTraNear = velTraNear(1:end/2);
velyTraNear = velTraNear(end/2+1:end);


%% now calculate with neural network
% Xlow = [interpft(X(1:end/2),16);interpft(X(end/2+1:end),16)];
% vesicleLow = capsules(Xlow,[],[],1,1,0);
% [~,NearVlow2T] = vesicleLow.getZone(tracers,2);
% 

load ./shannets/nearInterp_128modes_disth_params.mat

Nnet = 128;


% Normalize input
input_net = zeros(1,2*128,128);

[Xstand,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(X,128);

for ij = 1 : 128
input_net(1,(ij-1)*2+1,:) = (Xstand(1:end/2)-in_param(1,1))/in_param(1,2);
input_net(1,2*ij,:) = (Xstand(end/2+1:end)-in_param(1,3))/in_param(1,4);
end

modes = [(0:Nnet/2-1) (-Nnet/2:-1)];
modesInUse = 128; % DO THIS WITH MORE MODES
modeList = find(abs(modes)<=modesInUse);

input_conv = py.numpy.array(input_net);
[Xpredict] = pyrunfile("near_vel_allModesAth_predict.py","output_list",input_shape=input_conv);

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

% xlayers = zeros(128,3);
% ylayers = zeros(128,3);
% 
% for il = 1 : 3
%  Xl = destandardize(tracersX(:,il),trans,rotate,rotCent,scaling,sortIdx);
%  xlayers(:,il) = Xl(1:end/2);
%  ylayers(:,il) = Xl(end/2+1:end);
% end


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


maxLayerDist = 1/Nnet; % 2h
% Predictions on three layers
nlayers = 3;
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist;

% standardize input
tracersX = zeros(2*Nnet,nlayers);

[Xstand,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(X,128);
[~,tang] = oc.diffProp(Xstand);
nx = tang(Nnet+1:2*Nnet);
ny = -tang(1:Nnet);
% 
tracersX(:,1) = Xstand;
for il = 2 : nlayers 
tracersX(:,il) = [Xstand(1:end/2)+nx*dlayer(il); Xstand(end/2+1:end)+ny*dlayer(il)];
end



% outputs
% velx_real, vely_real, velx_imag, vely_imag, xlayers, ylayers, trans, rotate, rotCent, scaling, sortIdx

xlayers = zeros(128,nlayers);
ylayers = zeros(128,nlayers);
% 
for il = 1 : nlayers
 Xl = destandardize(tracersX(:,il),trans,rotate,rotCent,scaling,sortIdx);
 xlayers(:,il) = Xl(1:end/2);
 ylayers(:,il) = Xl(end/2+1:end);
end


% kernel = @op.exactStokesSL;
% kernelDirect = @op.exactStokesSL;
% SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
% 
% Xlayers = [xlayers(:,4);xlayers(:,5);ylayers(:,4);ylayers(:,5)];
% traLayers.N = numel(Xlayers)/2;
% traLayers.nv = 1;
% traLayers.X = Xlayers;
% [~,NearV2Tlayers] = vesicle.getZone(traLayers,2);
% 
% velNearLayer = op.nearSingInt(vesicle,bendF,SLP,[],NearV2Tlayers,kernel,kernelDirect,traLayers,false,false); 
% velx(:,4:5) = reshape(velNearLayer(1:end/2),Nnet,2);
% vely(:,4:5) = reshape(velNearLayer(end/2+1:end),Nnet,2);

% Now predict
Xup1 = [interpft(X(1:end/2),Nup1);interpft(X(end/2+1:end),Nup1)];
vesicleUp1 = capsules(Xup1,[],[],1,1,0);
bendFup1 = [interpft(bendF(1:end/2),Nup1);interpft(bendF(end/2+1:end),Nup1)];

[~,velTraDirect] = opNup1.exactStokesSL(vesicleUp1,bendFup1,[],tracers.X,1);

velDirectN128 = velTraDirect;
farField = velTraDirect;

dist = NearV2T.dist;
nearest = NearV2T.nearest;
icp = NearV2T.icp;
argnear = NearV2T.argnear;
interpMat = lagrangeInterp;
interpOrder = size(interpMat,1);
p = ceil((interpOrder+1)/2);
beta = 1.1;
vself = [velx(:,1);vely(:,1)]; % velocity on vesicle itself
hves = vesicle.length/vesicle.N;

Xin = [reshape(xlayers,1,nlayers*numel(xlayers(:,1))); reshape(ylayers,1,nlayers*numel(xlayers(:,1)))];
velXInput = reshape(velx, 1, nlayers*numel(xlayers(:,1))); 
velYInput = reshape(vely, 1, nlayers*numel(xlayers(:,1)));  

opX = rbfcreate(Xin,velXInput,'RBFFunction','linear');
opY = rbfcreate(Xin,velYInput,'RBFFunction','linear');

nearField = zeros(size(farField));

zone = NearV2T.zone;

J = find(zone{1}(:,1) == 1);
if numel(J) ~= 0
  [~,potTar] = opNup1.exactStokesSL(vesicleUp1, bendFup1, [], [Xtra(J,1);Xtra(J+Ntra,1)],1);
  nearField(J) = nearField(J) - potTar(1:numel(J));
  nearField(J+Ntra) = nearField(J+Ntra) - potTar(numel(J)+1:end);
  
  % now interpolate     
  for i = 1 : numel(J)
    pointsIn = [Xtra(J(i),1);Xtra(J(i)+Ntra,1)];
    % interpolate for the k2th vesicle's points near to the k1th vesicle
    rbfVelX = rbfinterp(pointsIn, opX);
    rbfVelY = rbfinterp(pointsIn, opY);    

    nearField(J(i)) = nearField(J(i)) + rbfVelX;
    nearField(J(i)+Ntra) = nearField(J(i)+Ntra) + rbfVelY; 
  end
  
end

velTraMLARMRBF = farField + nearField;
velxTraMLARMRBF = velTraMLARMRBF(1:end/2);
velyTraMLARMRBF = velTraMLARMRBF(end/2+1:end);


%% now calculate with neural network, N = 32
Xlow = [interpft(X(1:end/2),32);interpft(X(end/2+1:end),32)];
vesicleLow = capsules(Xlow,[],[],1,1,0);


load ./shannets/nearInterp_32modes_in_param.mat
load ./shannets/nearInterp_32modes_out_param.mat
Nnet = 32;

bendF32 = zeros(2*Nnet,1);
theta = (0:Nnet-1)'/Nnet * 2 * pi;
bendF32(1:end/2) = sin(theta); bendF32(end/2+1:end) = cos(theta);

Xup2 = [interpft(Xlow(1:end/2),Nup2);interpft(Xlow(end/2+1:end),Nup2)];
vesicleUp2 = capsules(Xup2,[],[],1,1,0);
bendFup2 = [interpft(bendF32(1:end/2),Nup2);interpft(bendF32(end/2+1:end),Nup2)];


G32 = op32.stokesSLmatrix(vesicleLow);
[~,NearV322T] = vesicleLow.getZone(tracers,2);

[Xstand,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(X,32);

% Normalize input
input_net = zeros(1,2,32);

input_net(1,1,:) = (Xstand(1:end/2)-in_param(1,1))/in_param(1,2);
input_net(1,2,:) = (Xstand(end/2+1:end)-in_param(1,3))/in_param(1,4);

modes = [(0:Nnet/2-1) (-Nnet/2:-1)];
modesInUse = 32; % DO THIS WITH MORE MODES
modeList = find(abs(modes)<=modesInUse);

nv = 1;
input_conv = py.numpy.array(input_net);
[Xpredict] = pyrunfile("32modes_near_vel_predict.py","output_list",input_shape=input_conv,num_ves=py.int(nv),modesInUse=py.int(modesInUse));


velx_real = zeros(Nnet,Nnet,3);
vely_real = zeros(Nnet,Nnet,3);
velx_imag = zeros(Nnet,Nnet,3);
vely_imag = zeros(Nnet,Nnet,3);

% denormalize output
for ij = 1 : numel(modeList)
  imode = modeList(ij);
  pred = double(Xpredict{ij});
 
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

fstand = standardize(bendF32,[0;0], rotate, [0;0], 1, sortIdx);
z = fstand(1:end/2) + 1i*fstand(end/2+1:end);
zh = fft(z);
fstandRe = real(zh); 
fstandIm = imag(zh);

velx32 = zeros(Nnet,3); 
vely32 = zeros(Nnet,3);


velx32_stand = zeros(Nnet,3);
vely32_stand = zeros(Nnet,3);
for il = 1 : 3
     velx32_stand(:,il) = velx_real(:,:,il) * fstandRe + velx_imag(:,:,il)*fstandIm;
     vely32_stand(:,il) = vely_real(:,:,il) * fstandRe + vely_imag(:,:,il)*fstandIm;
     
     vx = zeros(Nnet,1);
     vy = zeros(Nnet,1);
     
     vx(sortIdx) = velx32_stand(:,il);
     vy(sortIdx) = vely32_stand(:,il);

     VelBefRot = [vx;vy];
     
     VelRot = rotationOperator(VelBefRot, -rotate, [0;0]);
     velx32(:,il) = VelRot(1:end/2); vely32(:,il) = VelRot(end/2+1:end);
end


maxLayerDist = 1/Nnet; % 2h

% Predictions on three layers
nlayers = 3;
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist;

% standardize input
tracersX = zeros(2*Nnet,nlayers);

[Xstand,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(Xlow,Nnet);
[~,tang] = oc.diffProp(Xstand);
nx = tang(Nnet+1:2*Nnet);
ny = -tang(1:Nnet);
% 
tracersX(:,1) = Xstand;
for il = 2 : nlayers 
tracersX(:,il) = [Xstand(1:end/2)+nx*dlayer(il); Xstand(end/2+1:end)+ny*dlayer(il)];
end

% outputs

xlayers32 = zeros(Nnet,nlayers);
ylayers32 = zeros(Nnet,nlayers);

% 
for il = 1 : nlayers
 Xl = destandardize(tracersX(:,il),trans,rotate,rotCent,scaling,sortIdx);
 xlayers32(:,il) = Xl(1:end/2);
 ylayers32(:,il) = Xl(end/2+1:end);
end


% kernel = @op.exactStokesSL;
% kernelDirect = @op.exactStokesSL;
% SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
% 
% Xlayers = [xlayers32_128(:,4);xlayers32_128(:,5);ylayers32_128(:,4);ylayers32_128(:,5)];
% traLayers.N = numel(Xlayers)/2;
% traLayers.nv = 1;
% traLayers.X = Xlayers;
% [~,NearV2Tlayers] = vesicle.getZone(traLayers,2);
% 
% velNearLayer = op.nearSingInt(vesicle,bendF,SLP,[],NearV2Tlayers,kernel,kernelDirect,traLayers,false,false); 
% velx32_128(:,4:5) = reshape(velNearLayer(1:end/2),128,2);
% vely32_128(:,4:5) = reshape(velNearLayer(end/2+1:end),128,2);

Xin = [reshape(xlayers32,1,nlayers*numel(xlayers32(:,1))); reshape(ylayers32,1,nlayers*numel(xlayers32(:,1)))];
velXInput = reshape(velx32, 1, nlayers*numel(xlayers32(:,1))); 
velYInput = reshape(vely32, 1, nlayers*numel(xlayers32(:,1)));  

opX = rbfcreate(Xin,velXInput,'RBFFunction','linear');
opY = rbfcreate(Xin,velYInput,'RBFFunction','linear');


% Now predict
[~,velTraDirect] = opNup2.exactStokesSL(vesicleUp2,bendFup2,[],tracers.X,1);
velDirectN32 = velTraDirect;

farField = velTraDirect;
nearField = zeros(size(farField));

dist = NearV322T.dist;
nearest = NearV322T.nearest;
icp = NearV322T.icp;
argnear = NearV322T.argnear;
interpMat = lagrangeInterp;
interpOrder = size(interpMat,1);
p = ceil((interpOrder+1)/2);
beta = 1.1;
vself = [velx32(:,1);vely32(:,1)]; % velocity on vesicle itself
hves = vesicleLow.length/vesicleLow.N;

zone = NearV322T.zone;

J = find(zone{1}(:,1) == 1);
if numel(J) ~= 0
  [~,potTar] = opNup2.exactStokesSL(vesicleUp2, bendFup2, [], [Xtra(J,1);Xtra(J+Ntra,1)],1);
  nearField(J) = nearField(J) - potTar(1:numel(J));
  nearField(J+Ntra) = nearField(J+Ntra) - potTar(numel(J)+1:end);
  
  % now interpolate
  
  for i = 1 : numel(J)
    
    pointsIn = [Xtra(J(i),1);Xtra(J(i)+Ntra,1)];
    % interpolate for the k2th vesicle's points near to the k1th vesicle
    rbfVelX = rbfinterp(pointsIn, opX);
    rbfVelY = rbfinterp(pointsIn, opY);    

    nearField(J(i)) = nearField(J(i)) + rbfVelX;
    nearField(J(i)+Ntra) = nearField(J(i)+Ntra) + rbfVelY; 
  end
end

velTraMLARMRBF32 = farField + nearField;
velxTraMLARMRBF32 = velTraMLARMRBF32(1:end/2);
velyTraMLARMRBF32 = velTraMLARMRBF32(end/2+1:end);
%%
if 1

velTraNear = sqrt(velxTraNear.^2 + velyTraNear.^2);
velTraMLARMRBF = sqrt(velxTraMLARMRBF.^2 + velyTraMLARMRBF.^2);
velTraMLARMRBF32 = sqrt(velxTraMLARMRBF32.^2 + velyTraMLARMRBF32.^2);
velTraDirectN128 = sqrt(velDirectN128(1:end/2).^2 + velDirectN128(end/2+1:end).^2);
velTraDirectN32 = sqrt(velDirectN32(1:end/2).^2 + velDirectN32(end/2+1:end).^2);

figure(1);clf; % N = 128
plot(d2ves/h,velTraNear/max(abs(velTraNear)),'Color',[94 60 153]/255,'linewidth',2)
hold on
% plot(d2ves/h,velTraDirectN128/max(abs(velTraNear)),'Color',[178 171 210]/255,'linewidth',2)
plot(d2ves/h,velTraMLARMRBF/max(abs(velTraNear)),'--','Color',[253 184 99]/255,'linewidth',2)

% axis square
% xlim([-2 2])
% ylim([0 1])
% xticks([-2, -1, 0, 1, 2])
% yticks([0.9, 1, 1.1, 1.2])
% xlabel('Distance to vesicle (h)')
% ylabel('Velocity magnitude')
% grid on
% box on
% figName = ['~/Desktop/IdP' num2str(idP) '_3Layers.png'];
% ax = gca;
% exportgraphics(ax,figName,'Resolution',300)
% 
% legend('Ground truth','Ves-Net (N = 128)','Ves-Net (N = 32)','location','southeast')
% legend boxoff


%  
figure(2); clf; % N = 32
plot(d2ves/h,velTraNear/max(abs(velTraNear)),'Color',[94 60 153]/255,'linewidth',2)
hold on
% plot(d2ves/h,velTraDirectN32/max(abs(velTraNear)),'Color',[178 171 210]/255,'linewidth',2)
plot(d2ves/h,velTraMLARMRBF32/max(abs(velTraNear)),'Color',[253 184 99]/255,'linewidth',2)


% plot(X(1:end/2),X(end/2+1:end),'k')
% axis equal
% hold on
% quiver(X(1:end/2),X(end/2+1:end),vselfTrue(1:end/2),vselfTrue(end/2+1:end),'b')
% quiver(X(1:end/2),X(end/2+1:end),vself(1:end/2),vself(end/2+1:end),'r')
% legend('Vesicle','Ground truth','Ves-Net')
% legend boxoff

% xlim([0 0.01])
% ylim([0.95 1.05])
% figName = ['~/Desktop/IC' num2str(idIC) '_ZoomInAlongLine.png'];
% ax = gca;
% exportgraphics(ax,figName,'Resolution',300)
end



if 0
figure(2);clf;
errX = max(abs(velxTraNear-velxTraMLARMRBF))/max(abs(velTraNear))
errX2 = max(abs(velxTraNear-velxTraMLARMRBF32))/max(abs(velTraNear))

plot(d2ves/h,velxTraNear/max(abs(velTraNear)),'Color',[.5 .5 .5],'linewidth',2)
velDirectAllPoints
hold on
plot(d2ves/h,velxTraMLARMRBF/max(abs(velTraNear)),'--','Color',[69,117,180]/255,'linewidth',2)
% plot(d2ves/h,velxTraMLARMRBF32/max(abs(velTraNear)),'Color',[215,48,39]/255,'linewidth',2)
axis square
xlim([0 2])
% ylim([0 1])
xticks([0, 1, 2])
% yticks([0.9, 1, 1.1, 1.2])
% xlabel('Distance to vesicle (h)')
% ylabel('Velocity magnitude')
grid on
box on
% legend('Ground truth','Ves-Net (N = 128)','Ves-Net (N = 32)')
% legend boxoff
figName = ['~/Desktop/IC' num2str(idIC) '_velx_AlongLine.png'];
ax = gca;
exportgraphics(ax,figName,'Resolution',300)

figure(3); clf;
errY = max(abs(velyTraNear-velyTraMLARMRBF))/max(abs(velyTraNear))
errY2 = max(abs(velyTraNear-velyTraMLARMRBF32))/max(abs(velyTraNear))


plot(d2ves/h,velyTraNear/max(abs(velTraNear)),'Color',[.5 .5 .5],'linewidth',2)
hold on
plot(d2ves/h,velyTraMLARMRBF/max(abs(velTraNear)),'--','Color',[69,117,180]/255,'linewidth',2)
% plot(d2ves/h,velyTraMLARMRBF32/max(abs(velTraNear)),'Color',[215,48,39]/255,'linewidth',2)
axis square
xlim([0 2])
% ylim([0 1])
xticks([0, 1, 2])
% yticks([0.9, 1, 1.1, 1.2])
% xlabel('Distance to vesicle (h)')
% ylabel('Velocity magnitude')
grid on
box on
% legend('Near-singular int.','Direct evalaution','Proposed scheme')
% legend boxoff
figName = ['~/Desktop/IC' num2str(idIC) '_vely_AlongLine.png'];
ax = gca;
exportgraphics(ax,figName,'Resolution',300)


figure(4); clf;
velSelfMag = sqrt(velx(:,1).^2 + vely(:,1).^2);
velSelf32Mag = sqrt(velx32_128(:,1).^2 + vely32_128(:,1).^2);
errSelfX = max(abs(velx(:,1)-velx32_128(:,1)))./max(abs(velSelfMag))
errSelfY = max(abs(vely(:,1)-vely32_128(:,1)))./max(abs(velSelfMag))
errSelf = max(abs(velSelfMag-velSelf32Mag))./max(abs(velSelfMag))
plot(xlayers(:,1),ylayers(:,1),'Color',[.5 .5 .5],'linewidth',2)
hold on
quiver(xlayers(:,1),ylayers(:,1),velx(:,1),vely(:,1),'r')
quiver(xlayers(:,1),ylayers(:,1),velx32_128(:,1),vely32_128(:,1),'b')
axis equal
grid on
box on



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