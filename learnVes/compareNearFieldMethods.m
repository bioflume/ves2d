clear; 
iCalcGT = false;
fname = 'compareNearMeths_VesID8.mat';
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
Nup = 512;
op = poten(N);
opUp = poten(Nup);

% load('True_speed6000_width0.6455_FinalIC.mat')
% X = [Xic(1:end/2)-mean(Xic(1:end/2)); Xic(end/2+1:end)-mean(Xic(end/2+1:end))];

% load('finalShearXclose.mat')
% X = [Xf(1:end/2,1)-mean(Xf(1:end/2,1)); Xf(end/2+1:end,1)-mean(Xf(end/2+1:end,1))];

% load('True_speed3000_width0.32275_FinalIC.mat')
% X = [Xic(1:end/2)-mean(Xic(1:end/2)); Xic(end/2+1:end)-mean(Xic(end/2+1:end))];

% load('tayGreenStep140ic.mat')
% X = [Xic(1:end/2,1)-mean(Xic(1:end/2,1)); Xic(end/2+1:end,1)-mean(Xic(end/2+1:end,1))];

% Weird shape from set-2
load('weirdICfromSet2.mat')
X = Xic(:,3);

X = oc.upsThenFilterShape(X,512,16);

Xup = [interpft(X(1:end/2),Nup);interpft(X(end/2+1:end),Nup)];

[xgrid, ygrid] = meshgrid(linspace(-0.25,0.25,300)',linspace(-0.25,0.25,300)');
Xtra = [xgrid(:);ygrid(:)];
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;
Ntra = tracers.N;

vesicle = capsules(X,[],[],1,1,0);
bendF = vesicle.tracJump(X,zeros(N,1));
G = op.stokesSLmatrix(vesicle);
[~,NearV2T] = vesicle.getZone(tracers,2);

if iCalcGT
vesicleUp = capsules(Xup,[],[],1,1,0);
bendFup = vesicleUp.tracJump(Xup,zeros(Nup,1));
Gup = opUp.stokesSLmatrix(vesicleUp);
% Form the ground truth
[~,NearVup2T] = vesicleUp.getZone(tracers,2);
kernel = @opUp.exactStokesSL;
kernelDirect = @opUp.exactStokesSL;
SLP = @(X) opUp.exactStokesSLdiag(vesicleUp,Gup,X);
velTraUp = opUp.nearSingInt(vesicleUp,bendFup,SLP,[],NearVup2T,kernel,kernelDirect,tracers,false,false); 
% [~,velTraUp] = op.exactStokesSL(vesicleUp,bendFup,[],tracers.X,1);
% check which points are inside
f = [ones(Nup,1);zeros(Nup,1)];
kernel = @op.exactLaplaceDL;
DLP = @(X) zeros(2*Nup,1);

Fdlp = opUp.nearSingInt(vesicleUp,f,DLP,[],NearVup2T,kernel,kernel,tracers,false,false);
idcs = abs(Fdlp) > 1e-4;

velxTraUp = velTraUp(1:end/2);
velyTraUp = velTraUp(end/2+1:end);
velxTraUp(idcs) = NaN;
velyTraUp(idcs) = NaN;

velxTraUp = reshape(velxTraUp,size(xgrid));
velyTraUp = reshape(velyTraUp,size(xgrid));

% figure(1);clf;
% plot([Xup(1:end/2);Xup(1)], [Xup(end/2+1:end);Xup(end/2+1)],'r','linewidth',2)
% hold on
% pcolor(xgrid, ygrid, sqrt(velxTraUp.^2 + velyTraUp.^2))
% colorbar
% shading interp
% axis equal

save(fname,'Xup','velxTraUp','velyTraUp','idcs')
else

load(fname);
% loading Xup velxTraUp velyTraUp idcs (which ones inside)

end

%% now calculate with near-singular
kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
velTraNear = opUp.nearSingInt(vesicle,bendF,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false); 

velxTraNear = velTraNear(1:end/2);
velyTraNear = velTraNear(end/2+1:end);
velxTraNear(idcs) = NaN;
velyTraNear(idcs) = NaN;

velxTraNear = reshape(velxTraNear,size(xgrid));
velyTraNear = reshape(velyTraNear,size(xgrid));

%% now calculate without near-singular
[~,velTraNoNear] = op.exactStokesSL(vesicle,bendF,[],tracers.X,1);

velxTraNoNear = velTraNoNear(1:end/2);
velyTraNoNear = velTraNoNear(end/2+1:end);
velxTraNoNear(idcs) = NaN;
velyTraNoNear(idcs) = NaN;

velxTraNoNear = reshape(velxTraNoNear,size(xgrid));
velyTraNoNear = reshape(velyTraNoNear,size(xgrid));

%% now calculate with neural network
Xlow = [interpft(X(1:end/2),16);interpft(X(end/2+1:end),16)];
vesicleLow = capsules(Xlow,[],[],1,1,0);
[~,NearVlow2T] = vesicleLow.getZone(tracers,2);


load ./shannets/near_vel_allModes_normParams/nearInterp_allModes_in_param.mat
load ./shannets/near_vel_allModes_normParams/nearInterp_allModes_out_param.mat
Nnet = 128;

maxLayerDist = sqrt(1/Nnet); % length = 1, h = 1/Nnet;
% Predictions on three layers
nlayers = 3;
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist;

% standardize input
tracersX = zeros(2*Nnet,3);

[Xstand,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(X,128);
[~,tang] = oc.diffProp(Xstand);
nx = tang(Nnet+1:2*Nnet);
ny = -tang(1:Nnet);

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


% layers around the vesicle j 
Xin = [reshape(xlayers,1,3*N); reshape(ylayers,1,3*N)];
velXInput = reshape(velx, 1, 3*N); 
velYInput = reshape(vely, 1, 3*N);  
  
opX = rbfcreate(Xin,velXInput,'RBFFunction','linear');
opY = rbfcreate(Xin,velYInput,'RBFFunction','linear');

% Now predict
farField = velTraNoNear;
nearField = zeros(size(farField));

zone = NearV2T.zone;

J = find(zone{1}(:,1) == 1);
if numel(J) ~= 0
  [~,potTar] = op.exactStokesSL(vesicle, bendF, [], [Xtra(J,1);Xtra(J+Ntra,1)],1);
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

velTraMLARM = farField + nearField;
velxTraMLARM = velTraMLARM(1:end/2);
velyTraMLARM = velTraMLARM(end/2+1:end);
velxTraMLARM(idcs) = NaN;
velyTraMLARM(idcs) = NaN;

velxTraMLARM = reshape(velxTraMLARM,size(xgrid));
velyTraMLARM = reshape(velyTraMLARM,size(xgrid));


%%
figure(1);clf;
plot([Xup(1:end/2);Xup(1)], [Xup(end/2+1:end);Xup(end/2+1)],'r','linewidth',2)
hold on
pcolor(xgrid, ygrid, sqrt(velxTraUp.^2 + velyTraUp.^2))
plot([Xup(1:end/2);Xup(1)], [Xup(end/2+1:end);Xup(end/2+1)],'r','linewidth',2)
c = colorbar;
c.TickLabelInterpreter = 'latex';
shading interp
axis equal
xlim([-0.25 0.25])
ylim([-0.25 0.25])
box on
set(gcf, 'renderer', 'zbuffer')

% exportgraphics(ax,figName,'Resolution',300)
title('Ground Truth')


errNear =  sqrt((velxTraUp-velxTraNear).^2 + (velyTraUp-velyTraNear).^2) ./ (1e-6+sqrt(velxTraUp.^2 + velyTraUp.^2));
figure(2);clf;
plot([Xup(1:end/2);Xup(1)], [Xup(end/2+1:end);Xup(end/2+1)],'r','linewidth',2)
hold on
pcolor(xgrid, ygrid,errNear)
plot([Xup(1:end/2);Xup(1)], [Xup(end/2+1:end);Xup(end/2+1)],'r','linewidth',2)
c = colorbar;
c.TickLabelInterpreter = 'latex';
shading interp
axis equal
box on
xlim([-0.25 0.25])
ylim([-0.25 0.25])
clim([0 1])

% exportgraphics(ax,figName,'Resolution',300)
title('Near-singular Integration')

errNoNear =  sqrt((velxTraUp-velxTraNoNear).^2 + (velyTraUp-velyTraNoNear).^2) ./ (1e-6+sqrt(velxTraUp.^2 + velyTraUp.^2));
figure(3);clf;
plot([Xup(1:end/2);Xup(1)], [Xup(end/2+1:end);Xup(end/2+1)],'r','linewidth',2)
hold on
pcolor(xgrid, ygrid, errNoNear)
plot([Xup(1:end/2);Xup(1)], [Xup(end/2+1:end);Xup(end/2+1)],'r','linewidth',2)
c = colorbar;
c.TickLabelInterpreter = 'latex';
shading interp
axis equal
box on
xlim([-0.25 0.25])
ylim([-0.25 0.25])
clim([0 1])

% exportgraphics(ax,figName,'Resolution',300)
title('No Near-singular Integration')

errMLARM =  sqrt((velxTraUp-velxTraMLARM).^2 + (velyTraUp-velyTraMLARM).^2) ./ (1e-6+sqrt(velxTraUp.^2 + velyTraUp.^2));
figure(4);clf;
plot([Xup(1:end/2);Xup(1)], [Xup(end/2+1:end);Xup(end/2+1)],'r','linewidth',2)
hold on
% plot(xlayers,ylayers,'r','linewidth',3)
pcolor(xgrid, ygrid, errMLARM)
plot([Xup(1:end/2);Xup(1)], [Xup(end/2+1:end);Xup(end/2+1)],'r','linewidth',2)
box on
c = colorbar;
c.TickLabelInterpreter = 'latex';
shading interp
axis equal
xlim([-0.25 0.25])
ylim([-0.25 0.25])
clim([0 1])
set(gcf, 'renderer', 'zbuffer')

% exportgraphics(ax,figName,'Resolution',300)
title('MLARM')


vSelf = G * bendF;
figure(5); clf;
plot(X(1:end/2),X(end/2+1:end),'r','linewidth',2)
hold on
quiver(X(1:end/2),X(end/2+1:end),vSelf(1:end/2),vSelf(end/2+1:end),'r')
quiver(xlayers(:,1),ylayers(:,1),velx(:,1),vely(:,1),'b')
% quiver(xlayers(:,2),ylayers(:,2),velx(:,2),vely(:,2),'k')
% quiver(xlayers(:,3),ylayers(:,3),velx(:,3),vely(:,3),'k')
axis equal
legend('Vesicle','True Vel.','Predicted')
box on
ax = gca;
% exportgraphics(ax,figName,'Resolution',300)


max(errNear(~isnan(errNear(:))))
min(errNear(~isnan(errNear(:))))
mean(errNear(~isnan(errNear(:))))
std(errNear(~isnan(errNear(:))))


max(errNoNear(~isnan(errNoNear(:))))
min(errNoNear(~isnan(errNoNear(:))))
mean(errNoNear(~isnan(errNoNear(:))))
std(errNoNear(~isnan(errNoNear(:))))


max(errMLARM(~isnan(errMLARM(:))))
min(errMLARM(~isnan(errMLARM(:))))
mean(errMLARM(~isnan(errMLARM(:))))
std(errMLARM(~isnan(errMLARM(:))))


errVself = sqrt((vSelf(1:end/2)-velx(:,1)).^2 + (vSelf(end/2+1:end)-vely(:,1)).^2)./(1E-4+sqrt(vSelf(1:end/2).^2 + vSelf(end/2+1:end).^2));
max(errVself)
min(errVself)
mean(errVself)
std(errVself)

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