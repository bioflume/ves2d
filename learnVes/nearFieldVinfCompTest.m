clear; clc;

addpath ../src/

oc = curve;

Nup = 128;
N = 128;
nlayers = 3;

% load n128Dt1e-05RelaxMirrdDataSet.mat
% k = 11212;
% Xinit =  [interpft(XstandStore(1:end/2,k),Nup); interpft(XstandStore(end/2+1:end,k),Nup)];
load testIC

maxLayerDist = @(h) sqrt(h);

nmodes = 128;
theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N * exp(1i*theta*ks');

% activeModes = [(1:nmodes/2)';(N-nmodes/2+1:N)'];

op = poten(Nup);

[Xinit,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(Xinit,128);
[~,area,len] = oc.geomProp(Xinit);


% Generate the grid for velocity
[~,tang] = oc.diffProp(Xinit);
% get x and y components of normal vector at each point
nx = tang(Nup+1:2*Nup);
ny = -tang(1:Nup);

% Points where velocity is calculated involve the points on vesicle
tracersX = zeros(2*Nup, nlayers-1);
tracersX(:,1) = Xinit;

% initialize vesicle
vesicle = capsules(Xinit, [], [], 1, 1, 1);
vesicle.setUpRate();

% Generate tracers
h = vesicle.length/vesicle.N;  % arc-length spacing
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist(h);
for il = 2 : nlayers
tracersX(:,il) = [Xinit(1:end/2)+nx*dlayer(il);Xinit(end/2+1:end)+ny*dlayer(il)];
end

tracers.N = Nup;
tracers.nv = nlayers-1;
tracers.X = tracersX(:,2:nlayers);


Vinf = @(x,y) [y;zeros(size(x))];

% SLP 
G = op.stokesSLmatrix(vesicle);

% Get the near zone
[~,NearV2T] = vesicle.getZone(tracers,2);

% CALCULATE VELOCITY ON THE FOURIER MODES FIRST
VelOnGridModesReal = zeros(2*Nup,nlayers-1,nmodes);
selfVelModesReal = zeros(2*Nup,nmodes);

VelOnGridModesImag = zeros(2*Nup,nlayers-1,nmodes);
selfVelModesImag = zeros(2*Nup,nmodes);


Br = real(basis); Bi = imag(basis);

kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);

for imode = 1 : nmodes
  forRealVels = [Br(:,imode); Bi(:,imode)];
  forImagVels = [-Bi(:,imode); Br(:,imode)];

  VelOnGridModesReal(:,:,imode) = op.nearSingInt(vesicle,forRealVels,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);
  selfVelModesReal(:,imode) = G*forRealVels;

  VelOnGridModesImag(:,:,imode) = op.nearSingInt(vesicle,forImagVels,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);
  selfVelModesImag(:,imode) = G*forImagVels;
end

% Now reconstruct from the fourier coeff of velocity
% Choose vinf
vinf = [Xinit(end/2+1:end); 0*Xinit(1:end/2)];
 
vz = vinf(1:end/2)+1i*vinf(end/2+1:end);
coeffs = fft(vz);
Vre = real(coeffs); Vim = imag(coeffs);

% Reconstruct the self-velocity
selfVelX_recon = selfVelModesReal(1:end/2,:)*Vre + selfVelModesImag(1:end/2,:)*Vim;
selfVelY_recon = selfVelModesReal(end/2+1:end,:)*Vre + selfVelModesImag(end/2+1:end,:)*Vim;

% Actual function eval
selfVelactual = G*vinf;

errSelfVel = norm(selfVelactual-[selfVelX_recon;selfVelY_recon]);


% Loop over layers
gridVelX_recon = zeros(128,2);
gridVelY_recon = zeros(128,2);
for il = 1 : 2
  velModesReal = reshape(VelOnGridModesReal(:,il,:),256,128);
  velModesImag = reshape(VelOnGridModesImag(:,il,:),256,128);

  gridVelX_recon(:,il) = velModesReal(1:end/2,:)*Vre + velModesImag(1:end/2,:)*Vim;
  gridVelY_recon(:,il) = velModesReal(end/2+1:end,:)*Vre + velModesImag(end/2+1:end,:)*Vim;

end

gridVelActual = op.nearSingInt(vesicle,vinf,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);

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
