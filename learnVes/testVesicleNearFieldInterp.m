clear; clc;

addpath ../src/

oc = curve;

Nup = 128;
N = 128;
nlayers = 3;


load testICmulti
X = Xinit(:,1);

maxLayerDist = @(h) sqrt(h);

nmodes = 128;
theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N * exp(1i*theta*ks');

op = poten(Nup);

%[Xinit,scaling,rotate,rotCent,trans,sortIdx] = standardizationStep(Xinit,128);
%[~,area,len] = oc.geomProp(Xinit);


% Generate the grid for velocity
[~,tang] = oc.diffProp(X);
% get x and y components of normal vector at each point
nx = tang(Nup+1:2*Nup);
ny = -tang(1:Nup);

% Points where velocity is calculated involve the points on vesicle
tracersX = zeros(2*Nup, nlayers-1);
tracersX(:,1) = X;

% initialize vesicle
vesicle = capsules(X, [], [], 1, 1, 0);

% Generate tracers
h = vesicle.length/vesicle.N;  % arc-length spacing
dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist(h);
for il = 2 : nlayers
tracersX(:,il) = [X(1:end/2)+nx*dlayer(il);X(end/2+1:end)+ny*dlayer(il)];
end

tracers.N = Nup;
tracers.nv = nlayers-1;
tracers.X = tracersX(:,2:nlayers);


Vinf = @(x,y) [y;-x];

% SLP 
G = op.stokesSLmatrix(vesicle);

% Get the near zone
[~,NearV2T] = vesicle.getZone(tracers,2);

% Calculate velocity on the layers
kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);

% Now reconstruct from the fourier coeff of velocity
% Choose vinf
vinf = [X(end/2+1:end); 0*X(1:end/2)];
 

% Actual function eval
selfVelactual = G*vinf;
gridVelActual = op.nearSingInt(vesicle,vinf,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);


figure(1); clf;
plot(X(1:end/2), X(end/2+1:end),'linewidth',2)
hold on
axis equal
plot(tracersX(1:end/2,:),tracersX(end/2+1:end,:),'k.','markersize',8)
quiver(tracersX(1:end/2,1),tracersX(end/2+1:end,1),selfVelactual(1:end/2),selfVelactual(end/2+1:end))
quiver(tracersX(1:end/2,2:3),tracersX(end/2+1:end,2:3),gridVelActual(1:end/2,:),gridVelActual(end/2+1:end,:))


% query points
% queryX(:,1) = [X(2)+nx(2)*sqrt(h)/3;X(2+end/2)+ny(2)*sqrt(h)/3];
% queryX(:,2) = [X(40)+nx(40)*sqrt(h)/5;X(40+end/2)+ny(4)*sqrt(h)/5];
% queryX(:,3) = [X(80)+nx(80)*sqrt(h)*6;X(80+end/2)+ny(80)*sqrt(h)*6];
% queryX(:,4) = [X(10)+nx(10)*sqrt(h)*1.00001;X(10+end/2)+ny(10)*sqrt(h)*1.00001];
% queryX(:,5) = [X(100)+nx(100)*sqrt(h)*4;X(100+end/2)+ny(100)*sqrt(h)*4];
% 1) in, 2) in, 3) out, 4) nearly out, 5) out

Nrand = 100;
queryX = zeros(2,Nrand);
randIds = randperm(128,Nrand);
randDists = 0.1 + 2*rand(Nrand,1);
queryXx = X(randIds)+nx(randIds).*sqrt(h).*randDists;
queryXy = X(randIds+end/2)+ny(randIds).*sqrt(h).*randDists;
queryX = [queryXx';queryXy'];



figure(1);
plot(queryX(1,:),queryX(2,:),'ro','markersize',8,'markerfacecolor','r')

% Check the near zone
Xlarge = [X(1:end/2)+nx*sqrt(h); X(end/2+1:end)+ny*sqrt(h)];
vesicleLarge = capsules(Xlarge, [], [], 1, 1, 0);
fCheck = [ones(vesicle.N,1);zeros(vesicle.N,1)];

tracersCheck.N = 2;
tracersCheck.nv = numel(queryX(1,:));
tracersCheck.X = queryX;

[DLP,laplaceDLPtar] = op.exactLaplaceDL(vesicleLarge,fCheck,[],queryX,1);
% laplaceDLPtar == 1 or positive -- inside
% laplaceDLPtar == 0 or negative -- outside

buffer = 1e-4;
idsIn = abs(laplaceDLPtar(1,:)) > buffer;
idsOut = abs(laplaceDLPtar(1,:)) <= buffer;

slpTracers = zeros(2,tracers.nv);
%% CALCULATE THOSE OUTSIDE USING EXACT KERNEL
[slpSelf,slpTracers(:,idsOut)] = op.exactStokesSL(vesicle,vinf,[],queryX(:,idsOut),1);

%% CALCULATE THOSE INSIDE USING INTERPOLATION
xxInput = tracersX(1:end/2,:);
yyInput = tracersX(end/2+1:end,:);
selfDensityX = selfVelactual(1:end/2);
selfDensityY = selfVelactual(end/2+1:end);
layersDensityX = gridVelActual(1:end/2,:);
layersDensityY = gridVelActual(end/2+1:end,:);

opX = rbfcreate([xxInput(:)';yyInput(:)'],[selfDensityX;layersDensityX(:)]','RBFFunction','linear');
opY = rbfcreate([xxInput(:)';yyInput(:)'],[selfDensityY;layersDensityY(:)]','RBFFunction','linear');
rbfVelX = rbfinterp(queryX(:,idsIn), opX);
rbfVelY = rbfinterp(queryX(:,idsIn), opY);
slpTracers(1,idsIn) = rbfVelX;
slpTracers(2,idsIn) = rbfVelY;

%% COMPARE THAT WITH OUR NEAR SINGULAR SCHEME
[~,NearV2T] = vesicle.getZone(tracersCheck,2);
slpNearSing = op.nearSingInt(vesicle,vinf,SLP,[],NearV2T,kernel,kernelDirect,tracersCheck,false,false);


%% Errors
errVelx = max(abs(slpNearSing(1,:)-slpTracers(1,:))./abs(slpNearSing(1,:)));
errVely = max(abs(slpNearSing(2,:)-slpTracers(2,:))./abs(slpNearSing(2,:)));

errVel = sqrt(1/Nrand*sum(((slpNearSing(1,:)-slpTracers(1,:)).^2 + (slpNearSing(2,:)-slpTracers(2,:)).^2)./...
  (slpNearSing(1,:).^2 + slpNearSing(2,:).^2)));
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
