clear; clc;
addpath ../src/
addpath ../examples/
oc = curve;

load ./output/nv81N32DNN5thOrderNEWNETS_Kb2e1Dt5E5_VF35_bgFlowcouette_speed100.mat
N = numel(Xhist(:,1,1))/2; nv = numel(Xhist(1,:,1));
Xves = zeros(2*N,it*nv);
for k = 1 : it
  Xves(:,(k-1)*nv+1:k*nv) = Xhist(:,:,k);
end

% randomly choose nRandVes vesicles
nRandVes = 100;
vesIDs = randperm(it*nv,nRandVes);
X = Xves(:,vesIDs);

Nnet = 256;
nVelModes = 24;
nComp = 16;

% Fourier modes
velActiveModes = [(1:nVelModes/2)';(Nnet-nVelModes/2+1:Nnet)'];

% initialize mu,sdev. for PCA net, we do not output any of them.
muChan1 = []; sdevChan1 = []; scale = []; offset = []; outputSize = [];
netFilesFolder = './networks/NEWn256Mtimes24modesFFTNets/velPredPCAin_mode';

for imode = 1 : nVelModes
  pmode = velActiveModes(imode);
  load([netFilesFolder num2str(pmode) '_net_w1step.mat'])
  FCnets{imode} = net;
end

MVnets = FCnets; muChan_MV = muChan1; sdevChan_MV = sdevChan1;
scale_MV = scale; offset_MV = offset; MVoutSize = outputSize;

load necessaryMatFiles/pcaBasisNewest.mat % colMeans, evects
% Free-couette background flow
vBack = @(X) [X(end/2+1:end); -X(1:end/2)]...
    .*[(1-4.84./(X(1:end/2).^2+X(end/2+1:end).^2)); ...
    (1-4.84./(X(1:end/2).^2+X(end/2+1:end).^2))]/3.84;


% Evaluate MVinf exactly and approximately for 100 vesicles
MVinfTrue = zeros(2*N,nRandVes);
MVinfApp = zeros(2*N,nRandVes);
err = zeros(nRandVes,1);

op = poten(N);

for ives = 1 : nRandVes
  disp([num2str(ives) 'th vesicle out of ' num2str(nRandVes) ])
  % standardize  
  [Xstand,scaling,rotate,trans,sortIdx] = standardizationStep(X(:,ives),Nnet,oc);    
  % Prepare input for advection network
  Xinput = prepareInputForNet(Xstand,'advection',scale_MV,muChan_MV,...
      sdevChan_MV,offset_MV,colMeans,evects,nComp,[]);
  % Approximate the multiplication M*(FFTBasis)     
  Z11r = zeros(Nnet,numel(velActiveModes)); Z12r = Z11r;
  Z21r = Z11r; Z22r = Z11r;

  for k = 1 : numel(velActiveModes)
    pred = predict(MVnets{k},Xinput)';
    Z11r(:,k) = interpft(pred(1:MVoutSize/4),Nnet);
    Z21r(:,k) = interpft(pred(MVoutSize/4+1:MVoutSize/2),Nnet);
    Z12r(:,k) = interpft(pred(MVoutSize/2+1:3*MVoutSize/4),Nnet);
    Z22r(:,k) = interpft(pred(3*MVoutSize/4+1:MVoutSize),Nnet);
  end  
  
  % upsample vinf
  vinf = vBack(X(:,ives));
  vinfUp = [interpft(vinf(1:end/2),Nnet);interpft(vinf(end/2+1:end),Nnet)];
  
  % only sort points and rotate to pi/2 (no translation, no scaling)
  vinfStand = standardize(vinfUp,[0;0],rotate,1,sortIdx);
  z = vinfStand(1:end/2)+1i*vinfStand(end/2+1:end);

  zh = fft(z);
  V1 = real(zh(velActiveModes)); V2 = imag(zh(velActiveModes));
  % Compute the approximate value of the term M*vinf
  MVinfFull = [Z11r*V1+Z12r*V2; Z21r*V1+Z22r*V2];
  % Need to destandardize MVinf (take sorting and rotation back)
  MVinf = zeros(size(MVinfFull));
  MVinf([sortIdx;sortIdx+Nnet]) = MVinfFull;
  MVinf = rotationOperator(MVinf,-rotate);
  % downsample MVinf
  MVinfApp(:,ives) = [interpft(MVinf(1:end/2),N);interpft(MVinf(end/2+1:end),N)];
  
  % Compute the true one
  vesicle = capsules(X(:,ives),[],[],1,1,1);
  vesicle.setUpRate();
  
  G = op.stokesSLmatrix(vesicle);
  [~,Ten,Div] = vesicle.computeDerivs;
  
  M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
  MVinfTrue(:,ives) = M*vinf;
  
  err(ives) = sqrt(sum((MVinfTrue(1:end/2,ives)-MVinfApp(1:end/2,ives)).^2+...
      (MVinfTrue(1+end/2:end,ives)-MVinfApp(1+end/2:end,ives)).^2))./...
      sqrt(sum(MVinfTrue(:,ives).^2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xinput = prepareInputForNet(X,netType,scale,muChan,sdevChan,offset,...
    colMeans,evects,nComp,nCompRelax)
% Normalize input

if strcmp(netType,'relaxation')
  % find PCA coefficients  
  Xinput(:,1,1,1) = (X'-colMeans)*o.evects(:,1:nCompRelax);  
  Xinput(1:16,1,1,1) = scale(1)*(Xinput(1:16,1,1,1)-muChan(1))/...
    sdevChan(1)+offset(1);
  if nCompRelax > 16
    Xinput(17:32,1,1,1) = scale(2)*(Xinput(17:32,1,1,1)-muChan(2))/...
        sdevChan(2)+offset(2);
  end
else
  % find PCA coefficients  
  Xinput(:,1,1,1) = (X'-colMeans)*evects(:,1:nComp);  
  Xinput(:,1,1,1) = scale*(Xinput(:,1,1,1)-muChan)/sdevChan+offset;        
end     
end % prepareInputForNet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,scaling,rotate,trans,sortIdx] = standardizationStep(Xin,Nnet,oc)


N = numel(Xin)/2;
if Nnet ~= N
  Xin = [interpft(Xin(1:end/2),Nnet);interpft(Xin(end/2+1:end),Nnet)];    
end

X = Xin;
% Equally distribute points in arc-length
for iter = 1 : 5
  [X,~,~] = oc.redistributeArcLength(X);
end
% Fix misalignment in center and angle due to reparametrization
X = oc.alignCenterAngle(Xin,X);

% standardize angle, center, scaling and point order
[trans,rotate,scaling,sortIdx] = referenceValues(X,oc);
X = standardize(X,trans,rotate,scaling,sortIdx);
end % standardizationStep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XrotSort = standardize(X,translation,rotation,scaling,sortIdx)
N = numel(sortIdx);

% translate, rotate and scale configuration
Xrotated = scaling*rotationOperator(translateOp(X,translation),rotation);   

% now order the points
XrotSort = [Xrotated(sortIdx);Xrotated(sortIdx+N)];

end % standardize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = destandardize(XrotSort,translation,rotation,scaling,sortIdx)

N = numel(sortIdx);    
    
% change ordering back 
X = zeros(size(XrotSort));
X([sortIdx;sortIdx+N]) = XrotSort;

% scaling back
X = X/scaling;

% take rotation back
cx = mean(X(1:end/2)); cy = mean(X(end/2+1:end));
X = rotationOperator([X(1:end/2)-cx;X(end/2+1:end)-cy],-rotation);
X = [X(1:end/2)+cx; X(end/2+1:end)+cy];

% take translation back
X = translateOp(X,-translation);

end % destandardize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [translation,rotation,scaling,sortIdx] = referenceValues(Xref,oc)
N = numel(Xref)/2;

% find translation, rotation and scaling
translation = [-mean(Xref(1:end/2));-mean(Xref(end/2+1:end))];
rotation = pi/2-oc.getIncAngle2(Xref);
    
% amount of scaling
[~,~,length] = oc.geomProp(Xref);
scaling = 1/length;
    
% find the ordering of the points
Xref = scaling*rotationOperator(translateOp(Xref,translation),rotation);

firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

end % referenceValues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xrot = rotationOperator(X,theta)
% Get x-y coordinates
Xrot = zeros(size(X));
x = X(1:end/2); y = X(end/2+1:end);

% Rotated shape
xrot = (x)*cos(theta) - (y)*sin(theta);
yrot = (x)*sin(theta) + (y)*cos(theta);

Xrot(1:end/2) = xrot;
Xrot(end/2+1:end) = yrot;
end % rotationOperator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateOp(X,transXY)
Xnew = zeros(size(X));
Xnew(1:end/2) = X(1:end/2)+transXY(1);
Xnew(end/2+1:end) = X(end/2+1:end)+transXY(2);
end  % translateOp  
