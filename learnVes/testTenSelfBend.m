clear; clc;
addpath ../src/
addpath ../examples/
oc = curve;

kappa = 1e-2;

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
nTenModes = 32;
nComp = 16;

% Fourier modes
activeModes = [(1:nTenModes/2)';(Nnet-nTenModes/2+1:Nnet)'];

% initialize mu,sdev. for PCA net, we do not output any of them.
muChan1 = []; sdevChan1 = []; scale = []; offset = []; outputSize = [];
load('./networks/NEWfcPCAtenMaTonBendN256_32modes_FCNet_w1step.mat')

tenBendNets = net; muChan_tenBend = muChan1; 
sdevChan_tenBend = sdevChan1; scale_tenBend = scale; 
offset_tenBend = offset; 

load necessaryMatFiles/pcaBasisNewest.mat % colMeans, evects


% Evaluate MVinf exactly and approximately for 100 vesicles
selfBendTrue = zeros(N,nRandVes);
selfBendApp = zeros(N,nRandVes);
err = zeros(nRandVes,1);

op = poten(N);

for ives = 1 : nRandVes
  disp([num2str(ives) 'th vesicle out of ' num2str(nRandVes) ])
  % standardize  
  [Xstand,scaling,rotate,trans,sortIdx] = standardizationStep(X(:,ives),Nnet,oc);    
  % Prepare input for advection network
  Xinput = prepareInputForNet(Xstand,'tensionOnFourier',scale_tenBend,muChan_tenBend,...
      sdevChan_tenBend,offset_tenBend,colMeans,evects,nComp,[]);
  % Approximate the multiplication M*(FFTBasis)     
  realImagActive = predict(tenBendNets,Xinput)';
  % initialize full output
  zFull = zeros(Nnet,1);
  zFull(activeModes) = realImagActive(1:end/2)+1i*realImagActive(end/2+1:end);

  % ifft and find upsampled output
  upsampOutStand = real(ifft(zFull)*Nnet); 
  
  upsampOut = zeros(Nnet,1);  
  upsampOut(sortIdx) = kappa*upsampOutStand/scaling^2;
  % downsample
  selfBendApp(:,ives) = interpft(upsampOut,N);
  
  % Compute the true one
  vesicle = capsules(X(:,ives),[],[],kappa,1,1);
  vesicle.setUpRate();
  
  G = op.stokesSLmatrix(vesicle);
  [Ben,Ten,Div] = vesicle.computeDerivs;
  M = ((Div*G*Ten)\eye(vesicle.N))*Div;
  selfBendTrue(:,ives) = kappa*M*G*(-Ben)*X(:,ives);
  
  err(ives) = sqrt(sum((selfBendTrue(:,ives)-selfBendApp(:,ives)).^2))./...
      sqrt(sum(selfBendTrue(:,ives).^2));
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
