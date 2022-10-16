classdef dnnTools

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = relaxWNN(o,Xinput,nets,nComp,evects,colMeans,scaling,...
        rotate,trans,sortIdx,muChan1,sdevChan1,offset,scale,iRes)
% Input is equally distributed in arc-length,
% output is a time step on equally distributed arc-length
N = numel(sortIdx);
% predict PCA coefficients for Xnew
if ~iRes % if not resNet, then FC net, so predict all coeffs. together
  Ypred = predict(nets,Xinput);
  % we use normalized output so take that back
  Ypred = (Ypred-offset)*sdevChan1/scale+muChan1; 
else
  Ypred = zeros(1,nComp);
  for k = 1 : nComp
    Ypred(k) = predict(nets{k},Xinput); 
  end
end % ~iRes

% reconstruct Xnew using PCA basis
Xpred = (Ypred*evects(:,1:nComp)'+colMeans)';
    
% destandardize
Xnew = o.destandardize(Xpred,trans,rotate,scaling,sortIdx,N);

end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateVinfwNN(o,Xinput,FCnets,activeModes,...
    vinf,dt,Xold,outputSize,rotate,sortIdx)
% Xinput is equally distributed in arc-length
% Xold as well. So, we add up coordinates of the same points.

N = numel(Xold)/2;

% Approximate the multiplication M*(FFTBasis)     
Z11r = zeros(N,numel(activeModes)); Z12r = Z11r;
Z21r = Z11r; Z22r = Z11r;

for k = 1 : numel(activeModes)
  pred = predict(FCnets{k},Xinput);
  Z11r(:,k) = interpft(pred(1:outputSize/4),N);
  Z21r(:,k) = interpft(pred(outputSize/4+1:outputSize/2),N);
  Z12r(:,k) = interpft(pred(outputSize/2+1:3*outputSize/4),N);
  Z22r(:,k) = interpft(pred(3*outputSize/4+1:outputSize),N);
end

% Take fft of the velocity (should be standardized velocity)
% only sort points and rotate to pi/2 (no translation, no scaling)
vinfStand = o.standardize(vinf,[0;0],rotate,1,sortIdx,N);
z = vinfStand(1:end/2)+1i*vinfStand(end/2+1:end);

zh = fft(z);
V1 = real(zh(activeModes)); V2 = imag(zh(activeModes));

% Compute the approximate value of the term M*vinf
MVinfFull = [Z11r*V1+Z12r*V2; Z21r*V1+Z22r*V2];

% Need to destandardize MVinf (take sorting and rotation back)
MVinf = zeros(size(MVinfFull));
MVinf([sortIdx;sortIdx+N]) = MVinfFull;
MVinf = o.rotationOperator(MVinf,-rotate);

% Update the position
Xnew = Xold + dt*vinf-dt*MVinf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xinput = prepareInputForNet(o,iPCA,X,colMeans,evects,nComp,...
        scale,muChan,sdevChan,offset,shapeModes)
N = numel(X)/2;
if iPCA % input to a PCA net
  % find PCA coefficients  
  Xinput(:,1,1,1) = (X'-colMeans)*evects(:,1:nComp);
else % input to a FFT net
  % FFT modes  
  z = X(1:end/2)+1i*X(end/2+1:end);
  zh = fft(z)/N;
  Xinput = zeros(2*numel(shapeModes),1,1,1);
  Xinput(1:end/2,1,1,1) = real(zh(shapeModes));
  Xinput(end/2+1:end,1,1,1) = imag(zh(shapeModes));    
end
  
% Normalize input
Xinput(:,1,1,1) = scale*(Xinput(:,1,1,1)-muChan)/sdevChan+offset;   
    
    
end % prepareInputForNet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,scaling,rotate,trans,sortIdx] = standardizationStep(o,Xin,oc)
N = numel(Xin)/2;
X = Xin;
% Equally distribute points in arc-length
for iter = 1 : 3
  [X,~,~] = oc.redistributeArcLength(X);
end
% Fix misalignment in center and angle due to reparametrization
X = oc.alignCenterAngle(Xin,X);

% standardize angle, center, scaling and point order
[trans,rotate,scaling,sortIdx] = o.referenceValues(X,oc);
X = o.standardize(X,trans,rotate,scaling,sortIdx,N);
end % standardizationStep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XrotSort = standardize(o,X,translation,rotation,scaling,sortIdx,N)

% translate, rotate and scale configuration
Xrotated = scaling*o.rotationOperator(o.translateOp(X,translation),rotation);   

% now order the points
XrotSort = [Xrotated(sortIdx);Xrotated(sortIdx+N)];

end % standardize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = destandardize(o,XrotSort,translation,rotation,scaling,sortIdx,N)

% change ordering back 
X = zeros(size(XrotSort));
X([sortIdx;sortIdx+N]) = XrotSort;

% scaling back
X = X/scaling;

% take rotation back
cx = mean(X(1:end/2)); cy = mean(X(end/2+1:end));
X = o.rotationOperator([X(1:end/2)-cx;X(end/2+1:end)-cy],-rotation);
X = [X(1:end/2)+cx; X(end/2+1:end)+cy];

% take translation back
X = o.translateOp(X,-translation);

end % destandardize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [translation,rotation,scaling,sortIdx] = referenceValues(o,Xref,oc)

N = numel(Xref)/2;

% find translation, rotation and scaling
translation = [-mean(Xref(1:end/2));-mean(Xref(end/2+1:end))];
rotation = pi/2-oc.getIncAngle2(Xref);
    
% amount of scaling
[~,~,length] = oc.geomProp(Xref);
scaling = 1/length;
    
% find the ordering of the points
Xref = scaling*o.rotationOperator(o.translateOp(Xref,translation),rotation);

firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

end % referenceValues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xrot = rotationOperator(o,X,theta)
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
function Xnew = translateOp(o,X,transXY)
Xnew = zeros(size(X));
Xnew(1:end/2) = X(1:end/2)+transXY(1);
Xnew(end/2+1:end) = X(end/2+1:end)+transXY(2);
end  % translateOp  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,errAL] = exactlySolve(o,Xold,vinf,dt,area0,len0,oc,op,kappa)
% Semi-implicit time stepping w/o splitting
vesicle = capsules(Xold,[],[],kappa,1,1);
vesicle.setUpRate();

% first reparameterize 
Xold = oc.reparametrize(Xold, [], 6, 50);

% Take a time step
Xnew = o.relaxExactSolve(vesicle,vinf(Xold),dt,Xold,op);

% Area-length correction
[Xnew2,ifail] = oc.correctAreaAndLength3(Xnew,area0,len0,32);
if ifail
  disp('Area-length correction failed!')    
end
Xnew = oc.alignCenterAngle(Xnew,Xnew2);

% Error in area-length
[~,area,len] = oc.geomProp(Xnew);
errAL = max(abs(area-area0)/area0,abs(len-len0)/len0);
end
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
vesicle = capsules(Xold,[],[],1,1,1);
vesicle.setUpRate();
 
G = op.stokesSLmatrix(vesicle);
[~,Ten,Div] = vesicle.computeDerivs;

M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
Xnew = Xold + dt*(eye(2*vesicle.N)-M)*vinf;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vinf = setBgFlow(o,bgFlow,speed)

vinf = @(X) zeros(size(X));      
if strcmp(bgFlow,'relax')
  vinf = @(X) zeros(size(X));  % relaxation
elseif strcmp(bgFlow,'shear') 
  vinf = @(X) speed*[X(end/2+1:end);zeros(size(X(1:end/2)))]; 
elseif strcmp(bgFlow,'tayGreen')
  vinf = @(X) speed*[sin(X(1:end/2)).*cos(X(end/2+1:end));-...
    cos(X(1:end/2)).*sin(X(end/2+1:end))]; % Taylor-Green
elseif strcmp(bgFlow,'parabolic')
  vinf = @(X) [speed*(1-(X(end/2+1:end)/0.2).^2);...
      zeros(size(X(1:end/2)))];
elseif strcmp(bgFlow,'rotation')
  vinf = @(X) [-sin(atan2(X(end/2+1:end),X(1:end/2)))./sqrt(X(1:end/2).^2+X(end/2+1:end).^2);...
    cos(atan2(X(end/2+1:end),X(1:end/2)))./sqrt(X(1:end/2).^2+X(end/2+1:end).^2)]*speed;
end
    
end % setBgFlow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % methods

end % dnnTools
