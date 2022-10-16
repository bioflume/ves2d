%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD ROUTINES FROM poten.m  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [x,w] = regularWeights(o,a)
% [x,w] = regularWeights(a) gets quadrature rules for regular functions
% (Alpert's quadrature rules).  These are needed for integrating away
% from singularities This code is not currently used.  It is not needed
% to integrate functions with logarithmic singularities and the
% double-layer potential is accurate with the trapezoid rule

par = [2;3;4;5;7;10];
for i = 1:length(par)
  if par(i)==a
    key = i;
  end
end

lenth = [2;3;4;6;8;12];
xp = load('nodesRegular.dat');
if key==1
  starting = 1;
else
  starting = sum(lenth(1:key-1))+1; 
end

x = xp(starting:starting+lenth(key)-1,1);
w = xp(starting:starting+lenth(key)-1,2);

end % regularWeights


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function G = stokesRegularizedSLmatrix(o,vesicle,eps)
% G = stokesRegularizedSLtmatrix(vesicle,eps) generates the regularized
% single-layer potential for Stokes vesicle is a data structure defined
% as in the capsules class Gtrap is (2N,2N,nv) array where N is the
% number of points per curve and nv is the number of curves in X 


oc = curve;
[x,y] = oc.getXY(vesicle.X);
% seperate x and y coordinates of vesicle

G = zeros(2*o.N,2*o.N,vesicle.nv);
% initalize single-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:o.N % Loop over targets
    rho2 = (x(j,k)-x(:,k)).^2 + (y(j,k)-y(:,k)).^2;
    rho2 = rho2 + eps^2;
    % regularize the distance
    br = 1./rho2;
    logpart = -1/2*log(rho2);

    G(j,1:o.N,k) = logpart + (x(j,k) - x(:,k)).^2 .* br;
    G(j,o.N+1:2*o.N,k) = (x(j,k)-x(:,k)).*(y(j,k)-y(:,k)) .* br;
    G(j+o.N,1:o.N,k) = G(j,o.N+1:2*o.N,k);
    G(j+o.N,1+o.N:2*o.N,k) = logpart + (y(j,k) - y(:,k)).^2 .* br;
  end % j
  sak = repmat(vesicle.sa(:,k)',2*o.N,2);
  % 2 copies of Jacobian
  
  G(:,:,k) = 2*pi/o.N*G(:,:,k).*sak;
  % multiply both components of G by the Jacobian
  % don't need to divide by 4*pi as it is stored in o.qw
end % k

G = G/4/pi;

end % stokesRegularizedSLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function Gtrap = stokesSLtrapmatrix(o,vesicle)
% Gtrap = stokesSLtrapmatrix(vesicle) generates double-layer potential
% for Stokes vesicle is a data structure defined as in the capsules
% class Gtrap is (2N,2N,nv) array where N is the number of points per
% curve and nv is the number of curves in X 

Gtrap = [];
%Nup = vesicle.N*ceil(sqrt(vesicle.N));
%% upsample to N^(3/2).  
%if Nup < 1024
%  X = vesicle.X; % Vesicle positions
%  Xup = [interpft(X(1:vesicle.N,:),Nup);...
%         interpft(X(vesicle.N+1:2*vesicle.N,:),Nup)];
%
%  oc = curve;
%  [x,y] = oc.getXY(Xup);
%
%  Gtrap = zeros(2*Nup,2*Nup,vesicle.nv);
%  for k=1:vesicle.nv  % Loop over curves
%    for j=1:Nup % Loop over targets
%      rho2 = (x(j,k)-x(:,k)).^2 + (y(j,k)-y(:,k)).^2;
%      rho2(j) = 1;
%      br = 1./rho2;
%      % Set diagonal term to one to avoid dividing by zero
%      logpart = -1/2*log(rho2);
%
%      Gtrap(j,1:Nup,k) = logpart + (x(j,k) - x(:,k)).^2 .* br;
%      Gtrap(j,Nup+1:2*Nup,k) = (x(j,k)-x(:,k)).*(y(j,k)-y(:,k)) .* br;
%      Gtrap(j+Nup,1:Nup,k) = Gtrap(j,Nup+1:2*Nup,k);
%      Gtrap(j+Nup,1+Nup:2*Nup,k) = logpart + (y(j,k) - y(:,k)).^2 .* br;
%    end % j
%  end % k
%
%  Gtrap = Gtrap/4/pi;
%  % constant that multiplies for viscosity contrast.  For solid walls, 
%  % simply set viscCont to 0 to obtain the right scaling
%else
%  Gtrap = [];
%end


end % stokesSLtrapmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function G = stokesSLmatrixOld(o,vesicle)
% G = stokesSLmatrixOld(vesicle) generates the single-layer potential
% for Stokes vesicle is a data structure defined as in the curve class
% G is (2N,2N,nv) array where N is the number of points per curve and
% nv is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% seperate x and y coordinates of vesicle

G = zeros(2*o.N,2*o.N,vesicle.nv);
% initalize single-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:o.N % Loop over targets
    ind = 1 + mod(j-1 + (0:o.N-1),o.N);
    xin = o.qp*x(ind,k);
    yin = o.qp*y(ind,k);

    rho = (xin-x(j,k)).^2 + (yin-y(j,k)).^2;
    br = 1./rho;

    logpart = -1/2 * o.qw .* log(rho);
    
    G(j,ind,k) = (logpart+...
        (o.qw.*( ( x(j,k) - xin).^2 .* br)))'*o.qp;
    G(j, o.N+ind, k) = ...
        (o.qw.*(x(j,k)-xin).*(y(j,k)-yin) .* br)'*o.qp;
    G(j+o.N, 1:o.N ,k) = G(j,o.N+1:2*o.N,k);
    G(j+o.N,o.N+ind,k) = (logpart+...
        (o.qw.*( ( y(j,k) - yin).^2 .* br)))'*o.qp;
  end % j
  sak = repmat(vesicle.sa(:,k)',2*o.N,2);
  % 2 copies of Jacobian
  
  G(:,:,k) = G(:,:,k).*sak;
  % multiply both components of G by the Jacobian
  % don't need to divide by 4*pi as it is stored in o.qw
end % k

end % stokesSLmatrixOld

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function D = stokesDLmatrixOld(o,vesicle)
% D = stokesDLmatrix(vesicle), generate double-layer potential for 
% Stokes vesicle is a data structure defined as in the capsules class
% D is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions

normal = [vesicle.xt(vesicle.N+1:2*vesicle.N,:);...
         -vesicle.xt(1:vesicle.N,:)]; % Normal vector

D = zeros(2*vesicle.N,2*vesicle.N,vesicle.nv);
for k=1:vesicle.nv  % Loop over curves
  if (vesicle.viscCont(k) ~= 1)
    for j=1:vesicle.N % Loop over targets
      rho2 = (x(j,k)-x(:,k)).^2 + (y(j,k)-y(:,k)).^2;
      rho2(j) = 1;
      % Set diagonal term to one to avoid dividing by zero

      coeff = ((x(j,k) - x(:,k)).*normal(1:vesicle.N,k) + ...
          (y(j,k) - y(:,k)).*normal(vesicle.N+1:2*vesicle.N,k)).*...
          vesicle.sa(:,k)./rho2.^2/pi;
      % part of kernel of the double-layer potential

      D(j,:,k) = 2*pi/vesicle.N*[coeff.*(x(j,k) - x(:,k)).^2; ...
        coeff.*(x(j,k) - x(:,k)).*(y(j,k) - y(:,k))]';
      D(j+vesicle.N,:,k) = 2*pi/vesicle.N*[...
        coeff.*(y(j,k) - y(:,k)).*(x(j,k)-x(:,k)); ...
        coeff.*(y(j,k) - y(:,k)).^2]';
      % Build double-layer potential matrix D

      rc = [j j+vesicle.N];
      D(rc,rc,k) = -2*pi/vesicle.N*vesicle.sa(j,k)*...
        vesicle.cur(j,k)/2/pi*...
        [vesicle.xt(j,k)^2 ...
            vesicle.xt(j,k)*vesicle.xt(j+vesicle.N,k);...
        vesicle.xt(j+vesicle.N,k)*vesicle.xt(j,k) ...
            vesicle.xt(j+vesicle.N,k)^2];
      % Diagonal term requires the limiting value.  Otherwise, above
      % formula would divide by zero
    end % j
    D(:,:,k) = D(:,:,k)*(1-vesicle.viscCont(k));
    % constant that multiplies for viscosity contrast.  For solid walls,
    % simply set viscCont to 0 to obtain the right scaling
  end
end % k

end % stokesDLmatrixOld

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function N0 = stokesN0matrixOld(o,vesicle)
% N0 = stokesN0matrixOld(vesicle) generates the the integral operator with
% kernel normal(x) \otimes normal(y) which removes the rank one defficiency of
% the double-layer potential.  Need this operator for solid walls

oc = curve;
[x,y] = oc.getXY(vesicle.X); % Vesicle positions

normal = [vesicle.xt(vesicle.N+1:2*vesicle.N,:);...
         -vesicle.xt(1:vesicle.N,:)]; % Normal vector

N0 = zeros(2*vesicle.N,2*vesicle.N,vesicle.nv);
for k = 1:1 % Loop over curves
  % Only want to form N0 for the outer boundary.
  % Rotlets and Stokeslets take care of null space for inner 
  % boundaries
  for j = 1:2*vesicle.N % Loop over targets
    N0(j,:,k) = normal(j,k)*normal(:,k).*...
        [vesicle.sa(:,k);vesicle.sa(:,k)]*2*pi/vesicle.N;
    % compute the columns of the modification to the double-layer
    % potential.
  end
end
% Use N0 if solving (-1/2 + DLP)\eta = f where f has no flux through
% the boundary.  By solving (-1/2 + DLP + N0)\eta = f, we guarantee
% that \eta also has no flux through the boundary.  This is not
% required, but it means we're enforcing one addition condition on eta
% which removes the rank one kernel.  DLP is the double-layer potential
% for stokes equation

end % stokesN0matrixOld



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function G = laplaceRegularizedSLmatrix(o,vesicle,eps,order)
% G = laplaceRegularizedSLmatrix(vesicle,eps) generates the regularized
% single layer potential for laplace.  vesicle is a data structure
% defined as in the capsules class G is (N,N,nv) array where N is the
% number of points per curve and nv is the number of curves in X

oc = curve;
[x,y] = oc.getXY(vesicle.X); % Vesicle positions

G = zeros(o.N,o.N,vesicle.nv);
% initalize single-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:o.N % Loop over targets
    rho2 = (x(j,k)-x(:,k)).^2 + (y(j,k)-y(:,k)).^2;
    logpart = -1/2*log(rho2 + eps^2);
    % regularized kernel

    higherOrderTerms = zeros(o.N,1);
    for i = 1:order
      higherOrderTerms = higherOrderTerms + 0.5/i*eps^(2*i)*(rho2+eps^2).^(-i);
    end

    G(j,1:o.N,k) = logpart + higherOrderTerms;

  end % j
  sak = repmat(vesicle.sa(:,k)',o.N,1);
  
  G(:,:,k) = 1/4/pi*G(:,:,k).*sak*2*pi/o.N;
  % multiply G by the Jacobian
end % k


end % laplaceRegularizedSLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function G = laplaceD4SLmatrix(o,vesicle)
% G = laplaceD4SLmatrix(vesicle), generate the fourth derivative of the
% single-layer potential for laplace.  This can be used when the
% derivatives are moved off of the bending term and onto the layer
% potential.  vesicle is a data structure defined as in the capsules
% class G is (N,N,nv) array where N is the number of points per curve
% and nv is the number of curves in X

oc = curve;
[x,y] = oc.getXY(vesicle.X); % Vesicle positions
[tanx,tany] = oc.getXY(vesicle.xt); % Vesicle tangent
norx = tany;
nory = -tanx;
cur = vesicle.cur;
isa = 1./vesicle.sa;
Dcur = curve.arcDeriv(cur,1,isa,vesicle.IK);
D2cur = curve.arcDeriv(Dcur,1,isa,vesicle.IK);

G = zeros(o.N,o.N,vesicle.nv);
% initalize single-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:o.N % Loop over targets
    ind = 1 + mod(j-1 + (0:o.N-1),o.N);
    xin = o.qp*x(ind,k);
    yin = o.qp*y(ind,k);
    curin = o.qp*cur(ind,k);
    Dcurin = o.qp*Dcur(ind,k);
    D2curin = o.qp*D2cur(ind,k);

    rx = xin - x(j,k); ry = yin - y(j,k);
    rdotn = rx.*(o.qp*norx(ind,k)) + ry.*(o.qp*nory(ind,k));
    rdott = rx.*(o.qp*tanx(ind,k)) + ry.*(o.qp*tany(ind,k));
    rho2 = rx.^2 + ry.^2;

    if 1 % one derivative
      termOrder1 = rdott./rho2;
      termOrder2 = zeros(size(termOrder1));
      termOrder3 = zeros(size(termOrder1));
      termOrder4 = zeros(size(termOrder1));
    end
    if 0
      termOrder1 = (1+curin.*rdotn)./rho2;
      termOrder2 = -2*(rdott).^2./rho2.^2;
      termOrder3 = zeros(size(termOrder1));
      termOrder4 = zeros(size(termOrder1));
    end



    if 0 % four derivatives
      termOrder1 = (D2curin.*rdotn - ...
          3*curin.*Dcurin.*rdott - curin.^2 - ...
          curin.^3.*rdotn)./rho2;
      termOrder2 = (-8*Dcurin.*rdott.*rdotn + ...
          8*curin.^2.*rdott.^2 - 6*curin.^2.*rdotn.^2 -...
          12*curin.*rdotn - 6)./rho2.^2;
      termOrder3 = (48*rdott.^2 + ...
          48*curin.*rdott.^2.*rdotn)./rho2.^3;
      termOrder4 = -48*rdott.^4./rho2.^4;
    end

    G(j,ind,k) = -1*(o.qw.*(termOrder1 + termOrder2 + ...
        termOrder3 + termOrder4))'*o.qp;
    
  end % j
  sak = repmat(vesicle.sa(:,k)',o.N,1);
  
  G(:,:,k) = G(:,:,k).*sak;
  % multiply G by the Jacobian
end % k

end % laplaceD4SLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function G = FourthDerivLogMatrix(o,vesicle)
% G = FourthDerivLogMatrix(vesicle) generates the matrix for the kernel
% that is the fourth derivative of log multiplied by rho^4 which is the
% strength of the singularity
oc = curve;
[x,y] = oc.getXY(vesicle.X); % Vesicle positions
[tanx,tany] = oc.getXY(vesicle.xt); % Vesicle tangent
norx = tany;
nory = -tanx;
cur = vesicle.cur;
isa = 1./vesicle.sa;
Dcur = curve.arcDeriv(cur,1,isa,vesicle.IK);
D2cur = curve.arcDeriv(Dcur,1,isa,vesicle.IK);

eps = 0.1;

G = zeros(o.N,o.N,vesicle.nv);
% initalize layer potential to zero

for k=1:vesicle.nv  % Loop over curves
  for j=1:vesicle.N % Loop over targets
    rdotn = (x(:,k) - x(j,k)).*norx(:,k) + ...
        (y(:,k) - y(j,k)).*nory(:,k);
    rdott = (x(:,k) - x(j,k)).*tanx(:,k) + ...
        (y(:,k) - y(j,k)).*tany(:,k);
    rho2 = ((x(:,k)-x(j,k)).^2 + (y(:,k)-y(j,k)).^2);

    term1 = (1*D2cur.*rdotn - 1*3.*cur.*Dcur.*rdott - ...
        1*cur.^2 - 1*cur.^3.*rdotn).*rho2;
    term2 = (-8*Dcur.*rdott.*rdotn + 8*cur.^2.*rdott.^2 - ...
        6*cur.^2.*rdotn.^2 - 12*cur.*rdotn - 6);
    term3 = (48*rdott.^2 + 48*cur.*rdott.^2.*rdotn)./rho2;
    term3(j) = 48;
    term4 = -48*rdott.^4./rho2.^2;
    term4(j) = -48;

    coeff = term1 + term2 + term3 + term4;
    G(j,:,k) = coeff;
    % kernel of fourth-derivative of log

  end % j
end % k


end % FourthDerivLogMatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function G = FourthDerivLogRegularizedMatrix(o,vesicle,eps)
% G = FourthDerivLogRegularizedMatrix(vesicle) generates the matrix for
% the kernel that is the fourth derivative of 1/2*log(rho^2 + eps^2)
% multiplied by (rho^2+eps^2)^2 which is the strength of the
% singularity as eps goes to zero

oc = curve;
[x,y] = oc.getXY(vesicle.X); % Vesicle positions
[tanx,tany] = oc.getXY(vesicle.xt); % Vesicle tangent
norx = tany;
nory = -tanx;
cur = vesicle.cur;
isa = 1./vesicle.sa;
Dcur = curve.arcDeriv(cur,1,isa,vesicle.IK);
D2cur = curve.arcDeriv(Dcur,1,isa,vesicle.IK);

G = zeros(o.N,o.N,vesicle.nv);
% initalize layer potential to zero

for k=1:vesicle.nv  % Loop over curves
  for j=1:vesicle.N % Loop over targets
    rdotn = (x(:,k) - x(j,k)).*norx(:,k) + ...
        (y(:,k) - y(j,k)).*nory(:,k);
    rdott = (x(:,k) - x(j,k)).*tanx(:,k) + ...
        (y(:,k) - y(j,k)).*tany(:,k);
    rho2 = ((x(:,k)-x(j,k)).^2 + (y(:,k)-y(j,k)).^2);
    rho2 = rho2 + eps^2;

    term1 = (1*D2cur.*rdotn - 1*3.*cur.*Dcur.*rdott - ...
        1*cur.^2 - 1*cur.^3.*rdotn).*rho2;
    term2 = (-8*Dcur.*rdott.*rdotn + 8*cur.^2.*rdott.^2 - ...
        6*cur.^2.*rdotn.^2 - 12*cur.*rdotn - 6);
    term3 = (48*rdott.^2 + 48*cur.*rdott.^2.*rdotn)./rho2;
    term4 = -48*rdott.^4./rho2.^2;

    coeff = term1 + term2 + term3 + term4;
%    G(j,:,k) = coeff;
    G(j,:,k) = coeff./rho2.^2;
    % kernel of fourth-derivative of log

  end % j
end % k


end % FourthDerivLogRegularizedMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function G = distanceRatioMatrix(o,vesicle)
% function G = distanceRatioMatrix(vesicle) computes the matrix for the
% layer potential with the kernel rho_{c}^{4}/rho^{4} where rho_{c} is
% the distance measured on a circle of the same length as vesicle.

theta = (0:vesicle.N-1)'*2*pi/vesicle.N;
rad = vesicle.length/2/pi;
xc = rad*cos(theta); yc = rad*sin(theta);

oc = curve;
[x,y] = oc.getXY(vesicle.X);

G = zeros(o.N,o.N,vesicle.nv);
% initialize layer potential to zero

for k=1:vesicle.nv  % Loop over curves
  for j=1:vesicle.N % Loop over targets
    rho2 = (x(:,k)-x(j,k)).^2 + (y(:,k)-y(j,k)).^2;
    rho2c = (xc(:,k)-xc(j,k)).^2 + (yc(:,k)-yc(j,k)).^2;
    rho2(j) = 1;
    % Set diagonal term to one to avoid dividing by zero

    coeff = rho2c.^2./rho2.^2;
    % kernel 
    coeff(j) = rad.^4./vesicle.sa(j).^4;
    % diagonal term
    G(j,:,k) = coeff;
    % FourthOrderSingularityMatrix has the arclength and jacobian terms in it

  end % j
end % k

end % distanceRatioMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function G = rhocMatrix(o,vesicle)
% function G = rhocoMatrix(vesicle) computes the matrix for the layer
% potential with the kernel 1/rho_{c}^{4} where rho_{c} is the distance
% measured on a circle of the same length as vesicle.

theta = (0:vesicle.N-1)'*2*pi/vesicle.N;
rad = vesicle.length/2/pi;

G = zeros(o.N,o.N,vesicle.nv);
% initialize layer potential to zero

for k = 1:vesicle.nv
  for j = 1:vesicle.N
    rhoc2 = 4*rad^4 * (1-cos(theta - theta(j)));
    G(j,:,k) = 1./rhoc2.^2;
    G(j,j,k) = 0;
  end % j
end % k


end % rhocMatrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function LP = nearSingInt_New(o,vesicleSou,f,selfMat,selfMatTrap,...
%    NearStruct,kernel,kernelDirect,vesicleTar,tEqualS,trash)
% LP = nearSingIntNew(vesicle,f,selfMat,selfMatTrap,...
% NearStruct,kernel,kernelDirect,vesicleTar,tEqualS,trash) computes a
% layer potential due to f at all points in vesicleTar.X.  If
% tEqualS==true, then the vesicleTar == vesicleSou and the self-vesicle
% interaction is skipped.  selfMat is the diagonal of the potential
% needed to compute the layer potential of each vesicle indepenedent of
% all others.  kernel and kernelDirect are two (possibly the same)
% routines that compute the layer potential.  kernelDirect always uses
% the direct method whereas kernel may use an FMM-accelerated method.
% NearStruct is a structure containing the variables
% zone,dist,nearest,icp,argnear which are required by near-singular
% integration (they keep everything sorted and precomputed) Everything
% is in the 2*N x nv format Can pass a final argument if desired so
% that plots of the near-singular integration algorithm are displayed

dist = NearStruct.dist;
zone = NearStruct.zone;
nearest = NearStruct.nearest;
icp = NearStruct.icp;
argnear = NearStruct.argnear;

Xsou = vesicleSou.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
nvSou = size(Xsou,2); % number of source 'vesicles'
Xtar = vesicleTar.X; % target positions
Ntar = size(Xtar,1)/2; % number of target points
nvTar = size(Xtar,2); % number of target 'vesicles'

h = vesicleSou.length/Nsou; % arclength term

%Nup = Nsou*2^ceil(1/2*log2(Nsou));
Nup = Nsou*ceil(sqrt(Nsou));
% upsample to N^(3/2).  

vself = selfMat(f);
%vself = zeros(2*Nsou,nvSou);
%for k=1:nvSou
%  vself(:,k) = selfMat(:,:,k)*f(:,k);
%end
% Compute velocity due to each vesicle independent of others.
% This is needed when using near-singular integration since
% we require a value for the layer-potential on the vesicle of 
% sources 

Xup = [interpft(Xsou(1:Nsou,:),Nup);...
       interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];
fup = [interpft(f(1:Nsou,:),Nup);...
       interpft(f(Nsou+1:2*Nsou,:),Nup)];
% upsample positions, traction jump

vesicleUp = capsules(Xup,[],[],...
    vesicleSou.kappa,vesicleSou.viscCont,false);
% Build an object with the upsampled vesicle

interpOrder = size(o.interpMat,1);
% lagrange interpolation order
p = ceil((interpOrder+1)/2);
% want half of the lagrange interpolation points to the left of the
% closest point and the other half to the right
vel = zeros(2*Ntar,nvTar,nvSou);
% allocate space for storing velocity at intermediate points needed
% by near-singular integration

if tEqualS % sources == targets
  if nvSou > 1
    if strfind(char(kernel),'fmm');
%      if numel(fup) >= 256
        farField = kernel(vesicleUp,fup,selfMatTrap);
%      else
%        farField = kernelDirect(vesicleUp,fup,[]);
%      end
      farField = farField(1:Nup/Ntar:end,:);
      % evaluate layer potential at all targets except ignore the
      % diagonal term
    else
      for k = 1:nvSou
        K = [(1:k-1) (k+1:nvSou)];
        [~,farField(:,k)] = kernel(vesicleUp,fup,[],Xtar(:,k),K);
      end
%      % This is a huge savings if we are using a direct method rather
%      % than the fmm to evaluate the layer potential.  The speedup is
%      % more than N^{1/2}, where N is the resolution of the vesicles
%      that % we are computing on
    end
  else
    farField = zeros(2*Ntar,nvTar);
  end

else % sources ~= targets
%  if (numel(fup) + numel(Xtar)/2) >= 256
    [~,farField] = kernel(vesicleUp,fup,[],Xtar,1:nvSou);
%  else
%    [~,farField] = kernelDirect(vesicleUp,fup,[],Xtar,1:nvSou);
%  end
  % evaluate layer potential due to all 'vesicles' at all points
  % in Xtar;
end
% Use upsampled trapezoid rule to compute layer potential

nearField = zeros(2*Ntar,nvTar);
% Initialize potential at near points to zero

beta = 1.1;
% small buffer to make sure Lagrange interpolation points are
% not in the near zone
for k1 = 1:nvSou
  if tEqualS % sources == targets
    K = [(1:k1-1) (k1+1:nvTar)];
    % skip diagonal vesicle
  else % sources ~= targets
    K = (1:nvTar);
    % consider all vesicles
  end
  for k2 = K
    J = find(zone{k1}(:,k2) == 1);
    % set of points on vesicle k2 close to vesicle k1
    if (numel(J) ~= 0)
      indcp = icp{k1}(J,k2);
      % closest point on vesicle k1 to each point on vesicle k2 
      % that is close to vesicle k1
      for j = 1:numel(J)
        pn = mod((indcp(j)-p+1:indcp(j)-p+interpOrder)' - 1,Nsou) + 1;
        % index of points to the left and right of the closest point
        v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
          o.interpMat*vself(pn,k1));
        vel(J(j),k2,k1) = v(end);
        % x-component of the velocity at the closest point
        v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
          o.interpMat*vself(pn+Nsou,k1));
        vel(J(j)+Ntar,k2,k1) = v(end);
        % y-component of the velocity at the closest point
      end
%     compute values of velocity at required intermediate points
%     using local interpolant

%      if numel(J) >= 256
      [~,potTar] = kernel(vesicleUp,fup,[],...
         [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
%      else
%        [~,potTar] = kernelDirect(vesicleUp,fup,[],...
%           [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
%      end
      % Need to subtract off contribution due to vesicle k1 since
      % its layer potential will be evaulted using Lagrange
      % interpolant of nearby points
      nearField(J,k2) =  - potTar(1:numel(J));
      nearField(J+Ntar,k2) =  - potTar(numel(J)+1:end);

      XLag = zeros(2*numel(J),interpOrder - 1);
      nxOld = 0; nyOld = 0; indOld = 0;
      % need to compare the (approximate) normal direction with the
      % previous one so that we aren't ever putting two copies of the
      % same target point through the FMM
      Jind = zeros(size(J));
      for i = 1:numel(J)
        nx = (Xtar(J(i),k2) - nearest{k1}(J(i),k2))/...
            dist{k1}(J(i),k2);
        ny = (Xtar(J(i)+Ntar,k2) - nearest{k1}(J(i)+Ntar,k2))/...
            dist{k1}(J(i),k2);
        if ((nx - nxOld).^2 + (ny - nyOld).^2 > 1e-10)
          XLag(i,:) = nearest{k1}(J(i),k2) + ...
              beta*h*nx*(1:interpOrder-1);
          XLag(i+numel(J),:) = nearest{k1}(J(i)+Ntar,k2) + ...
              beta*h*ny*(1:interpOrder-1);
          nxOld = nx;
          nyOld = ny;
          indOld = indOld+1;
          Jind(i) = indOld;
        else
          Jind(i) = indOld;
        end
        % Lagrange interpolation points coming off of vesicle k1 All
        % points are behind Xtar(J(i),k2) and are sufficiently far from
        % vesicle k1 so that the Nup-trapezoid rule gives sufficient
        % accuracy
      end
      XLag = XLag(sum(XLag.^2,2)~=0,:);
      % remove the columns that all contain zeros since they correspond
      % to Lagrange interpolation points that will computed for another
      % target points

%      if numel(XLag)/2 >= 256
      [~,lagrangePtsCompressed] = kernel(vesicleUp,fup,[],XLag,k1);
      % evaluate velocity at the lagrange interpolation points
      lagrangePts = zeros(2*numel(J),interpOrder - 1);
      lagrangePts(1:numel(J),:) = lagrangePtsCompressed(Jind,:);
      lagrangePts(numel(J)+1:2*numel(J),:) = ...
          lagrangePtsCompressed(Jind+max(Jind),:);

      for i = 1:numel(J)
        Px = o.interpMat*[vel(J(i),k2,k1) ...
            lagrangePts(i,:)]';
        Py = o.interpMat*[vel(J(i)+Ntar,k2,k1) ...
            lagrangePts(i+numel(J),:)]';
        % Build polynomial interpolant along the one-dimensional
        % points coming out of the vesicle
        dscaled = full(dist{k1}(J(i),k2)/(beta*h*(interpOrder-1)));
        % Point where interpolant needs to be evaluated

        v = filter(1,[1 -dscaled],Px);
        nearField(J(i),k2) = nearField(J(i),k2) + ...
            v(end);
        v = filter(1,[1 -dscaled],Py);
        nearField(J(i)+Ntar,k2) = nearField(J(i)+Ntar,k2) + ...
            v(end);
        % Assign higher-order results coming from Lagrange 
        % integration to velocity at near point.  Filter is faster
        % than polyval

        if (nargin == 11)
          figure(2); clf; hold on;
          plot(Xsou(1:Nsou,:),Xsou(Nsou+1:end,:),'r.')
          plot(Xtar(1:Ntar,:),Xtar(Ntar+1:end,:),'k.')
          plot(Xtar(J,k2),Xtar(Ntar+J,k2),'b.')
          plot(XLag(1:numel(J),:),XLag(numel(J)+1:end,:),'kx')
          plot(XLag(i,:),XLag(numel(J)+i,:),'gx')
          axis equal

          figure(1); clf; hold on
          plot((0:interpOrder-1)*beta*h,...
              real([vel(J(i),k2,k1) lagrangePts(i,:)]),'g-o')
          plot((0:interpOrder-1)*beta*h,...
              real([vel(J(i)+Ntar,k2,k1) lagrangePts(i+numel(J),:)]),'r--o')
          pause
        end
        % DEBUG: PASS IN A DUMMY VARIABLE INTO THIS ROUTINE AND THEN
        % YOU CAN SEE THE INTERPOLATION POINTS AND CHECK THE SMOOTHNESS
        % OF THE INTERPOLANT

      end % i
    end % numel(J) ~= 0
    % Evaluate layer potential at Lagrange interpolation
    % points if there are any
  end % k2
end % k1

%if tEqualS % vesicleSou == vesicleTar
%  LP = farField(1:Nup/Ntar:end,:) + nearField;
%else % vesicleSou ~= vesicleTar
%  LP = farField + nearField;
%end
LP = farField + nearField;
% Add kernel due to far points and near points.  Far points were
% upsampled if source==vesicle so need to truncate here.  We are 
% only using Ntar target points.  Note that it is only the sources 
% that were upsampled


end % nearSingInt_New


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [laplaceDLP,laplaceDLPtar] = ...
%    exactLaplaceDLfmmOld(o,vesicle,f,DLPtrap,Xtar,K)
% [laplaceDLP,laplaceDLPtar] = exactLaplaceDLfmm(vesicle,f) uses 
% the FMM to compute the double-layer laplace potential due to all 
% vesicles except itself vesicle is a class of object capsules and 
% f is the density function NOT scaled by arclength term.  Xtar is 
% a set of points where the single-layer potential due to all 
% vesicles in index set K

oc = curve;
[x,y] = oc.getXY(vesicle.X); % seperate x and y coordinates
nx = vesicle.xt(vesicle.N+1:2*vesicle.N,:);
ny = -vesicle.xt(1:vesicle.N,:);

[f,~] = oc.getXY(f);
den = f.*vesicle.sa*2*pi/vesicle.N;

if nargin == 6
  laplaceDLP = [];
else
  [~,rfield,~] = fmm_laplace(nx(:).*den(:),x(:),y(:),1);
  [~,~,cfield] = fmm_laplace(ny(:).*den(:),x(:),y(:),1);
  % double-layer potential due to all vesicles
  potTot = rfield + cfield;

  laplaceDLP = zeros(vesicle.N,vesicle.nv); % initialize
  for k = 1:vesicle.nv
    is = (k-1)*vesicle.N+1;
    ie = k*vesicle.N;
    laplaceDLP(1:vesicle.N,k) = potTot(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format

  for k = 1:vesicle.nv
    [~,rfield,~] = fmm_laplace(nx(:,k).*den(:,k),x(:,k),y(:,k),1);
    [~,~,cfield] = fmm_laplace(ny(:,k).*den(:,k),x(:,k),y(:,k),1);
    laplaceDLP(:,k) = laplaceDLP(:,k) - rfield - cfield;
  end
  % Subtract potential due to each vesicle on its own.  Nothing
  % to change here for potential at Xtar
end
laplaceDLP = -[laplaceDLP;zeros(vesicle.N,vesicle.nv)]/(2*pi);
% pad with zeros so that it is compatible with Stokes near-singular
% integration

if nargin == 4
  laplaceDLPtar = [];
else
  x = x(:,K);
  y = y(:,K);
  nx = nx(:,K);
  ny = ny(:,K);
  den = den(:,K);
  % only care about points on vesicles indexed by K
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  den = [den(:);zeros(Ntar*ncol,1)];
  nx = [nx(:);zeros(Ntar*ncol,1)];
  ny = [ny(:);zeros(Ntar*ncol,1)];
  % pad density function and normal with zeros so that Xtar doesn't
  % affect the single-layer potential
  [~,rfield,~] = fmm_laplace(nx.*den,x(:),y(:),1);
  [~,~,cfield] = fmm_laplace(ny.*den,x(:),y(:),1);
  % evalaute single-layer potential using fmm

  laplaceDLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = vesicle.N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    laplaceDLPtar(1:Ntar,k) = rfield(is:ie)+cfield(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
  laplaceDLPtar = -laplaceDLPtar/(2*pi);
end


end % exactLaplaceDLfmmOld



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD ROUTINES FROM capsules.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [D4,Ten,Div] = computeDerivsOld(o)
% [D4,Ten,Div] = compute4Deriv computes the matricies that 
% takes a periodic function and maps it to the fourth 
% derivative, tension, and surface divergence all with 
% respect to arclength.
global derivs

persistent FF
% persistent variable for matrix that takes physical space to
% frequency space

derivs = derivs + 1;
if isempty(FF);
  FF = fft1.fourierInt(o.N);
end
% Build Fourier interpolation matrix if it is not constructed
% yet.  It can be used throught the whole simulation as it only
% depends on the spatial discretization N


theta = (0:o.N-1)'*2*pi/o.N;
modes = (-o.N/2:o.N/2-1)'; % Fourier modes

D4 = zeros(2*o.N,2*o.N,o.nv); % Bending
Ten = zeros(2*o.N,o.N,o.nv); % Tension
Div = zeros(o.N,2*o.N,o.nv); % Surface divergence

for k=1:o.nv % Loop over vesicles
  x1 = [o.xt(1:o.N,k) o.xt(o.N+1:2*o.N,k)];
  % tangent vector
  o.IK(o.N/2+1,k) = -o.N/2*1i;
  % this makes sure that the preconditioner is invertible so
  % that we can precompute its inverse with inv which is much
  % faster than pinv
  IK = o.IK(:,k);
%  tic
  for i=1:o.N % Loop over modes
%    rDx1 = curve.arcDeriv(cos(modes(i)*theta),4,1./o.sa(:,k),IK);
%    iDx1 = curve.arcDeriv(sin(modes(i)*theta),4,1./o.sa(:,k),IK);
    rDx1 = curve.arcDeriv(cos(modes(i)*theta),4,o.isa(:,k),IK);
    iDx1 = curve.arcDeriv(sin(modes(i)*theta),4,o.isa(:,k),IK);
    D4(1:o.N,i,k) = rDx1 + 1i*iDx1;
    % Fourth arclength derivative
  end % i

%  D4(o.N+1:2*o.N,o.N+1:2*o.N,k) = D4(1:o.N,1:o.N,k);
%  D4(:,:,k) = D4(:,:,k) * [FF zeros(o.N); zeros(o.N) FF];
%  % This is considerabely slower.  I'm not sure why I ever did this

  D4(1:o.N,1:o.N,k) = D4(1:o.N,1:o.N,k)*FF;
  D4(o.N+1:2*o.N,o.N+1:2*o.N,k) = D4(1:o.N,1:o.N,k);
  % Move to physical space
%  toc

  % two components because of tension times both x and y coordinates
%  tic
  for i=1:o.N % Loop over modes
    Drsig = curve.arcDeriv(...
        [cos(modes(i)*theta).*x1(:,1) cos(modes(i)*theta).*x1(:,2)],...
        1,[o.isa(:,k) o.isa(:,k)],[IK IK]);
    Disig = curve.arcDeriv(...
        [sin(modes(i)*theta).*x1(:,1) sin(modes(i)*theta).*x1(:,2)],...
        1,[o.isa(:,k) o.isa(:,k)],[IK IK]);
    Ten(1:2*o.N,i,k) = Drsig(:) + 1i*Disig(:);
  end % i

  Ten(:,:,k) = Ten(:,:,k) * FF;
  % Move to physical space
%  toc

%  tic
  for i=1:o.N % Loop over modes
    rDx1 = curve.arcDeriv(cos(modes(i)*theta),1,o.isa(:,k),IK);
    iDx1 = curve.arcDeriv(sin(modes(i)*theta),1,o.isa(:,k),IK);
    Div(1:o.N,i,k) = (rDx1 + 1i*iDx1) .* o.xt(1:o.N,k);
    Div(1:o.N,i+o.N,k) = (rDx1 + 1i*iDx1) .* o.xt(o.N+1:2*o.N,k);
  end % i

  Div(:,1:o.N,k) = Div(:,1:o.N,k)*FF;
  Div(:,o.N+1:end,k) = Div(:,o.N+1:end,k)*FF;
  % Move to physical spacce
%  toc
end % k

D4 = real(D4); 
Ten = real(Ten);
Div = real(Div);
% Imaginary part should be 0 since we are preforming a 
% real operation
% When we build the preconditioner, it doesn't matter that 
% D4, Ten, and Div zero high frequencies because of the identity
% term that is in the linear system.

end % computeDerivsOld

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [sigma,eta,RS,velVes,iter] = computeSigAndEtaOld(vesicle,tt,walls)
% [sigma,eta,RS,velVes,iter] = computeSigAndEta(vesicle,tt,walls)
% computes the tension, density function (if flow is confined), and
% velocity of the vesicle for the vesicle configuration given by
% vesicle.X and solid wall configuration given by walls.X

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles
if tt.confined
  Nbd = walls.N; % points per solid wall
  nvbd = walls.nv; % number of solid walls
else
  Nbd = 0;
  nvbd = 0;
end

op = tt.op;
[Ben,Ten,Div] = vesicle.computeDerivs;
% compute self bending, tension, and divegence terms
tt.Galpert = op.stokesSLmatrix(vesicle);
%tt.Gbarnett = op.laplaceSLcomplexMatrix(vesicle);
tt.Gbarnett = [];
% single-layer potential
tt.D = op.stokesDLmatrix(vesicle);
% double-layer potential

if tt.near
  if tt.confined
    [tt.NearV2V,tt.NearV2W] = vesicle.getZone(walls,3);
    [tt.NearW2W,tt.NearW2V] = walls.getZone(vesicle,3);
  else
    tt.NearV2V = vesicle.getZone([],1);
    tt.NearW2V = [];
    tt.NearV2W = [];
    tt.NearW2W = [];
  end
  % near-singular integration structures
else
  tt.NearV2V = [];
  tt.NearW2V = [];
  tt.NearV2W = [];
  tt.NearW2W = [];
  % empty since we are not doing near-singular integration
end

rhs = zeros(N*nv + 2*Nbd*nvbd + 3*(nvbd-1),1);
% initalize the right-hand side for the tension followed by the density
% function followed by the rotlets and stokeslets

if ~tt.fmm
  kernel = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
end
% kernel for single-layer potential.  Only difference is if the FMM is
% used or not

%f = vesicle.bendingTerm(vesicle.X);
f = vesicle.tracJump(vesicle.X,zeros(N,nv));
% traction jump due only to the position
if ~tt.near
  Fslp = kernel(vesicle,f,[]);
  if tt.confined
    [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
  end
else
  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);
  SLPtrap = SLP;
  kernelDirect = kernel;
  Fslp = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
      tt.NearV2V,kernel,kernelDirect,vesicle,true);

  if tt.confined
    FSLPwall = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
      tt.NearV2W,kernel,kernelDirect,walls,false);
    % Evaluate single-layer potential due to all vesicles on
    % the solid walls WITH near-singular integration
  end
end
% single-layer potential due to all vesicles except the diagonal one

if ~tt.confined
  Ffar = tt.farField(vesicle.X);
else
  Ffar = zeros(2*N,nv);
  U = tt.farField(walls.X);
  % no slip boundary condition for velocity on solid walls
end
% If doing unbounded flow, add in the background velocity

velBen = Fslp + op.exactStokesSLdiag(vesicle,tt.Galpert,f) + Ffar;

if any(vesicle.viscCont ~= 1)
  velBen = tt.solveIminusD(velBen,vesicle);
end
% velocity due to bending term.  The position is known and we are
% solving for the other two variables (eta and sigma)

DivVelBen = vesicle.surfaceDiv(velBen);
for k = 1:nv
  istart = (k-1)*N + 1;
  iend = istart + N - 1;
  rhs(istart:iend) = -DivVelBen(:,k);
end

for k = 1:nvbd
  istart = nv*N+1 + (k-1)*2*Nbd;
  iend = istart + 2*Nbd - 1;
  rhs(istart:iend) = U(:,k) - FSLPwall(:,k);
end

if any(vesicle.viscCont ~= 1) & tt.confined
  if ~tt.fmmDLP
    kernel = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLfmm;
  end

  if ~tt.near
    [~,FwallDLP] = kernel(vesicle,velBen,[],walls.X,(1:nv));
  else
    jump = -0.5; % jump condition of double-layer potential
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,tt.D,X);
    DLPtrap = DLP;
    kernelDirect = kernel;
    FwallDLP = op.nearSingInt(vesicle,velBen,DLP,DLPtrap,...
        tt.NearV2W,kernel,kernelDirect,walls,false);
  end
  % contribution from the double-layer potential due to the velocity
  % due to bending with viscosity contrast evaluated along the solid
  % walls
  rhs(nv*N+1:nv*N+2*nvbd*Nbd) = rhs(nv*N+1:nv*N+2*nvbd*Nbd) - ...
      FwallDLP(:);
end

for k = 1:nv
  alpha = 0.5*(1 + vesicle.viscCont(k));
  tt.bdiagTen(:,:,k) = inv(Div(:,:,k)*...
      inv(alpha*eye(2*N) - tt.D(:,:,k))*...
      tt.Galpert(:,:,k)*Ten(:,:,k));
end
% build part of the block-diagonal preconditioner


[sigDen,F,R,I] = gmres(@(X) tt.sigDenMatVec(X,vesicle,walls),rhs,...
    [],tt.gmresTol,min(tt.gmresMaxIter,N*nv+2*Nbd*nvbd + 3*(nvbd-1)),...
    @tt.preconditionerTen);
% solve for tension and density function with block-diagonal 
% preconditioning
iter = I(2);

sigma = zeros(N,nv); % tension
eta = zeros(2*Nbd,nvbd); % solid wall density
RS = zeros(3,nvbd); % rotlets and stokeslets
for k = 1:nv
  sigma(:,k) = sigDen((k-1)*N+1:k*N);
end
for k = 1:nvbd
  eta(:,k) = sigDen(2*(k-1)*Nbd+nv*N+1:2*k*Nbd+nv*N);
end
for k = 2:nvbd
  istart = nv*N + 2*nvbd*Nbd + 3*(k-2) + 1;
  iend = istart + 2;
  RS(:,k) = sigDen(istart:iend);
end
% unstack the tension and density function from the GMRES solver

f = vesicle.tracJump(vesicle.X,sigma);
% traction jump due to the position and tension
if ~tt.near
  Fslp = kernel(vesicle,f,[]);
  if tt.confined
    [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
  end
else
  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);
  SLPtrap = SLP;
  kernelDirect = kernel;
  Fslp = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
      tt.NearV2V,kernel,kernelDirect,vesicle,true);
end
% single-layer potential due to all vesicles except the diagonal one

if tt.confined
  if ~tt.fmmDLP
    kernel = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLfmm;
  end

  if ~tt.near
    [~,Fwall2Ves] = kernel(walls,eta,[],vesicle.X,1:nvbd);
  else
    jump = -1/2;
    DLP = @(X) jump*X + op.exactStokesDLdiag(walls,tt.wallDLP,X);
    kernelDirect = @op.exactStokesDL;
    Fwall2Ves = op.nearSingInt(walls,eta,DLP,DLP,...
        tt.NearW2V,kernel,kernelDirect,vesicle,false);
  end
  
  LetsVes = zeros(2*N,nv);
  for k = 2:nvbd
    LetsVes = LetsVes + tt.RSlets(vesicle.X,walls.center(:,k),...
        RS(1:2,k),RS(3,k));
  end
  % contribution from the rotlets and stokeslets
else
  Fwall2Ves = zeros(2*N,nv);
  LetsVes = zeros(2*N,nv);
end

velVes = Ffar + op.exactStokesSLdiag(vesicle,tt.Galpert,f) + ...
    Fslp + Fwall2Ves + LetsVes;
% velocity of the vesicle if there is no viscosity contrast.
% Otherwise, need to apply inv(alpha*I - D) to this variable.

if any(vesicle.viscCont ~= 1)
  velVes = tt.solveIminusD(velVes,vesicle);
end
% velVes is the velocity along the vesicle

end % computeSigAndEtaOld


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [eta,RS,iter] = computeEta(walls,tt)
% [eta,RS,iter] = computeEta(tt,walls) computes the density function
% for the solid wall configuration given by walls.X

Nbd = walls.N; % points per solid wall
nvbd = walls.nv; % number of solid walls
U = tt.farField(walls.X); % no-slip on solid walls

op = tt.op;

if tt.near
  tt.NearW2W = walls.getZone([],1);
  % near-singular integration structure
else
  tt.NearW2W = [];
  % empty since we are not doing near-singular integration
end

rhs = zeros(2*Nbd*nvbd + 3*(nvbd-1),1);
% initalize the right-hand side for the density function 
% followed by the rotlets and stokeslets

for k = 1:nvbd
  istart = (k-1)*2*Nbd + 1;
  iend = istart + 2*Nbd - 1;
  rhs(istart:iend) = U(:,k);
end

Pre = tt.bdiagWall;
% Build the block-diagonal preconditioner for solving for the tension
% and density function given a configuration


if 1
  [den,F,R,I] = gmres(@(X) tt.denMatVec(X,walls),...
      rhs,[],tt.gmresTol,min(200,2*Nbd*nvbd + 3*(nvbd-1)),...
      @(X) Pre*X);
  iter = I(2);
end
% solve for tension and density function with block-diagonal
% preconditioning


eta = zeros(2*Nbd,nvbd);
RS = zeros(3,nvbd);
for k = 1:nvbd
  eta(:,k) = den(2*(k-1)*Nbd+1:2*k*Nbd);
end
for k = 2:nvbd
  istart = 2*nvbd*Nbd + 3*(k-2) + 1;
  iend = istart + 2;
  RS(:,k) = den(istart:iend);
end
% unstack the density function from the GMRES solver


%kernel = @op.exactStokesDL;
%valWalls = zeros(2*Nbd,nvbd);
%FDLPwall2wall = kernel(walls,eta);
%% compute the velocity on the solid walls due to the
%% solid walls without near-singular integration
%% END OF EVALUATING WALL TO WALL INTERACTIONS
%
%
%% START OF EVALUATING VELOCITY ON WALLS
%for k = 1:nvbd
%  valWalls(:,k) = valWalls(:,k) -1/2*eta(:,k) + ...
%      (tt.wallDLP(:,:,k) + tt.wallN0(:,:,k))*eta(:,k);
%end

end % computeEta


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD ROUTINES FROM curve.m  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [x,y] = diffuserGeom(o,N)
% make geometry of a diffuser so that we can look at sorting vesicles
% (maybe by viscosity contrast)
%
%
%               ((H-h)/m,H)   ((H-h)/m+L,H)
%                   __________
%                  /          |
% (-l,h)   (0,h) /            |
%  ____________/              |
% |                           |
% |          C| (0,0)         |
% |          C|               |
% |____________               |
%               \             |
%                 \           |
%                   \_________|
%
% everything is symmetric.  parameters are l, L, h, H, m
% r is the radius of the semi-circle
%
%

r = 1;
m = 1;
l = 3;
L = 2;
h = 2;
H = 8;
P = 2*l + 2*L + 2*h + 2*H + 2*(H-h)*sqrt(1+1/m^2);
% total perimeter

modes = (-N/2:N/2-1)';
N1 = round(2*h/P*N);
N2 = round(l/P*N);
N3 = round(L/P*N);
N4 = round(2*H/P*N);
N5 = round(N-N1-2*N2-2*N3-N4);

if (N5/2 == round(N5/2))
  N5 = N5/2;
else
  N5 = (N5-1)/2;
  N4 = N4 + 1;
end


x1 = -l + (0:N2-1)/N2*l;
y1 = h*ones(1,N2);
x2 = 0 + (0:N5-1)/N5*(H-h)/m;
y2 = h + (0:N5-1)/N5*(H-h);
x3 = (H-h)/m + (0:N3-1)/N3*L;
y3 = H*ones(1,N3);
x4 = (H-h)/m+L*ones(1,N4);
y4 = H - (0:N4-1)/N4*2*H;
x5 = (H-h)/m+L - (0:N3-1)/N3*L;
y5 = -H*ones(1,N3);
x6 = (H-h)/m - (0:N5-1)/N5*(H-h)/m;
y6 = -H + (0:N5-1)/N5*(H-h)/m;
x7 = 0 - (0:N2-1)/N2*l;
y7 = -h*ones(1,N2);
x8 = -l*ones(1,N1);
y8 = -h + (0:N1-1)/N1*2*h;
xx = [x1 x2 x3 x4 x5 x6 x7 x8];
yy = [y1 y2 y3 y4 y5 y6 y7 y8];

thresh = N;
z1 = xx + 1i*yy;
z1h = fftshift(fft(z1));
z1h(abs(modes) > thresh) = 0;
z1 = ifft(ifftshift(z1h));
% end of the smoothed exterior to the geometry

P = pi*r + 2*r;
N1 = round((pi*r)/P*N);
N2 = N - N1;
theta = linspace(pi/2,3*pi/2,N1+1); theta = theta(1:end-1);
x1 = r*cos(theta);
y1 = r*sin(theta);
x2 = zeros(1,N2);
y2 = -r + (0:N2-1)/N2*2*r;
xx = [x1 x2];
yy = [y1 y2];
% semi-circle

theta = (0:N-1)*2*pi/N;
x1 = r*cos(theta);
y1 = r*sin(theta);
xx = [x1];
yy = [y1];
% full circle

z2 = xx + 1i*yy;
z2h = fftshift(fft(z2));
z2h(abs(modes) > thresh) = 0;
z2 = ifft(ifftshift(z2h));
% end of smoothed deflector

x = [real(fliplr(z1))' real(fliplr(z2))'];
y = [imag(fliplr(z1))' imag(fliplr(z2))'];

end % diffuserGeom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function z = fourierFilter(z,IK,freq)
% z = fourierFilter(z,IK,freq) thresholds all the frequencies of f
% whose frequency is greater than freq.  IK stores the frequencies

zh = fft(z);
zh(abs(IK) > freq) = 0;
zh(end/2+1) = 0;
z = ifft(zh);

end % fourierFilter



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD ROUTINES FROM tstep.m  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [accept,dtScale] = newTimeStepSize2(o,...
%      aFinal,lFinal,aInit,lInit,accept,om)
% [accept,dtScale] = newTimeStepSize2(...
%      aFinal,lFinal,aInit,lInit,accept,om)
% finds a new time step size based on the change in the area and length
% from the first to final Gauss- Lobatto points.  The output is a flag
% indicating acceptance or rejection, and the amount the time step is
% scaled.  To increase the likelihood that the a new time step size is
% accepted, the time step size is never increased if the previous time
% step was rejected.  The only difference ebtween this routine and the
% routine newTimeStepSize is that this one looks at the global error
% commited thus far and allows continously adjusts the time step size so
% that it is aiming for the desired global tolerance.  newTimeStepSize
% always tries to commit the same amount of local truncation error per
% time step regardless of the error commited thus far

alpha = o.alpha;
betaUp = o.betaUp;
betaDown = o.betaDown;
T = o.finalTime; % time horizon
t = o.currentTime; % current time

errArea = abs(aInit - aFinal);
errLength = abs(lInit - lFinal);

tauArea = max(aInit,aFinal)*o.dt/(T-t).* ...
    (o.rtolArea*T - abs(aInit - om.area)./om.area);
tauLength = max(lInit,lFinal)*o.dt/(T-t).* ...
    (o.rtolLength*T - abs(lInit - om.length)./om.length);
err = max(max(errArea./tauArea),max(errLength./tauLength));

actualOrder = o.expectedOrder;
% order of time stepping method
dtOPT = err^(-1/(actualOrder))*o.dt;
% optimal time step size

dtOld = o.dt;
if accept
  o.dt = alpha^(1/actualOrder) * ...
      min(betaUp*o.dt,max(dtOPT,betaDown*o.dt));
else
  o.dt = alpha^(1/(actualOrder)) * ...
      min(o.dt,max(dtOPT,betaDown*o.dt));
  % don't want to scale up this time step if it was
  % previously rejected.  This hopefully gets rid of
  % pattern where a solution alternates between being
  % accepted and rejected.
end
% safety factor added to the optimal time step size also, time step size
% is not scaled up or down too fast For safety factor, take 1/p root.
% In our formulation, this makes alpha the desired value for err
% regardless of the order of the method.
dtScale = o.dt/dtOld;
% time step scaling

if err > 1
  accept = false;
  % reject time step because the error is too large
  message = ['Time Step REJECTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
  om.writeMessage(' ','%s\n')
%  fprintf('\n')
else
  accept = true;
  % accept the solution because the error is small
  message = ['Time Step ACCEPTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
end


end % newTimeStepSize2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function val = denMatVec(o,eta,walls,NearW2W)
% val = denMatVec(eta,vesicle,walls,NearV2V,NearV2W) does the matvec
% multiply required to find the density function of a given solid wall
% configuration and boundary condition

Nbd = walls.N; % Number of points per wall
nvbd = walls.nv; % Number of walls
valWalls = zeros(2*Nbd,nvbd);

op = o.op;

etaM = zeros(2*Nbd,nvbd);
for k = 1:nvbd
  istart = (k-1)*2*Nbd + 1;
  iend = k*2*Nbd;
  etaM(:,k) = eta(istart:iend);
end
otlets = eta(2*nvbd*Nbd+1:end);

% START OF EVALUATING WALL TO WALL INTERACTIONS
if ~o.fmmDLP
  kernel = @op.exactStokesDL;
else
  kernel = @op.exactStokesDLfmm;
end
FDLPwall2wall = kernel(walls,etaM,[]);
% compute the velocity on the solid walls due to the solid walls without
% near-singular integration
% END OF EVALUATING WALL TO WALL INTERACTIONS


% START OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS
if nvbd > 1
  LetsWalls = zeros(2*Nbd,nvbd);
  for k = 2:nvbd
    stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
    rotlet = otlets(3*(k-1));

    LetsWalls = LetsWalls + o.RSlets(walls.X,walls.center(:,k),...
        stokeslet,rotlet);
    % compute velocity due to rotlets and stokeslets on the solid walls
  end
  valLets = o.letsIntegrals(otlets,etaM,walls);
  % Integral constraints on the density function eta related to the
  % weights of the stokeslets and rotlets
else
  LetsWalls = zeros(2*Nbd,nvbd);
  valLets = [];
end
% END OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS


valWalls = valWalls - 1/2*etaM + ...
    op.exactStokesDLdiag(walls,o.wallDLP,etaM);
valWalls(:,1) = valWalls(:,1) + ...
    op.exactStokesN0diag(walls,o.wallN0,etaM(:,1));
% self solid wall interaction
% evaluate velocity on solid walls due to the density function.

valWalls = valWalls + FDLPwall2wall;
% velocity on walls due to all other walls
valWalls = valWalls + LetsWalls;
% velocity on walls due to the rotlets and stokeslets
% END OF EVALUATING VELOCITY ON WALLS

val = [valWalls(:);valLets];


end % denMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function val = sigDenMatVecOld(o,sigma,vesicle,walls)
% val = sigDenMatVec(sigma,vesicle,walls,NearV2V,NearV2W) does the
% matvec multiply required to find the tension and density function of a
% given vesicle configuration and farfield or solid wall velocity field

N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles
if o.confined
  Nbd = walls.N; % Number of points per wall
  nvbd = walls.nv; % Number of walls
else
  Nbd = 0;
  nvbd = 0;
end

paddedSigma = zeros(3*N*nv + 2*Nbd*nvbd + 3*(nvbd-1),1); 
% initialize matvec multiplication including room for dynamic equation
% so that we can use TimeMatVec.  Then we can simply extract the part
% that we want

for k = 1:nv
  paddedSigma((k-1)*3*N+2*N+1:k*3*N) = sigma((k-1)*N+1:k*N);
end
% put in zeros corresponding to position so that TimeMatVec doesn't add
% any new contributions due to the vesicle positions

paddedSigma(3*N*nv+1:end) = sigma(nv*N+1:end);
% add the density function, stokeslets, and rotlets
% Now is ready to pass to TimeMatVec

vesves = o.vesves;
o.vesves = 'implicit';
% want interactions to be implicit so that the most accurate tension and
% density functions are found
inextens = o.solver;
o.solver = 'method2';
dt = o.dt;
o.dt = 1;

op = o.op;
paddedVal = o.TimeMatVec(paddedSigma,vesicle,walls);
% Do a matvec but let the incoming postitions be zero since they are
% handled by the initial condition

o.vesves = vesves;
o.solver = inextens;
o.dt = dt;
% change back to old vesicle-vesicle and vesicle-boundary interactions

valTen = zeros(N*nv,1);
for k = 1:nv
  valTen((k-1)*N+1:k*N) = paddedVal((3*k-1)*N+1:3*k*N);
end
% part that corresponds to the tension
valDen = paddedVal(3*nv*N+1:3*nv*N + 2*Nbd*nvbd);
% part that corresponds to the density function
valRS = paddedVal(3*nv*N+2*Nbd*nvbd+1:end);
% part that corresponds to the Stokeslets and Rotlets

if any(vesicle.viscCont ~= 1) && o.confined
  valDen = valDen + o.viscVelocity(sigma,vesicle,walls);
  % add in contribution coming from viscosity contrast evaluated along
  % the solid walls
end

val = [valTen;valDen;valRS];
% stack the different components coming from the inextensibility, solid
% walls, and rotlets/stokeslets

end % sigDenMatVecOld

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function viscVel = viscVelocity(o,densities,vesicle,walls)
% function val = o.viscVelocity(densities,vesicle,walls) computes the
% velocity due to the viscosity contrast along the solid walls.  This
% is the extra term that is not in TimeMatVec that arises when doing
% the Schur complement to eliminate the position x.

N = vesicle.N; nv = vesicle.nv;
Nbd = walls.N; nvbd = walls.nv;
op = o.op;

sigma = zeros(N,nv);
% tension
eta = zeros(2*Nbd,nvbd);
% density along solid walls
RS = zeros(3,(nvbd-1));
% rotlets and stokeslets
for k = 1:nv
  istart = (k-1)*N+1;
  iend = k*N;
  sigma(:,k) = densities(istart:iend);
end
for k = 1:nvbd
  istart = nv*N + (k-1)*2*Nbd + 1;
  iend = istart + 2*Nbd - 1;
  eta(:,k) = densities(istart:iend);
end
% unstack the input density into terms corresponding to the tension and
% terms corresponding to the density function
for k = 2:nvbd
  istart = nv*N + nvbd*2*Nbd +(k-2)*3 + 1;
  iend = istart + 2;
  RS(:,k) = densities(istart:iend);
end

if ~o.fmm
  kernel = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
end
% kernel for single-layer potential.  Only difference is if the FMM is
% used or not

%f = vesicle.tensionTerm(sigma);
f = vesicle.tracJump(zeros(2*N,nv),sigma);
% traction nump due only to the tension

if ~o.near
  Fslp = op.exactStokesSL(vesicle,f,[]);
else
  SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
  SLPtrap = SLP;
  kernelDirect = kernel;
  Fslp = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
      o.NearV2V,kernel,kernelDirect,vesicle,true);
end
% compute the single-layer potential on all other vesicles due to the
% tension term sigma

if ~o.fmmDLP
  kernel = @op.exactStokesDL;
else
  kernel = @op.exactStokesDLfmm;
end
% kernel for double-layer potential.  Only difference is if the FMM is
% used or not

if ~o.near
  [~,FDLPvesicle] = op.exactStokesDL(walls,eta,[],vesicle.X,(1:nvbd));
  FDLPvesicle1 = FDLPvesicle;
else
  jump = -0.5; % jump condition of double-layer potential
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(walls,o.wallDLP,X);
  DLPtrap = DLP;
  kernelDirect = kernel;
  FDLPvesicle = op.nearSingInt(walls,eta,DLP,DLPtrap,...
      o.NearW2V,kernel,kernelDirect,vesicle,false);
end
% evaluate the double-layer potential due to the density function on the
% solid walls along all of the vesicles

velTen = Fslp + op.exactStokesSLdiag(vesicle,o.Galpert,f);
% add in the contribution from the self interaction

velRS = zeros(2*N,nv);
for k = 2:nvbd
  velRS = velRS + ...
      o.RSlets(vesicle.X,walls.center(:,k),RS(1:2,k),RS(3,k));
end
% add in contribution from the stokeslets and rotlets

if o.profile
  tic
end
vesVel = o.solveIminusD(velTen + FDLPvesicle + velRS,vesicle);
if o.profile
  fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
end
% Apply the operator inv(alpha*I - D).  vesVel is the velocity of the
% vesicle which is required to compute the velocity due to the viscosity
% contrast along the solid walls

if ~o.near
  [~,viscVel] = op.exactStokesDL(vesicle,vesVel,[],walls.X,(1:nv));
else
  jump = 0.5*(1-vesicle.viscCont);
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
  DLPtrap = DLP;
  kernelDirect = @op.exactStokesDL;
  viscVel = op.nearSingInt(vesicle,vesVel,DLP,DLPtrap,...
      o.NearV2W,kernel,kernelDirect,walls,false);
end
% compute the double-layer potential due to all the 
viscVel = viscVel(:);

end % viscVelocity


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function val = newPreco(o,vesicle,x,Ac,C,Schurc)
% val = newPreco(vesicle,x,Ac,C,Schurc) will apply the preconditioner
% [Bc 0; Q Sc] where Bc and Sc are the bending and schur complement
% operators defined on the circle, and Q=-Div*SLP*Ben.  Having a few
% problems getting consistency with the analytic approach (Fourier) and
% using the matricies built by vesicle.computeDerivs

N = vesicle.N;
rad = vesicle.length/2/pi;
alpha = -o.dt*vesicle.kappa/8/rad^3;
scaling = 1./(1-2*alpha*[(0:N/2-1) (N/2:-1:1)]'.^3);

bx = x(1:2*N);
bsig = x(2*N+1:3*N);

x = pinv(Ac)*bx;

bschur = bsig - vesicle.surfaceDiv(...
    o.Galpert*vesicle.bendingTerm(x));

sig = pinv(Schurc)*bschur;

val = [x;sig];


end % newPreco


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function vel = tracersVelOnlyWall(o,walls,eta,RS,Xtra)
% vel = tracersVelOnlyWall(X,sigma,kappa,walls,eta,RS,Xtra) computes
% the velocity at the set of points Xtra due to the vesicle, viscosity
% contrast, solid walls, and stokeslets and rotlets.  The particles
% Xtra are treated passively

tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;
% build a structure for tracers which mimics the class capsules.
% Can't use capsules unless there is an even number of tracer
% points (because of FFTs), so we don't want to build the class
% as we do for vesicle and walls

[~,NearW2T] = walls.getZone(tracers,2);
% build near-singular integration structures for vesicle to tracer
% and wall to tracer interactions

Ntra = size(Xtra,1)/2; % Number of tracers
Nbd = walls.N; % Number of points on walls
nvbd = walls.nv; % Number of components to walls
vel = zeros(2*Ntra,1);

op = poten(Nbd);

if ~o.fmmDLP
  kernel = @op.exactStokesDL;
else
  kernel = @op.exactStokesDLfmm;
end

% START OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY
if ~o.near
  [~,Fwall2Tra] = kernel(walls,eta,Xtra,1:nvbd);
else
  DLP = o.wallDLP;
  jump = -1/2;
  for k = 1:nvbd
    DLP(:,:,k) = DLP(:,:,k) + jump*eye(2*Nbd);
  end
  kernelDirect = @op.exactStokesDL;
  Fwall2Tra = op.nearSingInt(walls,eta,DLP,DLP,...
      NearW2T,kernel,kernelDirect,tracers,false);
end
% compute the velocity on the tracers due to the
% solid walls

FLets2Tra = zeros(2*Ntra,1);
for k = 2:nvbd
  stokeslet = RS(1:2,k);
  rotlet = RS(3,k);
  FLets2Tra = FLets2Tra + o.RSlets(Xtra,walls.center(:,k),...
    stokeslet,rotlet);
end
% velocity due to the stokeslets and rotlets
% END OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY



vel = Fwall2Tra + FLets2Tra;
% velocity is the summy of velocities due to vesicles, solid walls,
% rotlet and stokeslets, and viscosity contrast

%vel = Fwall2Tra + FLets2Tra;
% DEBUG: FOR TESTING TRACERS WITH NO VESICLES
%r2 = Xtra(1).^2 + Xtra(2).^2;
%speed = -1/3 + 400/3/r2;
%velTrue = speed*[-Xtra(2);Xtra(1)];
% This is the true velocity for the couette apparatus when both
% boundaries are centered at the origin, the inner boundary rotates
% once every 2*pi time units and the outer boundary is fixed



end % tracersVelOnlyWalls



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function checkErrors(o,Xprov,sigmaProv,...
%    vesicle,Galpert,vesicleNew,Gnew,...
%    deltaX,deltaSigma,residual)
% checkErrors(Xprov,sigmaProv,...
%    vesicle,Galpert,vesicleNew,Gnew,...
%    deltaX,deltaSigma,residual)
% can be used to check how much error is left over in the residual when
% doing SDC

N = vesicle(1).N;
nv = vesicle(1).nv;
kappa = vesicle(1).kappa;
viscCont = vesicle(1).viscCont;
Ben = zeros(2*N,2*N,nv,o.orderGL);
Ten = zeros(2*N,N,nv,o.orderGL);
BenNew = zeros(2*N,2*N,nv,o.orderGL);
TenNew = zeros(2*N,N,nv,o.orderGL);

op = poten(N);
for k = 1:o.orderGL
  [Ben(:,:,:,k),Ten(:,:,:,k),~] = vesicle(k).computeDerivs;
  [BenNew(:,:,:,k),TenNew(:,:,:,k),~] = ...
      vesicleNew(k).computeDerivs;
end

integrand = zeros(2*N,nv,o.orderGL);
for k = 1:o.orderGL
  integrand(:,:,k) = integrand(:,:,k) - ...
      kappa*(Gnew(:,:,:,k)*BenNew(:,:,:,k) - ...
       Galpert(:,:,:,k)*Ben(:,:,:,k))*Xprov(:,:,k);
  integrand(:,:,k) = integrand(:,:,k) + ...
      (Gnew(:,:,:,k)*TenNew(:,:,:,k) - ...
       Galpert(:,:,:,k)*Ten(:,:,:,k))*sigmaProv(:,:,k);
  integrand(:,:,k) = 0*integrand(:,:,k) + ...
      (-kappa*Gnew(:,:,:,k)*BenNew(:,:,:,k)*deltaX(:,:,k) + ...
        Gnew(:,:,:,k)*TenNew(:,:,:,k)*deltaSigma(:,:,k));
end


integral = o.lobattoInt(integrand);
err = residual + integral - deltaX;
for k = o.orderGL:o.orderGL
  disp([norm(err(:,:,k),inf) norm(deltaX(:,:,k),inf)])
end


end % checkErrors
