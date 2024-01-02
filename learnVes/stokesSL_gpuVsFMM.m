clear; clc;

load ./ICs/VF35_81VesIC.mat; % 81 vesicles loaded
% X loaded, N = 64 
nv = size(X,2);
Ntests = [32;64;128];
% Ntests = [64];
sys_size = Ntests*81;

addpath ../src/
oc = curve;

gpuTimes = zeros(numel(Ntests),1);
cpuTimes = zeros(numel(Ntests),1);
fmmTimes = zeros(numel(Ntests),1);
valCPU = zeros(numel(Ntests),1);
valGPU = zeros(numel(Ntests),1);
valFMM = zeros(numel(Ntests),1);

for it = 1 : numel(Ntests)
  disp(it)
  % Choose discretization
  N = Ntests(it);
  % upsample the vesicle resolution
  Xv = [interpft(X(1:end/2,:),N); interpft(X(end/2+1:end,:),N)];
  % build vesicle object
  vesicle = capsules(Xv, [], [], 1, 1, 0);
  % calculate bending force
  fBend = vesicle.tracJump(Xv,zeros(N,nv));

  op = poten(N);

  % Direct calculation on single core
  stokesDirect = @() exactStokesSLDirect(vesicle,fBend);
  cpuTimes(it) = timeit(stokesDirect);
  valCPU = exactStokesSLDirect(vesicle,fBend);
  
  % FMM Calculation
  stokesFMM = @() exactStokesSLfmm(vesicle,fBend);
  fmmTimes(it) = timeit(stokesFMM);
  valFMM = exactStokesSLfmm(vesicle,fBend);

  % GPU Calculation
  saGPU = gpuArray(single(vesicle.sa));
  fGPU = gpuArray(single(fBend));
  XGPU = gpuArray(single(vesicle.X));
  stokesGPUarr = gpuArray(single(zeros(2*vesicle.N,vesicle.nv)));
  stokesGPU = @() exactStokesSLGPUVect(XGPU, saGPU, fGPU, stokesGPUarr);
  gpuTimes(it) = timeit(stokesGPU);
  stokesGPUarr = gpuArray(single(zeros(2*vesicle.N,vesicle.nv)));
  valGPU = exactStokesSLGPU(XGPU,saGPU,fGPU,stokesGPUarr);
end
save('comparison.mat','fmmTimes','gpuTimes','cpuTimes')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesSLP = exactStokesSLDirect(vesicle,f)

den = f.*[vesicle.sa;vesicle.sa]*2*pi/vesicle.N;
% multiply by arclength term


stokesSLP = zeros(2*vesicle.N,vesicle.nv); % Initialize to zero

for k = 1:vesicle.nv % vesicle of targets
  K = [(1:k-1) (k+1:vesicle.nv)];
  % Loop over all vesicles except k
  for j=1:vesicle.N
    dis2 = (vesicle.X(j,k) - vesicle.X(1:vesicle.N,K)).^2 + ...
        (vesicle.X(j+vesicle.N,k) - ...
          vesicle.X(vesicle.N+1:2*vesicle.N,K)).^2;
    diffxy = [vesicle.X(j,k) - vesicle.X(1:vesicle.N,K) ; ...
        vesicle.X(j+vesicle.N,k) - ...
          vesicle.X(vesicle.N+1:2*vesicle.N,K)];
    % distance squared and difference of source and target location

    coeff = 0.5*log(dis2);
    % first part of single-layer potential for Stokes

    val = coeff.*den(1:vesicle.N,K);
    stokesSLP(j,k) = -sum(val(:));
    val = coeff.*den(vesicle.N+1:2*vesicle.N,K);
    stokesSLP(j+vesicle.N,k) = -sum(val(:));
    % logarithm terms in the single-layer potential
  
    coeff = (diffxy(1:vesicle.N,:).*den(1:vesicle.N,K) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*...
        den(vesicle.N+1:2*vesicle.N,K))./dis2;
    % second part of single-layer potential for Stokes

    val = coeff.*diffxy(1:vesicle.N,:);
    stokesSLP(j,k) = stokesSLP(j,k) + sum(val(:));
    val = coeff.*diffxy(vesicle.N+1:2*vesicle.N,:);
    stokesSLP(j+vesicle.N,k) = stokesSLP(j+vesicle.N,k) + ...
        sum(val(:));
    % r \otimes r term of the single-layer potential

  end % j
end % k
% Evaluate single-layer potential at vesicles but oneself
  
end % exactStokesSLDirect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesSLP = exactStokesSLGPUVect(X, sa, f, stokesSLP)
N = size(X,1)/2;
nv = size(X,2);
Ntar = N;
fx = f(1:end/2,:);
fy = f(end/2+1:end,:);

% multiply by arclength term

for k = 1:nv % vesicle of targets
  K = [(1:k-1) (k+1:nv)];
  % Loop over all vesicles except k
  allX = X(1:N,K); allX = allX(:);
  allY = X(N+1:2*N,K); allY = allY(:);
  allFx = fx(:,K);
  allFy = fy(:,K);
  allsa = sa(:,K);

  denx = allFx.*allsa*2*pi/N; denx = denx(:);
  deny = allFy.*allsa*2*pi/N; deny = deny(:);
  
  ddenx = denx(:,ones(Ntar,1))';
  ddeny = deny(:,ones(Ntar,1))';
  
  xsou = allX(:,ones(Ntar,1))';
  ysou = allY(:,ones(Ntar,1))';

  Nsou = numel(allX);
  xx = X(1:N,k);
  yy = X(N+1:2*N,k);
  xtar = xx(:,ones(Nsou,1));
  ytar = yy(:,ones(Nsou,1));
  
  diffx = xtar-xsou;
  diffy = ytar-ysou;
  dis2 = diffx.^2 + diffy.^2;
  
  coeff = 0.5*log(dis2);
  
  val = coeff.*ddenx;
  stokesSLP(1:N,k) = -sum(val,2);
  val = coeff.*ddeny;
  stokesSLP(N+1:2*N,k) = -sum(val,2);

  coeff = (diffx.*ddenx + diffy.*ddeny)./dis2;

  val = coeff.*diffx;
  stokesSLP(1:N,k) = stokesSLP(1:N,k) + sum(val,2);
  val = coeff.*diffy;
  stokesSLP(N+1:2*N,k) = stokesSLP(N+1:2*N,k) + sum(val,2); 
end % k

end % exactStokesSLGPUVect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stokesSLP = exactStokesSLfmm(vesicle,f)
% [stokesSLP,stokeSLPtar] = exactStokesSLfmm(vesicle,f,Xtar,K) uses the
% FMM to compute the single-layer potential due to all vesicles except
% itself vesicle is a class of object capsules and f is the density
% function NOT scaled by arclength term.  Xtar is a set of points where
% the single-layer potential due to all vesicles in index set K needs
% to be evaulated


oc = curve;
[x,y] = oc.getXY(vesicle.X); % seperate x and y coordinates

den = f.*[vesicle.sa;vesicle.sa]*2*pi/vesicle.N;

iprec = 3;

[f1,f2] = oc.getXY(den);
% need to multiply by arclength term.  Seperate it into
% x and y coordinate

[u,v] = stokesSLPfmm(f1(:),f2(:),x(:),y(:),...
      x(:),y(:),1,iprec);
stokesSLP = zeros(2*vesicle.N,vesicle.nv); % initialize
SLPself = zeros(2*vesicle.N, vesicle.nv);
for k = 1:vesicle.nv
  is = (k-1)*vesicle.N+1;
  ie = k*vesicle.N;
  stokesSLP(1:vesicle.N,k) = u(is:ie);
  stokesSLP(vesicle.N+1:2*vesicle.N,k) = v(is:ie);
 
  SLPmat = stokesSLmatrixNoCorr(vesicle.X(:,k), vesicle.sa(:,k));
  SLPself(:,k) = SLPmat * f(:,k);
end
% Wrap the output of the FMM into the usual 
% [[x1;y1] [x2;y2] ...] format
  
% we need to remove the self-interaction terms, so compute those terms
% using the self-interaction matrix

stokesSLP = stokesSLP - SLPself;

end % exactStokesSLfmm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = stokesSLmatrixNoCorr(X, sa)
% G = stokesSLmatrixNoCorr(vesicle) computes stokes single layer potential
% matrix without any quadrature. This corresponds to the way SLP FMM
% computes the SLP. Its counterpart is stokesDLmatrixNoCorr

oc = curve;
[x,y] = oc.getXY(X);
% Vesicle positions
N = size(X,1)/2;
% points per vesicle
      
sa = sa';
sa = sa(ones(N,1),:);
% Jacobian

xtar = x(:,ones(N,1))'; 
ytar = y(:,ones(N,1))'; 
% target points

xsou = x(:,ones(N,1)); 
ysou = y(:,ones(N,1));
% source points

diffx = xtar - xsou;
diffy = ytar - ysou;
rho2 = (diffx.^2 + diffy.^2).^(-1);
% one over the distance squared

logpart = 0.5*log(rho2);
% sign changed to positive because rho2 is one over distance squared

G11 = (logpart + diffx.^2.*rho2).*sa;
G12 = (diffx.*diffy.*rho2).*sa;
G22 = (logpart + diffy.^2.*rho2).*sa;
% don't need negative sign in front of log since rho2 is one over
% distance squared
G11(1:N+1:N^2) = 0;
G12(1:N+1:N^2) = 0;
G22(1:N+1:N^2) = 0;
% zero out the diagonal term

G = [G11 G12; G12 G22];
% build matrix with 4 blocks
G = 1/4/pi*G*2*pi/N; 
% scale with the arclength spacing and divide by 4*pi


end % stokesSLmatrixNoCorr

