p1 = '../src';
p2 = '../src/neareval'; 
addpath(p1);
addpath(p2);
%load traj_example3
load tracers_example4.mat

XX = X(1,:);
YY = Y(1,:);
X = [posx(:,:,1) ; posy(:,:,1)];
M = 0;
fmm = 1;


N = size(X,1); N = N/2;
nv = size(X,2);
%Get number of vesicles and number of points on each vesicle
[kap2 tang2 sa2] = curveProp(X);
%Need the curvature, unit tangent vector and arclength

X1 = zeros(N*nv+ntrac,1);
X2 = zeros(N*nv+ntrac,1);
norx = zeros(N*nv+ntrac,1);
nory = zeros(N*nv+ntrac,1);
sa = zeros(N*nv+ntrac,1);
kap = zeros(N*nv+ntrac,1);
for k=1:nv
  X1((k-1)*N+1:k*N) = X(1:N,k);
  X2((k-1)*N+1:k*N) = X(N+1:2*N,k);
  norx((k-1)*N+1:k*N) = tang2(N+1:2*N,k);
  nory((k-1)*N+1:k*N) = -tang2(1:N,k);
  sa((k-1)*N+1:k*N) = sa2(1:N,k)/N;
  kap((k-1)*N+1:k*N) = kap2(1:N,k);
end
%Arrange everything in an appropriate manner and get the unit normal

X1(nv*N+1:end) = XX;
X2(nv*N+1:end) = YY;


if (fmm)
  size(norx)
  size(sa)
  size(X1)
  size(X2)
  pause
  [~,rfield,~] = fmm_laplace(norx.*sa,X1,X2,X1,X2);
  [~,~,cfield] = fmm_laplace(nory.*sa,X1,X2,X1,X2);
  pot = -(rfield + cfield - 1/2*kap.*sa);
else
  pot = zeros(N*nv + ntrac,1);
  for k = 1:N*nv + ntrac
    dist2 = (X1(k) - X1).^2 + (X2(k)-X2).^2;
    ind = [(1:k-1) (k+1:N*nv+ntrac)];
    pot(k) = -sum((((X1(k)-X1(ind)).*norx(ind) + ...
        (X2(k)-X2(ind)).*nory(ind)))./dist2(ind).*...
        sa(ind));
    pot(k) = pot(k) + 1/2*kap(k)*sa(k);
  end
end
%Two calls to the FMM and add the curavture term for the diagonal

pot = pot(N*nv+1:end);
% Only care about the tracers, not the vesicles 
