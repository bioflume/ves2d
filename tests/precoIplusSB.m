%clear all
addpath ../src
addpath ../examples
oc = curve;

N = 128;
theta = (0:N-1)'*2*pi/N;
[x,y] = geom(theta);
X = [x;y];
%X = oc.initConfig(N);
%X = oc.initConfig(N,'curly');
%X = oc.initConfig(N,'openStar',5,'amplitude',0);
%X = oc.initConfig(N,'openStar',20,'amplitude',0.4);
%X = [cos(theta);sin(theta)];
[ra,area,length] = oc.geomProp(X);
radius = length/2/pi;
%radius = sqrt(area/pi);
Xc = radius*[cos(theta);sin(theta)];

vesicle = capsules(X,[],[],1,1,false);
vesiclec = capsules(Xc,[],[],1,1,false);
op = poten(N);

order = 0;
eps = 1e-1;

G = op.stokesSLmatrix(vesicle);
Gc = op.stokoesSLmatrix(vesiclec);
%G = op.laplaceSLmatrix(vesicle);
%Gc = op.laplaceSLmatrix(vesiclec);
%Gc = kron(eye(2),Gc);
%if 1
%  Greg = op.laplaceRegularizedSLmatrix(vesicle,eps,order);
%else
%  Greg = op.FourthDerivLogRegularizedMatrix(vesicle,eps);
%end
%Greg = kron(eye(2),Greg);
% laplace SLP

Greg = zeros(N);
modes = (-N/2:N/2-1)';
fnct = 0;
%for j = 1:N
%  disp(j)
%  for k = 1:N
%    [Greg(k,j),nn] = quad(@(theta) regularizedKernel(theta,...
%        @(theta) exp(1i*modes(j)*theta),X(k),X(k+N),eps,order),...
%        0,2*pi,1e-8);
%    fnct = fnct + nn;
%  end
%end
%Greg = real(Greg * fft1.fourierInt(N));

%G = kron(eye(2),G);
Greg = kron(eye(2),Greg);
Ben = vesicle.computeDerivs;

dt = 1e-1;
tol = 1e-8;

A = eye(2*N) + dt*G*Ben;
Ac = eye(2*N) + dt*Gc*Ben;
Areg = eye(2*N) + dt*Greg*Ben;

rhs = X;
xTrue = A\rhs;

if 1
  [x0,flag,relres,iter,relresvec0] = gmres(A,rhs,[],tol,2*N);
  iterUnpreco = iter(2);
  fprintf('***************************************************\n')
  fprintf('Unpreconditioned:\n')
  fprintf('GMRES requires %d iterations\n',iter(2));
  res = A*x0 - rhs;
  err = xTrue - x0;
  fprintf('Norm of residual is %4.2e\n',...
      norm(res)/norm(rhs));
  fprintf('Norm of error is    %4.2e\n',...
      norm(err)/norm(xTrue));
  fprintf('***************************************************\n\n')
end


if 1
  P = inv(Areg);
  preco = @(X) P*X;
  [x1,flag,relres,iter,relresvec1] = gmres(A,rhs,[],...
      tol,2*N,preco);
  iterUnpreco = iter(2);
  fprintf('***************************************************\n')
  fprintf('Regularized preconditioned:\n')
  fprintf('GMRES requires %d iterations\n',iter(2));
  res = A*x1 - rhs;
  err = xTrue - x1;
  fprintf('Norm of residual is %4.2e\n',...
      norm(res)/norm(rhs));
  fprintf('Norm of error is    %4.2e\n',...
      norm(err)/norm(xTrue));
  fprintf('***************************************************\n\n')
end

if 0
  P = inv(Ac);
  preco = @(X) P*X;
  [x2,flag,relres,iter,relresvec2] = gmres(A,rhs,[],...
      tol,2*N,preco);
  iterUnpreco = iter(2);
  fprintf('***************************************************\n')
  fprintf('Circle preconditioned:\n')
  fprintf('GMRES requires %d iterations\n',iter(2));
  res = A*x2 - rhs(1:2*N);
  err = xTrue - x2;
  fprintf('Norm of residual is %4.2e\n',...
      norm(res)/norm(rhs));
  fprintf('Norm of error is    %4.2e\n',...
      norm(err)/norm(xTrue));
  fprintf('***************************************************\n\n')
end

fprintf('Average quadrature function evalutions per call is %4.2e\n',...
  fnct/N^2);


