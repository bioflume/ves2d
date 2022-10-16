% Build interaction matricies for two vesicles so that we can look at
% preconditioner

clear all;
addpath ../src

dt = 1e-1; % time step size
N = 32; % points per vesicle
modes = (-N/2:N/2-1)'; % fourier frequencies
theta = (0:N-1)'*2*pi/N;

x1 = cos(theta) - 1.1; y1 = 3*sin(theta);
% vesicle 1
x2 = cos(theta) + 1.1; y2 = 3*sin(theta);
% vesicle 2

X = [[x1;y1] [x2;y2]];
vesicle = capsules(X,[],[],1,1);

op = poten(N);
% for computing single-layer potentials
SLPdiag = op.stokesSLmatrix(vesicle);
% single-layer potential where source and targets are vesicle

[Ben,Ten,Div] = vesicle.computeDerivs;
% bending, tension, and surface divergence operators of both
% vesicles


V11 = [eye(2*N) + dt*SLPdiag(:,:,1)*Ben(:,:,1) ...
                  -dt*SLPdiag(:,:,1)*Ten(:,:,1); ...
                  Div(:,:,1) ...
                  zeros(N)];
V22 = [eye(2*N) + dt*SLPdiag(:,:,2)*Ben(:,:,2) ...
                  -dt*SLPdiag(:,:,2)*Ten(:,:,2); ...
                  Div(:,:,2) ...
                  zeros(N)];
% diagonal terms of time stepping matrix

FF = fft1.fourierInt(N);
% fourier differentiation matrix

SLP12 = zeros(2*N,2*N);
SLP21 = zeros(2*N,2*N);
for k = 1:N
  % loop over fourier modes to generate the off-diagonal
  % single-layer potentials in fourier space
  den = [[cos(modes(k)*theta);zeros(N,1)] ...
         [cos(modes(k)*theta);zeros(N,1)]];
  rstokesSLP = op.exactStokesSL(vesicle,den);
  % real part of Fourier mode

  den = [[sin(modes(k)*theta);zeros(N,1)] ...
         [sin(modes(k)*theta);zeros(N,1)]];
  istokesSLP = op.exactStokesSL(vesicle,den);
  % imaginary part of Fourier mode

  SLP12(:,k) = rstokesSLP(:,1) + 1i*istokesSLP(:,1);
  SLP21(:,k) = rstokesSLP(:,2) + 1i*istokesSLP(:,2);

  den = [[zeros(N,1);cos(modes(k)*theta)] ...
         [zeros(N,1);cos(modes(k)*theta)]];
  rstokesSLP = op.exactStokesSL(vesicle,den);
  % real part of Fourier mode

  den = [[zeros(N,1);sin(modes(k)*theta)] ...
         [zeros(N,1);sin(modes(k)*theta)]];
  istokesSLP = op.exactStokesSL(vesicle,den);
  % imaginary part of Fourier mode

  SLP12(:,k+N) = rstokesSLP(:,1) + 1i*istokesSLP(:,1);
  SLP21(:,k+N) = rstokesSLP(:,2) + 1i*istokesSLP(:,2);
end

SLP12 = real(SLP12 * [FF zeros(N); zeros(N) FF]);
SLP21 = real(SLP21 * [FF zeros(N); zeros(N) FF]);
% move off-diagonal single-layer potentials to physical space

V12 = [dt*SLP12*Ben(:,:,2) -dt*SLP12*Ten(:,:,2); zeros(N,3*N)];
V21 = [dt*SLP21*Ben(:,:,1) -dt*SLP21*Ten(:,:,1); zeros(N,3*N)];
% off-diagonal components of time stepping matrix

TimeStepper = [V11 V12; V21 V22];
% time stepping matrix
Preco = [inv(V11) zeros(3*N); zeros(3*N) inv(V22)];
% block diagonal preconditioner

zExact = [vesicle.X(:,1);ones(N,1);vesicle.X(:,2);ones(N,1)];
rhs = TimeStepper*zExact;

%rhs = [vesicle.X(:,1);Div(:,:,1)*vesicle.X(:,1); ...
%       vesicle.X(:,2);Div(:,:,2)*vesicle.X(:,2)] + ...
%      dt*[-vesicle.X(1:N,1);vesicle.X(N+1:2*N,1);zeros(N,1); ...
%       -vesicle.X(1:N,2);vesicle.X(N+1:2*N,2);zeros(N,1)];
%% right-hand side coming from extensional flow


tol = 1e-16; % gmres tolerance
maxit = 6*N; % maximum number of iterations
restart = []; % restart

warning off

[z,flag,relres1,iter] = gmres(TimeStepper,rhs,restart,tol,maxit);
% solve with unpreconditioned gmres
fprintf('\n*****************UNPRECONDITIONERD GMRES****************\n');
fprintf('Conditioner # of linear system is %8.4e\n',...
    cond(TimeStepper));
fprintf('Residual is                       %8.4e\n',...
    norm(TimeStepper*z-rhs)/norm(rhs))
fprintf('Error is                          %8.4e\n',...
    norm(z-zExact)/norm(zExact))
fprintf('Number of Iterations is           %d\n',...
    iter(2));
if (flag == 3) fprintf('GMRES STAGNATED!!!!!!\n'); end
fprintf('\n\n')

[z,flag,relres2,iter] = gmres(TimeStepper,rhs,restart,tol,maxit,...
    @(x) Preco*x);
% solve with preconditioned gmres
fprintf('*****************PRECONDITIONERD GMRES*****************\n');
fprintf('Conditioner # of preconditioned linear system is %8.4e\n',...
    cond(Preco*TimeStepper));
fprintf('Conditioner # of preconditioner is               %8.4e\n',...
    cond(Preco));
fprintf('Residual is                                      %8.4e\n',...
    norm(TimeStepper*z-rhs)/norm(rhs))
fprintf('Error is                                         %8.4e\n',...
    norm(z-zExact)/norm(zExact))
fprintf('Number of Iterations is                          %d\n',...
    iter(2));
if (flag == 3) fprintf('GMRES STAGNATED!!!!!!\n'); end
fprintf('\n\n')


%e1 = eig(TimeStepper);
%e2 = eig(Preco*TimeStepper);

%figure(1); clf;
%plot(real(e1),imag(e1),'r.')
%title('EigenValues of Unpreconditioned System')
%figure(2); clf;
%plot(real(e2),imag(e2),'r.')
%title('EigenValues of Preconditioned System')
%

warning on
