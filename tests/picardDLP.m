% look at using picard smoothing for dlp
% Need to do it in blocks

addpath ../examples
addpath ../src

prams.Nbd = 64;
prams.nvbd = 1;

options.farField = 'cylinder';
options.farField = 'choke';
options.confined = true;
[options,prams] = initVes2D(options,prams);

oc = curve;
%Xwalls = oc.initConfig(prams.Nbd,...
%    options.farField,'nv',prams.nvbd,'center',[[0;0] [0;0]]);
Xwalls = oc.initConfig(prams.Nbd,...
    options.farField,'nv',prams.nvbd,...
    'center',[0;0],'scale',1);

walls = capsules(Xwalls,[],[],0,0);

tt = tstep(options,prams);
op = poten(prams.Nbd);
tt.wallDLP = op.stokesDLmatrix(walls);
tt.wallN0 = op.stokesN0matrix(walls);

preco = tt.wallsPrecond(walls);
forwardMat = inv(preco); 
% build the matricies that we need
forwardMat(1:2*prams.Nbd,1:2*prams.Nbd) = ...
    forwardMat(1:2*prams.Nbd,1:2*prams.Nbd) - tt.wallN0(:,:,1);
% Don't want the extra N0 term.  It mucks up the spectrum of the
% smoother.  Hopefully, it will only need to be used in the coarse grid
% solve

U = tt.farField(walls.X);
U = rand(size(U));
% couette boundary condition
rhs = [U(:);zeros(3*(prams.nvbd-1),1)];
etaExact = preco*rhs;

lambda = -1/2;
splitLeft = zeros(size(forwardMat));
splitLeft(1:2*prams.Nbd*prams.nvbd,1:2*prams.Nbd*prams.nvbd) =  ...
    lambda*eye(2*prams.Nbd*prams.nvbd);
%splitLeft(1:2*prams.Nbd,1:2*prams.Nbd) = lambda*eye(2*prams.Nbd);
%splitLeft(2*prams.Nbd+1:4*prams.Nbd,2*prams.Nbd+1:4*prams.Nbd) = ...
%    lambda*eye(2*prams.Nbd);
splitLeft(4*prams.Nbd+1:end,4*prams.Nbd+1:end) = ...
      -2*pi*eye(3*(prams.nvbd-1));
%splitLeft(4*prams.Nbd+1:end,4*prams.Nbd+1:end) = ...
%      lambda*eye(3*(prams.nvbd-1));

splitRight = forwardMat - splitLeft;

eta = rhs;
iterMat = inv(splitLeft)*splitRight;
err = eta - etaExact;
plot(err)
pause

for k = 1:10
  eta = splitLeft\(rhs - splitRight*eta);
  err = eta - etaExact;
  plot(err)
  pause
end






