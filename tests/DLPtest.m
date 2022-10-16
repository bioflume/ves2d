clear all; clc;
% March 16, 2016

% This test routine is to compare the performances of the new (Manas) and the
% old (Zydrunas) FMMs.

prams.N = 12
prams.nv = 36
prams.T = 0.1;
prams.m = 1;
prams.kappa = 1e-1;
prams.viscCont = 10*ones(1,prams.nv);

options.farField = 'taylorGreen';
options.order = 1;
options.vesves = 'implicit';
options.near = true;
options.fmm = true;
options.fmmDLP = true;
options.orderGL = 2;
options.saveData = false;

options.repulsion = false;
options.nsdc = 0;
options.antiAlias = true;
options.timeAdap = false;

if 0
oc = curve;
angle = pi/2*ones(prams.nv,1);
scale = 0.5;
X = oc.initConfig(prams.N,'nv',prams.nv,'scale',scale,...
    'reducedArea',0.65,...
    'angle',angle,...
    'center',[[-6 -1 4 9];[0 0 0 0]]);
%Initial configuration
end

if 1
addpath ../examples/
load IC_36VesTG.mat
end

if 0
addpath ../examples/  
oc = curve;
load ICcouetteVF5.mat
prams.nv = size(X,2);
prams.Nbd = 256;
prams.Nbdcoarse = 256;
prams.nvbd = 2;
options.confined = true;
options.farField = 'couette';
Xwalls = oc.initConfig(prams.Nbd,options.farField,'nv',prams.nvbd,'center',...
  [[0;0] [0;0]]);
end
[options,prams] = initVes2D(options,prams);

op = poten(prams.N,options.fmmPrecision,12*prams.N);

kernelOld = @op.exactStokesDLfmm;
kernelNew = @op.exactStokesDLnewfmm;

Nup = prams.N * ceil(sqrt(prams.N));

if prams.N ~= 12
  X = [interpft(X(1:end/2,:),prams.N);interpft(X(end/2+1:end,:),prams.N)];
end
Xup = [interpft(X(1:prams.N,:),Nup);...
       interpft(X(prams.N+1:2*prams.N,:),Nup)];

vesicleUp = capsules(Xup,[],[],prams.kappa,prams.viscCont,options.antiAlias);

if options.confined
vesicle = capsules(X,[],[],prams.kappa,prams.viscCont,options.antiAlias);
walls = capsules(Xwalls,[],[],zeros(prams.nvbd,1),zeros(prams.nvbd,1),options.antiAlias);
end

tic
farFieldOld = kernelOld(vesicleUp,Xup);
t1 = toc;
tic
farFieldNew = kernelNew(vesicleUp,Xup);
t2 = toc;

if options.confined
tic 
[~,wallFarOld] = kernelOld(vesicle,X,Xwalls,(1:prams.nv));
t3 = toc;
tic
[~,wallFarNew] = kernelNew(vesicle,X,Xwalls,(1:prams.nv));
t4 = toc;
end

%[farFieldOld(1:Nup/prams.N:end,:) farFieldNew(1:Nup/prams.N:end,:)]
err = norm(farFieldOld-farFieldNew)/norm(farFieldOld)
if options.confined
errWall = norm(wallFarOld-wallFarNew)/norm(wallFarOld)
end

disp(['Old FMM: ' num2str(1/t1) '. New FMM: ' num2str(1/t2)])
if options.confined
disp(['Old FMM: ' num2str(1/t3) '. New FMM: ' num2str(1/t4)])
end
