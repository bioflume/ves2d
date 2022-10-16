clear all;clc;
addpath ../src/

nv = 2;
n = 64;
Ri = 10;            
Ro = 20;            
omegaIn = 1;        
omegaOut = 0;       
reducedArea = .9;   
vesSize = 1;        
M = [128 64];

%%-- Simulation parameters and options
prams.T = 200;        
prams.n = n;              
prams.ts= 1e-3;
prams.kappa = 1e-1;                           
prams.flowType = 'confined';             
prams.m = prams.T/prams.ts;                               
prams.order = 2;
prams.M = M;

options.usePlot = 0;                          
options.progressBar = 0;                      
options.useGPU = 0;                           
options.showError = true;                     
options.saveData = true;
options.dataStride = floor(prams.m/100);
options.fileName = 'couetteFlow.bin';

prams.bd = @(ind,m) sampleBd(ind,m,1,'couette','Ri',Ri,'Ro',Ro);
prams.bc = @(x) forcing(x,'couette','Ri',Ri,'Ro',Ro,'omegaIn',omegaIn, ...
                        'omegaOut',omegaOut);

prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,...
                                 prams,'direct',options.useGPU);
[t prams options] = monitor(prams,options);
domain = fixedBound(prams.M,prams.bd,1);

%% Reading file
fileId = fopen(options.fileName,'r');
Result = fread(fileId,'double');
fclose(fileId);

Result = reshape(Result,5*nv*n+1+2*sum(M)+3,[]); 
Result = Result(:,end);

Xv   = reshape(Result(1       :2*nv*n,:),[],nv);
sig  = reshape(Result(2*n*nv+1:3*n*nv,:),[],nv);
u    = reshape(Result(3*n*nv+1:5*n*nv,:),[],nv); 

%% streamline grid
[xg yg] = meshgrid(-20:.1:20); 
IN = inDomain(domain,xg,yg,.2);
Xg = [xg(IN) yg(IN)];

%% vInf
%vesicle perturbation
[F FX] = InteractionForce(Xv,u,sig,prams,options,1,Xg');
FX = FX';

uu = zeros(size(xg));
vv = zeros(size(yg));
uu(IN) = FX(:,1); uu = reshape(uu,size(xg));   
vv(IN) = FX(:,2); vv = reshape(vv,size(yg));   

subplot(3,1,1);
hold on;
viewer(Xv,[],prams, options);
h = streamslice(xg,yg,uu,vv,2,'cubic'); 
axis([-20 20 -20 20]);
hold off;
title('Streamlines due to the vesicle \newline(no boundary effect)');

%true streamlines
[trash XX] = prams.vInf([],[]); XX = XX';
[F Fbg] = InteractionForce(Xv,u,sig,prams,options,1,XX);
bgFlow = reshape(prams.vInf(Xg, prams.bc(XX')-Fbg'),[],2);
FX = bgFlow+FX;

uu = zeros(size(xg));
vv = zeros(size(yg));
uu(IN) = FX(:,1); uu = reshape(uu,size(xg));   
vv(IN) = FX(:,2); vv = reshape(vv,size(yg));   

subplot(3,1,2);
hold on;
viewer(Xv,[],prams, options);
h = streamslice(xg,yg,uu,vv,2,'cubic');
axis([-20 20 -20 20]);
hold off;

title('Couette flow streamlines');
%subtracting the couette flow
unperturbed = prams.bc(reshape(Xg,[],2));
FX=FX-unperturbed;

uu = zeros(size(xg));
vv = zeros(size(yg));
uu(IN) = FX(:,1); uu = reshape(uu,size(xg));   
vv(IN) = FX(:,2); vv = reshape(vv,size(yg));   

subplot(3,1,3);
hold on;
viewer(Xv,[],prams, options);
h = streamslice(xg,yg,uu,vv,2,'cubic');
axis([-20 20 -20 20]);
hold off;
title('The variation with respect to \newline an unperturbed flow');




