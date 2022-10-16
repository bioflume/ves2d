clear all;clc;
addpath ../src/

nv = 1;
omegaOut = 0;       
reducedArea = .9;   
vesSize = 1;        

%%-- Simulation parameters and options
options.usePlot = 0;                          
options.progressBar = 0;                      
options.saveData = 0;                         
options.useGPU = 0;                           
options.showError = true;                     
prams.kappa = 1e-1;                           
prams.flowType = 'confined';             
prams.T = 10;        
prams.n = 128;              
prams.m =  1600;
prams.ts = prams.T/prams.m;                               
prams.order = 1;

RI = 100;
for jj = 1:length(RI)
  Ri = RI(jj);
  Ro = Ri+10;            
  omegaIn = 10/Ri;        
  M0 = 128*[Ro Ri]/10;
  
  prams.bd = @(ind,m) sampleBd(ind,m,1,'couette','Ri',Ri,'Ro',Ro);
  prams.bc = @(x) forcing(x,'couette','Ri',Ri,'Ro',Ro,'omegaIn',omegaIn, ...
                          'omegaOut',omegaOut);
  cc = 2.^(-3:4);
  for ii = 1:length(cc)
    clear functions;
    
    %%-- Setting up the boundary
    prams.M = cc(ii)*M0;
    prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,...
                                     prams,'direct',options.useGPU);
    [t prams options] = monitor(prams,options);
    domain = fixedBound(prams.M,prams.bd,1);
    X = boundary(prams.n,'couette',domain,'uniform','nv',nv,'scale',vesSize, ...
                 'reducedArea',reducedArea,'angle',pi/2);
    [rv ai li] = reducedVolume(X);
    
    %%-- Convergence test
    domain = fixedBound(prams.M,prams.bd,1);
    [Xf status] = Ves2D(X,prams,options,@monitor);
    disp(status);
    [rv af lf] = reducedVolume(Xf);
    
    err_A(ii,jj) = max(abs(af./ai-1));
    err_L(ii,jj) = max(abs(lf./li-1));
    save('results/boundaryEffects.mat');
  end
end
