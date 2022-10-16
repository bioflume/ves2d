clear all; clc

VolFrac = [4 5 6 7];
ReducedArea = [65 80];

for kk=1:length(VolFrac)
  for jj=1:length(ReducedArea)
    fname  = 'couetteFlow_';
    fname = [fname num2str(VolFrac(kk)) '_' num2str(ReducedArea(jj)) '_200_1_64'];
    
    load(fname);
    fname = [fname '.bin'];
    fileId = fopen(fname,'r');
    Result = fread(fileId,'double');
    fclose(fileId);
    
    Result = reshape(Result,5*nv*n+1+2*sum(prams.M)+3,[]); 
    Time = Result(5*nv*n+1,:);


    idx = find(diff(Time) < 0);
    Result = Result(:,1:idx);
    
    fileId = fopen(fname,'w');
    Result = fwrite(fileId,Result,'double');
    fclose(fileId);
  end
end

% addpath ../src/

% getNv =@(volFrac, area, size) round(volFrac * area / pi / size^2);

% reducedArea = .8;             
% vesDist = ''; 
% vesSize = 1;         
% volFrac = .01;
% Ri = 10; Ro = 20;             


% prams.nv = getNv(volFrac, pi * (Ro^2-Ri^2), vesSize);
% prams.n = 64;                     
% prams.T = 200;                               
% prams.ts = 1e-3;                              
% prams.m = prams.T/prams.ts;                   
% prams.order = 1;                    
% prams.viscCont = ones(1,prams.nv);                           
% prams.kappa = 1e-1;                            

% fileName = ['./couetteFlow_' num2str(prams.nv) '_' num2str(100*reducedArea) ...
%             '_' num2str(prams.T) '_' num2str(1000*prams.ts) '_' ...
%             num2str(prams.n)];

% options.usePlot = 0;                          
% options.progressBar = 0;                      
% options.saveData = 1;                         
% options.dataStride = ceil(prams.m/1000);      
% options.fileName = [fileName '.bin'];         
% options.useGPU = false;                             
% options.showError = true;                     
% options.GmresTol = 1e-8;                     
% options.verbose = true;

% %%-- Setting up the fixed boundary
% Ri = 10; Ro = 20;             
% omegaIn = 1.0; omegaOut = 0;        
% prams.flowType = 'confined';             
% prams.M = [256 128];                     
% prams.bd = @(ind,m) sampleBd(ind,m,1,'couette','Ri',Ri,'Ro',Ro);
% prams.bc = @(x) forcing(x,'couette','Ri',Ri,'Ro',Ro,'omegaIn',omegaIn,'omegaOut',omegaOut);;
% prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams,'direct',options.useGPU);

% %%-- Generating the vesicles
% domain = fixedBound(prams.M,prams.bd,1);
% X = boundary(prams.n,'couette',domain,vesDist,'nv',prams.nv,'scale',vesSize,'reducedArea',reducedArea);
% X = makeGridUniform(X);

% %%-- Checking 
% [trash prams options] = monitor(prams,options);
% viewer(X,[],prams,options);
             
% save([fileName '.mat']);
% [XF status] = Ves2D(X,prams,options,@monitor);