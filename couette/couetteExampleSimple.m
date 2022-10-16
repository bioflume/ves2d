clear functions global;
clc;
addpath ../src/

n = 64;              %Number of discretization points on each vesicle                    
nv = 2;              %Number of vesicles
Ri = 10;             %Couette flow inner circle radius
Ro = 20;             %Couette flow outer circle radius
omegaIn = 1;         %Inner circle angular velocity
omegaOut = 0;        %Outer cylinder angular velocity
reducedArea = .8;    %Vesicles' reduced area
volFrac = .01;       %The volume fraction (overrides nv if vesDist is set to volFrac)
vesDist = 'uniform'; %The distribution of vesicles, 'uniform' or 'random' or 'volFrac'
vesSize = 1;         %Nondimensional size of the vesicle \sqrt(A/pi)
M = [256 128];       %Number of discretization point on the boundary [outer inner]
                     
%%-- Simulation parameters and options
prams.T = 10;                                 %Simulation time horizon
prams.ts = .04;                               %Time step
prams.m = prams.T/prams.ts;                   %Number of time steps       
prams.kappa = 1e-1;                           %Bending modulus of the vesicles
prams.order = 1;                              %Time-stepping order

options.usePlot = 1;                          %Whether to plot real-time
options.progressBar = 0;                      %Show progress
options.saveData = 0;                         %Save data as a binary file
options.dataStride = ceil(prams.m/200);       %Frequency of saving data
options.fileName = './couetteFlow.bin';       %Data file name
options.useGPU = 0;                             
options.showError = true;                     %Showing the error in area.

%%-- Setting up the boundary
prams.flowType = 'confined';             
prams.M = M;                     
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette','Ri',Ri,'Ro',Ro);
prams.bc = @(x) forcing(x,'couette','Ri',Ri,'Ro',Ro,'omegaIn',omegaIn,'omegaOut',omegaOut);;
prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,...
                                 prams,'direct',options.useGPU);

[t prams options] = monitor(prams,options);
%%-- Generating the distribution of vesicles
centers = [0 12;13 9];
domain = fixedBound(prams.M,prams.bd,1);
X = boundary(n,'couette',domain,vesDist,volFrac,'nv',nv,'scale',vesSize,'reducedArea',reducedArea,'center',centers);
viewer(X,[],prams,options);
nv = size(X,2);

%%-- Time simulation
[XF status] = Ves2D(X,prams,options,@monitor);

