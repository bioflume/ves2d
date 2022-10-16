function randomShearRuns(runName,N,Th,dt,speed,VC,kappa,freqRange,numRandFuncs)
rng('shuffle');
fprintf('Simple elliptical vesicle in a shear flow.\n');
options.farField = 'random_shear'; % background velocity
options.farFieldSpeed = speed; % scaling of background velocity
prams.nv = 1;    % number of vesicles

options.saveVinf = true;

% Spatial-resolution
prams.N = N;    % points per vesicle

% physical properties of vesicles
prams.kappa = kappa;   % bending coefficient
prams.viscCont = VC;   % viscosity contrast

% parameters for numerics
options.fmm = ~true;  % fmm for single-layer potentials
options.fmmDLP = ~true; % fmm for double-layer potentials

prams.gmresTol = 1e-10;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA

% Temporal Resolution (parameters for new implementation)
prams.T = Th;   % time horizon (or 100*lenScale/options.farFieldSpeed) etc.

options.timeAdap = ~true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 5e-3; % tolerance for errors in area-length for time-stepping


options.usePlot = false; % Plot on-the-fly
options.track = false; % trackers on membrane

% Save vesicle information and create a log file
% Name of log file for saving messages
options.logFile = ['./output/randomShears/' runName '.log'];
% Name of binary data file for storing vesicle information
options.dataFile = ['./output/randomShears/' runName 'Data.bin'];
prams.runName = runName;

% Set options and parameters that the user doesn't
% Also add src to path
[options,prams] = initVes2D(options,prams);

% Initial configuration of vesicles
oc = curve;
load relaxed64.mat
X = [interpft(X(1:end/2),prams.N);interpft(X(end/2+1:end),prams.N)];
% X = oc.initConfig(prams.N,'reducedArea',0.65,'angle',pi/2);
% Get the length and velocity scale of simulation to decide on params.
% Get the inclination angle
IA = oc.getIncAngle(X);
IAnew = pi*rand;

% Decenter and rotate
cy = -0.8 + 1.6*rand;
% cy = 0;
X = [X(1:end/2)-mean(interpft(X(1:end/2),256));...
    X(end/2+1:end)-mean(interpft(X(end/2+1:end),256))];
Xnew = zeros(size(X));
Xnew(1:end/2) = cos(-IA)*X(1:end/2)-sin(-IA)*X(end/2+1:end);
Xnew(end/2+1:end) = sin(IA)*X(1:end/2)+cos(-IA)*X(end/2+1:end);

X(1:end/2) = cos(IAnew)*Xnew(1:end/2)-sin(IAnew)*Xnew(end/2+1:end);
X(end/2+1:end) = sin(IAnew)*Xnew(1:end/2)+cos(IAnew)*Xnew(end/2+1:end)+cy;


[~,~,arcLen] = oc.geomProp(X); 
lenScale = max(arcLen); % length scale of the simulation

% time scale for simulation
% # of time steps (if constant time stepping) or Dt for the first step
% (adaptive)
prams.m = prams.T/dt;
prams.dtMax = 5*lenScale/options.farFieldSpeed; % maximum allowable time-step size
prams.dtMin = 1e-4*lenScale/options.farFieldSpeed; % minimum allowable time-step size

om = monitor(X,options,prams);
tt = tstep(options,prams,om);


tt.farField = @(X,Xint) tt.bgFlow(X,options.farField,'vInf',...
    [],'Speed',options.farFieldSpeed,'numRandFuncs',...
    numRandFuncs,'freqRange',freqRange);
% Run vesicle code
Xfinal = Ves2D(X,[],[],[],prams,options,tt,[]);
end


