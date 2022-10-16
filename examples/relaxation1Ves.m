clear; clc;
restoredefaultpath;

disp('Simple curly vesicle in a relaxation flow.');
options.farField = 'relaxation'; % background velocity
prams.nv = 1; % number of vesicles

% Spatial-resolution
prams.N = 64;    % points per vesicle

% physical properties of vesicles
prams.kappa = 1e-3;   % bending coefficient
prams.viscCont = 1*ones(prams.nv,1);   % viscosity contrast

% parameters for numerics
options.fmm = false;  % fmm for single-layer potentials
options.fmmDLP = false; % fmm for double-layer potentials

prams.gmresTol = 1e-10;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

% Axis for the plot
options.axis = [-0.5 0.5 -0.5 0.5]; % Axis for plots if usePlot=true

% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = ~true; % reparametrization of membrane
options.equiDistArcLength = true;
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA


% Temporal Resolution (parameters for new implementation)
options.timeAdap = false; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 5e-3; % tolerance for errors in area-length for time-stepping

options.usePlot = true; % Plot on-the-fly
options.track = false; % trackers on membrane

% Save vesicle information and create a log file
% Name of log file for saving messages
options.logFile = 'output/relaxation1Ves.log';
% Name of binary data file for storing vesicle information
options.dataFile = 'output/relaxation1VesData.bin';

% Set options and parameters that the user doesn't
% Also add src to path
[options,prams] = initVes2D(options,prams);

% Initial configuration of a vesicle
oc = curve;

% Curly shape
% X = oc.initConfig(prams.N,'curly');
% Open star
% X = oc.initConfig(prams.N,'openStar',8,'amplitude',0);
load vesShapesv7.mat
X = [interpft(Xshapes(1:end/2,200),prams.N);interpft(Xshapes(end/2+1:end,200),prams.N)];

% Get the length and velocity scale of simulation to decide on params.
[~,~,arcLen] = oc.geomProp(X); 
lenScale = max(arcLen); % length scale of the simulation

% time scale for simulation
% # of time steps (if constant time stepping) or Dt for the first step
% (adaptive)
prams.T = 20;   % time horizon (or 100*lenScale/options.farFieldSpeed) etc.
prams.m = prams.T/1e-1; 
prams.dtMax = 5*lenScale^2/(prams.kappa); % maximum allowable time-step size
prams.dtMin = 1e-4*lenScale^2/(prams.kappa); % minimum allowable time-step size

% Run vesicle code
Xfinal = Ves2D(X,[],[],[],prams,options,[],[]);

