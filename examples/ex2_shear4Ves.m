clear; clc;
restoredefaultpath;

fprintf('Four elliptical vesicles in a shear flow.\n');
options.farField = 'shear'; % background velocity
options.farFieldSpeed = 1; % scaling of background velocity
prams.nv = 4;    % number of vesicles

% Spatial-resolution
prams.N = 32;    % points per vesicle

% physical properties of vesicles
prams.kappa = 1e-1;   % bending coefficient
prams.viscCont = 1*ones(prams.nv,1);   % viscosity contrast

% parameters for numerics
options.fmm = false;  % fmm for single-layer potentials
options.fmmDLP = false; % fmm for double-layer potentials

prams.gmresTol = 1e-10;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

% Axis for the plot
% options.axis = [-6 6 -5 5];

% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA
options.repulsion = true; % repulsion between vesicles


% Temporal Resolution (parameters for new implementation)
prams.T = 20;   % time horizon (or 100*lenScale/options.farFieldSpeed) etc.

options.timeAdap = true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 5e-3; % tolerance for errors in area-length for time-stepping

options.usePlot = true; % Plot on-the-fly
options.track = false; % trackers on membrane

% Save vesicle information and create a log file
% Name of log file for saving messages
options.logFile = 'output/shear4Ves.log';
% Name of binary data file for storing vesicle information
options.dataFile = 'output/shear4VesData.bin';

% Set options and parameters that the user doesn't
% Also add src to path
[options,prams] = initVes2D(options,prams);

% Initial configuration of vesicles
oc = curve;
angle = pi/2*ones(prams.nv,1);
scale = 0.5;
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'scale',scale,...
    'reducedArea',0.65,...
    'angle',angle,...
    'center',[[ -3 -1 1 3];[ 0 0 0 0 ]]);

% Get the length and velocity scale of simulation to decide on params.
[~,~,arcLen] = oc.geomProp(X); 
lenScale = max(arcLen); % length scale of the simulation

% time scale for simulation
% # of time steps (if constant time stepping) or Dt for the first step
% (adaptive)
prams.m = ceil(prams.T/(0.01*lenScale/options.farFieldSpeed)); 
prams.dtMax = 5*lenScale/options.farFieldSpeed; % maximum allowable time-step size
prams.dtMin = 1e-4*lenScale/options.farFieldSpeed; % minimum allowable time-step size

% Run vesicle code
Xfinal = Ves2D(X,[],[],[],prams,options,[],[]);


