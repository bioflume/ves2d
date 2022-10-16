clear; clc;
restoredefaultpath;

fprintf('Multiple elliptical vesicles in a couette apparatus.\n');
options.farField = 'couette'; % background velocity
options.farFieldSpeed = 1; % scaling of background velocity
prams.nv = 125;    % number of vesicles
options.confined = true; % confined or unbounded geometry

% Spatial-resolution
prams.N = 96;    % points per vesicle
prams.Nbd = 256; % points on a wall
prams.nvbd = 2;  % number of solid walls

% physical properties of vesicles
prams.kappa = 1e-1;   % bending coefficient
prams.vesViscCont = 1;   % viscosity contrast

% parameters for numerics
options.fmm = ~false;  % fmm for single-layer potentials
options.fmmDLP = ~false; % fmm for double-layer potentials
options.matFreeWalls = false; % W2W interactions are done without a matrix

prams.gmresTol = 1e-8;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

% Axis for the plot
options.axis = [-20 20 -20 20];

% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA
options.repulsion = true; % repulsion between vesicles


% Temporal Resolution (parameters for new implementation)
prams.T = 100;   % time horizon (or 100*lenScale/options.farFieldSpeed) etc.

options.timeAdap = true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 2e-3; % tolerance for errors in area-length for time-stepping

options.usePlot = true; % Plot on-the-fly
options.track = false; % trackers on membrane

% Save vesicle information and create a log file
% Name of log file for saving messages
options.logFile = 'output/couetteRA75_2.log';
% Name of binary data file for storing vesicle information
options.dataFile = 'output/couetteRA75Data_2.bin';

% Set options and parameters that the user doesn't
% Also add src to path
[options,prams] = initVes2D(options,prams);

% Initial configuration of vesicles
oc = curve;
% Build solid walls
Xwalls = oc.initConfig(prams.Nbd,...
    options.farField,'nv',prams.nvbd,'center',[[0;0] [0;0]]);

% SINGLE VESICLE
if 0
  prams.nv = 1;
  X = oc.initConfig(prams.N,'nv',prams.nv,'angle',pi/2,...
    'center',[15;0],'scale',0.6);
end

% TWO VESICLES 
if 0
  prams.nv = 2;
  centers = [[15;0] [-15;0]];
  X = oc.initConfig(prams.N,'nv',prams.nv,'center',centers,'scale',0.5);
end

% INPUT VOLUME FRACTION OF VESICLES
% Randomly place a template vesicle until the desired volume
% fraction is achieved.  Uses collision detection to decide if they
% a randomly placed vesicle intersects a solid wall or another
% vesicle
if 0
  [X,prams.nv] = oc.initConfig(prams.N,'volFrac',0.01,...
     Xwalls,options.fmm,options,prams);
end

% 70% VOLUME FRACTION (260 vesicles)
if 1
  % Generate the vesicles 
  scale = 0.45;
  nOuter = ceil(0.41*prams.nv);
  nMiddle = ceil(0.33*prams.nv);
  nInner = prams.nv - nOuter - nMiddle;

  omega = (0:nOuter - 1)*2*pi/nOuter;
  center = 18.3*[cos(omega);sin(omega)];
  angle1 = omega;

  omega = (0:nMiddle - 1)*2*pi/nMiddle;
  center = [center 15*[cos(omega);sin(omega)]];
  angle2 = omega;

  omega = (0:nInner - 1)*2*pi/nInner;
  center = [center 11.7*[cos(omega);sin(omega)]];
  angle3 = omega;

  angle = [angle1 angle2 angle3];
  angle = angle+(-0.2+0.4*rand(1,prams.nv));

  X = oc.initConfig(prams.N,'nv',prams.nv,'angle',angle,...
      'scale',scale,'center',center,...
      'reducedArea',0.75);
end
prams.viscCont = prams.vesViscCont*ones(prams.nv,1);


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
Xfinal = Ves2D(X,Xwalls,[],[],prams,options,[],[]);



