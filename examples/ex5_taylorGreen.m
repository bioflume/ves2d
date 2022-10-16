clear; clc;
restoredefaultpath;

fprintf('Vesicles in a Taylor-Green flow.\n');
options.farField = 'taylorGreen'; % background velocity
options.farFieldSpeed = 1; % scaling of background velocity
prams.nv = 9;    % number of vesicles

% Spatial-resolution
prams.N = 16;    % points per vesicle

% physical properties of vesicles
prams.kappa = 1e-1;   % bending coefficient
prams.viscCont = 10*ones(prams.nv,1);   % viscosity contrast

% parameters for numerics
options.fmm = false;  % fmm for single-layer potentials
options.fmmDLP = false; % fmm for double-layer potentials

prams.gmresTol = 1e-10;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

% Axis for the plot
options.axis = [0 pi 0 pi] + [-0.1 0.1 -0.1 0.1];

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
options.logFile = 'output/taylorGreen9Ves.log';
% Name of binary data file for storing vesicle information
options.dataFile = 'output/taylorGreen9VesData.bin';

% Set options and parameters that the user doesn't
% Also add src to path
[options,prams] = initVes2D(options,prams);

% Initial configuration of vesicles
oc = curve;
% ############# 9 VESICLES ######################
sx = [.02 -.01 -.014 .12 .04 -.08 .08 .02 -.06]; 
sy = [-.02 -.03 -.011 .04 .02 .02 .03 .06 -.02]; 
cenx = kron(ones(1,3),[pi/4 pi/2 3*pi/4]);
ceny = kron([pi/4 pi/2 3*pi/4],ones(1,3));
scale = 0.15;

% ############# 16 VESICLES ######################
% sx = -0.015 + 0.03*rand(1,16);
% sy = -0.015 + 0.03*rand(1,16);
% cenx = kron(ones(1,4),(1:4)*pi/5);
% ceny = kron((1:4)*pi/5,ones(1,4));
% scale = 0.145;

% ############# 36 VESICLES ######################
% sx = -0.025 + 0.05*rand(1,36);
% sy = -0.025 + 0.05*rand(1,36);
% cenx = kron(ones(1,6),(1:6)*pi/7);
% ceny = kron((1:6)*pi/7,ones(1,6));
% scale = 0.09;

% ############# 81 VESICLES ######################
%sx = -0.02 + 0.04*rand(1,81);
%sy = -0.02 + 0.04*rand(1,81);
%cenx = kron(ones(1,9),(1:9)*pi/10);
%ceny = kron((1:9)*pi/10,ones(1,9));
%scale = 0.065;

cenx = cenx + sx;
ceny = ceny + sy;
% x and y coordinates of centers
% Have added pertubation to centers so that flow is 
% more interesting
angle = ones(prams.nv,1);

X = oc.initConfig(prams.N,'nv',prams.nv,...
  'reducedArea',0.65,...
  'angle',angle,...
  'center',[cenx;ceny],...
  'scale',scale);


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

