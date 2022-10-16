clear; clc;
restoredefaultpath;

fprintf('One vesicle in a parabolic flow.\n');
options.farField = 'parabolic'; % background velocity
options.farFieldSpeed = 1; % scaling of background velocity
prams.nv = 1;    % number of vesicles

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
options.axis = [-20 20 -10 30]; 

% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA

% Temporal Resolution (parameters for new implementation)
prams.T = 50;   % time horizon (or 100*lenScale/options.farFieldSpeed) etc.

options.timeAdap = true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 5e-3; % tolerance for errors in area-length for time-stepping

options.usePlot = true; % Plot on-the-fly
options.track = false; % trackers on membrane

% Save vesicle information and create a log file
% Name of log file for saving messages
options.logFile = 'output/parabolic1Ves.log';
% Name of binary data file for storing vesicle information
options.dataFile = 'output/parabolic1VesData.bin';

% Set options and parameters that the user doesn't
% Also add src to path
[options,prams] = initVes2D(options,prams);


% Initial configuration of vesicles
angle = pi/2;
%angle = 0;
scale = 2.75e0;
oc = curve;
cy = 3*1e0;
ra = 0.75;

% For 200 vesicles
% rng('shuffle');
% angle = pi*ones(prams.nv,1);
% scale = 1/2;
% oc = curve;
% xx = 40;
% yy = 8;
% cx = linspace(-xx,xx,40);
% cy = linspace(-yy,yy,5);
% [cx,cy] = meshgrid(cx,cy);
% center = [cx(:)'+4e-1*rand(1,prams.nv);cy(:)'+4e-1*rand(1,prams.nv)];
% center = center(:,1:prams.nv);

X = oc.initConfig(prams.N,'nv',prams.nv,'angle',angle,...
    'scale',scale,'center',[0;cy],'reducedArea',ra);

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


