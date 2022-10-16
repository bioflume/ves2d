clear; clc;
restoredefaultpath;

fprintf('Simple elliptical vesicle in a shear flow.\n');
options.farField = 'shear'; % background velocity
options.farFieldSpeed = 4; % scaling of background velocity
prams.nv = 1;    % number of vesicles

options.saveVinf = ~true;

% Spatial-resolution
prams.N = 32;    % points per vesicle

% physical properties of vesicles
prams.kappa = 1e-3;   % bending coefficient
prams.viscCont = 10;   % viscosity contrast

% parameters for numerics
options.fmm = false;  % fmm for single-layer potentials
options.fmmDLP = false; % fmm for double-layer potentials

prams.gmresTol = 1e-10;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA

% Temporal Resolution (parameters for new implementation)
prams.T = 20;   % time horizon (or 100*lenScale/options.farFieldSpeed) etc.

options.timeAdap = ~true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 5e-3; % tolerance for errors in area-length for time-stepping


options.usePlot = true; % Plot on-the-fly
options.track = false; % trackers on membrane

% Save vesicle information and create a log file
% Name of log file for saving messages
options.logFile = 'output/shear1VesTB2.log';
% Name of binary data file for storing vesicle information
options.dataFile = 'output/shear1VesTB2Data.bin';
prams.runName = 'shear1VesTB';

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

% Decenter and rotate
cy = 0;
X = [X(1:end/2)-mean(interpft(X(1:end/2),256));...
    X(end/2+1:end)-mean(interpft(X(end/2+1:end),256))];
Xnew = zeros(size(X));
Xnew(1:end/2) = cos(pi/2-IA)*X(1:end/2)-sin(pi/2-IA)*X(end/2+1:end);
Xnew(end/2+1:end) = sin(pi/2-IA)*X(1:end/2)+cos(pi/2-IA)*X(end/2+1:end)+cy;
X = Xnew;

[~,~,arcLen] = oc.geomProp(X); 
lenScale = max(arcLen); % length scale of the simulation

% time scale for simulation
% # of time steps (if constant time stepping) or Dt for the first step
% (adaptive)
prams.m = ceil(prams.T/(0.01*lenScale/options.farFieldSpeed)); 
prams.m = prams.T/0.01;
prams.dtMax = 5*lenScale/options.farFieldSpeed; % maximum allowable time-step size
prams.dtMin = 1e-4*lenScale/options.farFieldSpeed; % minimum allowable time-step size

om = monitor(X,options,prams);
tt = tstep(options,prams,om);

Xfinal = Ves2D(X,[],[],[],prams,options,tt,[]);
% Run vesicle code

