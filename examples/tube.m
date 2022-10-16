clear; clc;
restoredefaultpath;

disp('One elliptical vesicles in a constricted tube.');
options.farField = 'tube'; % background velocity
options.farFieldSpeed = 1.2; % scaling of background velocity
prams.nv = 1;    % number of vesicles
options.confined = true; % confined or unbounded geometry

options.putBackOrigin = true;

% Spatial-resolution
prams.N = 32;    % points per vesicle
prams.Nbd = 128; % points on a wall
prams.nvbd = 1;  % number of solid walls

% physical properties of vesicles
prams.kappa = 1e-3;   % bending coefficient
prams.viscCont = 2*ones(prams.nv,1);   % viscosity contrast

% parameters for numerics
options.fmm = false;  % fmm for single-layer potentials
options.fmmDLP = false; % fmm for double-layer potentials
options.matFreeWalls = true; % W2W interactions are done without a matrix

prams.gmresTol = 1e-10;  % tolerance for gmres
prams.errorTol = 1; % maximum error (used for constant time stepping)

% Axis for the plot
options.axis = [-1.5 1.5 -2 2]; 

% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA
options.repulsion = ~true; % repulsion between vesicles


% Temporal Resolution (parameters for new implementation)
prams.T = 200;   % time horizon (or 100*lenScale/options.farFieldSpeed) etc.

options.timeAdap = true; % time adaptivity
options.nsdc = 0;  % number of sdc corrections
prams.areaLenTol = 1e-3; % tolerance for errors in area-length for time-stepping

options.usePlot = true; % Plot on-the-fly
options.track = false; % trackers on membrane

% Save vesicle information and create a log file
% Name of log file for saving messages
options.logFile = 'output/tube1Ves.log';
% Name of binary data file for storing vesicle information
options.dataFile = 'output/tube1VesData2.bin';

% Set options and parameters that the user doesn't
% Also add src to path
[options,prams] = initVes2D(options,prams);

% Initial configuration of vesicles
% Initial vesicles has reduced area 0.65 and is a vertical ellipse
% so that it must distort to pass through the contriction
oc = curve;
load relaxed64.mat
X = [interpft(X(1:end/2),prams.N);interpft(X(end/2+1:end),prams.N)];
IA = oc.getIncAngle(X);
Xnew = zeros(size(X));
for i = 1 : prams.N
  Xnew(i) = cos(-IA)*X(i)-sin(-IA)*X(i+prams.N);
  Xnew(i+prams.N) = sin(-IA)*X(i)+cos(-IA)*X(i+prams.N);
end

X = [Xnew(1:end/2)-mean(Xnew(1:end/2));Xnew(end/2+1:end)-mean(Xnew(end/2+1:end))-1.5];

% X = oc.initConfig(prams.N,'nv',prams.nv,'angle',pi/2,...
%    'scale',0.75,...
%    'center',[[-6;0]],'reducedArea',0.65);
% centers = [-9 -5 -4.7 -4.4 -4.1 -3.8 -3.5 -3.2 -2.9 -2.6 -2.3 ...
%              -4.9 -4.3 -3.7 -4.8 -4.2 -3.6; ....
%            0 1 -0.8 0.9 -0.9 0.8 -1 1.1 -1.2 0.9 -0.8 ...
%               -2.3 -2.3 -2.3 2.3 2.3 2.3];
% X = oc.initConfig(prams.N,'nv',prams.nv,'angle',pi/2*ones(prams.nv,1),...
%     'scale',0.15,...
%     'center',centers,'reducedArea',0.65);


% Walls configuration
Xwalls = oc.initConfig(prams.Nbd,options.farField);

% Get the length and velocity scale of simulation to decide on params.
[~,~,arcLen] = oc.geomProp(X); 
lenScale = max(arcLen); % length scale of the simulation

% time scale for simulation
% # of time steps (if constant time stepping) or Dt for the first step
% (adaptive)
prams.m = ceil(prams.T/(0.01*lenScale/options.farFieldSpeed)); 
prams.dtMax = 0.5*lenScale/options.farFieldSpeed; % maximum allowable time-step size
prams.dtMin = 1e-4*lenScale/options.farFieldSpeed; % minimum allowable time-step size

% Run vesicle code
Xfinal = Ves2D(X,Xwalls,[],[],prams,options,[],[]);



