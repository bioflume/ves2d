clear all; %clc

disp('Vesicles in a pinched cylinder with tangential boundary condition.');

prams.N = 16;                 % points per vesicle
prams.nv = 2;                  % number of vesicles
prams.T = 10;                  % time horizon
prams.m = 1000;                 % number of time steps
prams.Nbd = 192;                % number of points on solid wall
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = 1e0;          % viscosity contrast
prams.gmresTol = 1e-10;        % GMRES tolerance
prams.errorTol = 5e-1;
% Maximum relative error in area and length before the simulation
% is stopped


options.farField = 'figureEight';    % Constricted domain
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 1;          % time stepping order
options.inextens = 'method1'; 
% Inextensibility condition can be treated as either
options.vesves = 'implicit';
% Discretization of vesicle-vesicle and vesicle-boundary 
options.near = true;        % near-singular integration
options.fmm = false;
options.fmmDLP = false;
% use FMM to compute single-layer and double-layer potentials
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = true;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-2 2 -2 2];
% Axis for plots if usePlot = true
options.saveData = true;
% Save vesicle information and create a log file
options.logFile = 'output/figureEight.log';
% Name of log file
options.dataFile = 'output/figureEightData.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.profile = false;    % Profile code
options.tracers = false;      % trace passive particles in flow
% For now, there will be a bug when using near-singular integration to
% compute the pressure at a point interior to a vesicle

options.orderGL = 2;
% number of Gauss-Lobatto points to use in the substep.  Letting


% ADD-ONS
options.nsdc = 0;
options.correctShape = true;
options.alignCenterAngle = true;
options.reparameterization = true;
prams.maxReparamIter = 50;

options.repulsion = true;
prams.minDist = 0.2;
prams.minSpecRatio = 90;
prams.repStrength = 900;

options.timeAdap = true;
prams.rtolArea = 1e-4;
prams.rtolLength = 1e-4;

prams.dtMax = 2;
prams.dtMin = 1e-4;


[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
if options.confined
  Xwalls = oc.initConfig(prams.Nbd,options.farField);
end
% Build solid walls



X = oc.initConfig(prams.N,'nv',prams.nv,...
    'scale',0.06,'center',[[0.2;0.2] [-0.2;-0.2]],...
    'reducedArea',0.65);
%centers = [0.2 -0.2 -1 1 -1.2  1.2 -1.2 1.2; ...
%           0.2 -0.2  0 0 -0.5 -0.5  0.5 0.5];
%X = oc.initConfig(prams.N,'nv',prams.nv,...
%    'scale',0.12,...
%    'center',centers,...
%    'reducedArea',0.65);
%[X,prams.nv] = oc.initConfig(prams.N,'volFrac',0.15,...
%    Xwalls,options.fmm,'scale',0.12);
% Initial vesicles has reduced area 0.65 and is a vertical ellipse
% so that it must distort to pass through the contriction

Xtra = [];
pressTar = [];
% target points to compute the pressure

Xfinal = Ves2D(X,Xwalls,prams,options,Xtra,pressTar);
% Run vesicle code


