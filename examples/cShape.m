clear all; clc

disp('Two vesicles in a relaxation flow.');
disp('Implicit vesicle-vesicle interactions.');

prams.N = 96;                     % points per vesicle
prams.nv = 2;                     % number of vesicles
prams.T = 1000;                     % time horizon
prams.m = 10000;                    % number of time steps
prams.kappa = 1e-1;               % bending coefficient
prams.viscCont = 10;               % viscosity contrast
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.gmresTol = 1e-10;           % GMRES tolerance
prams.rtolArea = 1e-1;            % allowable error in area
prams.rtolLength = 1e-1;          % allowable error in length
%prams.gmresMaxIter = 400;

options.farField = 'relaxation';   % background velocity
options.order = 1;          % time stepping order
% ONLY HAVE INCLUDED FIRST AND SECOND ORDER
options.inextens = 'method1'; 
% Inextensibility condition can be treated as either
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.vesves = 'implicit';
% Discretization of vesicle-vesicle intearctions.
% Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.fmm = false;      
options.fmmDLP = false;
% use FMM to compute single-layer potential
options.usePlot = true;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-3 3 -3 3]; % Axis for plots if usePlot=true
options.saveData = true;    
% Save vesicle information and create a log file
options.logFile = 'output/cShape.log';
% Name of log file
options.dataFile = 'output/cShapeData.bin';
% Name of binary data file for storing vesicle information
options.resFile = 'output/cShapeRes.dat';
% Name of dat file for storing residual
options.verbose = true;
% Decides how much info is written to console
options.profile = false;    % Profile code
options.timeAdap = true;  % Adaptive time stepping
options.orderGL = 2;       % Number of Gauss-Lobatto points
options.nsdc = 1;
options.collision = false;  % Collision detection
options.pressure = false;   % compute the pressure
options.redistributeArcLength = true;
options.correctShape = true;


[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'nv',1,...
    'cShape');
theta = (0:prams.N-1)'*2*pi/prams.N;
X = [X [0.105+0.5*sin(theta);0.6*cos(theta)]];

Xwalls = [];

Xfinal = Ves2D(X,Xwalls,prams,options);
% run vesicle simulation
