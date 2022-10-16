clear all; clc

disp('Restart a flow mid way through a simulation')
file = '/Users/quaife/Ves2D/docs/papers/sdc/makeplots/extensional/sdc2/m2400/extensional2VesData.bin';
file = 'output/extensional2VesDataRef.bin';

addpath output;

format long

[posx,posy,ten,wallx,wally,ea,el,time,n,nv] = loadFile(file);

prams.N = n;                   % points per vesicle
prams.nv = nv;                 % number of vesicles
prams.nvbd = 0;
prams.Nbd = 0;
prams.T = 12;                     % time horizon
prams.m = 1200;                   % number of time steps
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = 1;            % viscosity contrast
prams.gmresTol = 1e-12;         % GMRES tolerance
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped


options.farField = 'extensional';
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 2;          % time stepping order
options.inextens = 'method1'; 
% Inextensibility condition can be treated as either
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.vesves = 'implicit';
% Discretization of vesicle-vesicle and vesicle-boundary 
% intearctions.  Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.fmm = true;        
% use FMM to compute single-layer potential
options.confined =  false;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = true;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-3 3 -3 3]; % Axis for plots if usePlot=true
options.saveData = true;    
% Save vesicle information and create a log file
options.logFile = 'output/Restart.log';
% Name of log file
options.dataFile = 'output/RestartData.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.profile = false;    % Profile code
options.collision = false;   % Collision detection
options.timeAdap = false;
options.orderGL = 2;
options.nsdc = 0;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

X = [posx(:,:,2001);posy(:,:,2001)];

Xwalls = [wallx;wally];

Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code



