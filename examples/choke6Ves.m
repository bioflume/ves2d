clear all; clc

disp('Six elliptical vesicles in a constricted tube .');
disp('First-order semi-implicit time stepping.');
disp('Implicit vesicle-vesicle interactions.');
disp('Implicit vesicle-boundary interactions.');

format long
prams.N = 128;                 % points per vesicle
prams.nv = 6;                  % number of vesicles
prams.T = 8;                   % time horizon
prams.m = 3200;                % number of time steps
prams.Nbd = 256;               % number of points on solid wall
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = 1;            % viscosity contrast
prams.gmresTol = 1e-12;        % GMRES tolerance
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped


options.farField = 'choke';   
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 1;          % time stepping order
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
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = false;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-10.1 10.1 -3.1 3.1]; % Axis for plots if usePlot=true
options.saveData = true;    
% Save vesicle information and create a log file
options.logFile = 'output/choke6Ves.log';
% Name of log file
options.dataFile = 'output/choke6VesData.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.profile = false;    % Profile code
options.collision = false;   % Collision detection

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

spacing = 0.1;
xmax = -3;
xmin = xmax - (1+spacing)*prams.nv/2 + 1;
% parameters for the initial configuration

scale = 0.8;
angle = 0.4*pi;
oc = curve;
[cenx,ceny] = meshgrid(xmin:1+spacing:xmax,[-1.5 1.5]);
cenx(1:2:end) = cenx(1:2:end) - 1/2;
cen = [cenx(:)';ceny(:)'];

X = oc.initConfig(prams.N,'nv',prams.nv,'angle',angle,...
    'scale',scale,'center',cen,'reducedArea',0.65);

Xfinal = Ves2D(X,prams,options);
% Run vesicle code



