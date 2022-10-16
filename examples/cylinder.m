%clear all; %clc

disp('One elliptical vesicles in a cylinder.');
disp('First-order semi-implicit time stepping.');
disp('Implicit vesicle-boundary interactions.');

format short e
prams.N = 32;                  % points per vesicle
prams.nv = 1;                  % number of vesicles
prams.T = 1e-1;                 % time horizon
prams.m = 1;                  % number of time steps
prams.Nbd = 128;                % number of points on solid wall
prams.Nbdcoarse = 32;          % coarse grid
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = 1e0;          % viscosity contrast
prams.gmresTol = 1e-12;        % GMRES tolerance
prams.errorTol = 5e-1;
% Maximum relative error in area and length before the simulation
% is stopped


options.farField = 'cylinder';    % Constricted domain
%options.farField = 'relaxation';
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 1;          % time stepping order
options.inextens = 'method2'; 
% Inextensibility condition can be treated as either
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.vesves = 'implicit';
% Discretization of vesicle-vesicle and vesicle-boundary 
% intearctions.  Either 'explicit' or 'implicit'
% Boundary-boundary interactions are always handled 
% implicitly.  The preconditioner for the boundary-boundary
% interactions is exact
options.collision = false;
options.near = true;        % near-singular integration
options.fmm = false;
% use FMM to compute single-layer potential
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = true;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-11 11 -11 11];
% Axis for plots if usePlot = true
options.saveData = true;
% Save vesicle information and create a log file
options.logFile = 'output/cylinder.log';
% Name of log file
options.dataFile = 'output/cylinderData.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.profile = false;    % Profile code
options.timeAdap = false;
options.orderGL = 2;
options.nsdc = 0;



[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
Xwalls = oc.initConfig(prams.Nbd,...
    options.farField,'nv',prams.nvbd,...
    'center',[0;0],'scale',0.5);
% Build solid walls if they exist

%X = oc.initConfig(prams.N,'center',[0;0]);
theta = (0:prams.N-1)'*2*pi/prams.N;
X = [cos(theta);1.5*sin(theta)];
%X = [cos(theta)+8.99;1.5*sin(theta)];

Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code


