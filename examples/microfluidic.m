clear all; clc

disp('Multiple elliptical vesicles in a double-couette apparatus.');
disp('Second-order semi-implicit time stepping.');
disp('Implicit vesicle-boundary interactions.');

format long
prams.N = 96;                  % points per vesicle
prams.nv = 2;                  % number of vesicles
prams.T = 10*10;                   %time horizon
prams.m = 1000*10;                 % number of time steps
prams.Nbd = 256;                % number of points on solid wall
prams.Nbdcoarse = 256;
prams.nvbd = 5;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = 1e0;          % viscosity contrast
prams.gmresTol = 1e-10;        % GMRES tolerance
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.rtolArea = 1e-1;
prams.rtolLength = 1e-1;


options.farField = 'microfluidic';    % Constricted domain
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
% Boundary-boundary interactions are always handled 
% implicitly.  The preconditioner for the boundary-boundary
% interactions is exact
options.near = true;        % near-singular integration
options.fmm = true;        
options.fmmDLP = true;
% use FMM to compute single-layer potential
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = true;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = 2*[-1 1 -1 1]; 
% Axis for plots if usePlot = true
options.saveData = true;    
% Save vesicle information and create a log file
options.logFile = 'output/microfluidic.log';
% Name of log file
options.dataFile = 'output/microfluidicData.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.profile = false;    % Profile code
options.collision = false;   % Collision detection
options.tracers = false;
options.timeAdap = false;
options.orderGL = 2;
options.nsdc = 0;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
Xwalls = oc.initConfig(prams.Nbd,options.farField);

%X = oc.initConfig(prams.N,'center',[-1.8;0.05],...
%    'scale',0.03,'angle',pi/2,'reducedArea',0.65);
center = [-1.3 -0.3 0.5 1.5 0 0 0 0;...
           1e-2 -1e-2 1e-2 -1e-2 -1.5-1e-2 -0.5+1e-2 0.5-1e-2 1.5-1e-2];
angle = [pi/2*ones(prams.nv/2,1);zeros(prams.nv/2,1)];
center = [-0.3 0.5;...
           -1e-2 -1e-2];
angle = [pi/2*ones(prams.nv/2,1);zeros(prams.nv/2,1)];
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'center',center,'scale',0.02,...
    'angle',angle,'reducedArea',0.65);

Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code


