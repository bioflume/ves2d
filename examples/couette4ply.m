clear all; clc

fprintf('Multiple elliptical vesicles in a double-couette apparatus.\n');
fprintf('First-order semi-implicit time stepping.\n');
fprintf('Explicit vesicle-boundary interactions.\n');

%prams.N = 256;                      % points per vesicle
prams.N = 128;                      % points per vesicle
prams.nv = 12;
prams.T = 100;                     % time horizon
%prams.m = 1e4;                   % number of time steps
prams.m = 1e5;                   % number of time steps
prams.Nbd = 128;                    % number of points on solid wall
%prams.Nbd = 256;                    % number of points on solid wall
prams.nvbd = 2;                    % number of components to solid walls
prams.kappa = 1e-1;                % bending coefficient
prams.errorTol = 1e1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.rtolArea = 1e12;        % maximum allowable area error
prams.rtolLength = 1e12;      % maximum allowable length error

%options.farField = 'quadCouette'; % Constricted domain
options.farField = 'couette'; % Constricted domain
%options.farField = 'cylinder'; % Constricted domain
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 1;          % time stepping order
options.inextens = 'method2'; 
% Inextensibility condition can be treated as either
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.vesves = 'implicit';
% Discretization of vesicle-vesicle intearctions.
% Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.fmm = true;
% use FMM to compute single-layer potential
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.axis = [-21 21 -21 21];
% Axis for plots if usePlot = true
options.logFile = 'output/couette4ply.log';
% Name of log file for saving messages
options.dataFile = 'output/couette4plyData.bin';
% Name of binary data file for storing vesicle information
options.timeAdap = true;
% use time adaptivity
options.orderGL = -7;
% order of internal time steps

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
Xwalls = oc.initConfig(prams.Nbd,...
    options.farField,'nv',prams.nvbd,...
    'center',[[0;0] [0;-11] [-11;0] [11;0] [0;11]]);
% Build solid walls

%[X,prams.nv] = oc.initConfig(prams.N,'volFrac',0.15,...
%    Xwalls,options.fmm);
%X = oc.initConfig(prams.N,'nv',prams.nv,...
%    'center',[[3;3] [-3;3] [3;-3] [-3;-3]]+0.1*rand(2,prams.nv),...
%    'scale',1,'angle',pi*rand(prams.nv,1));
%X = oc.initConfig(prams.N,'nv',prams.nv,...
%    'center',[4;3],'scale',1,...
%    angle',3*pi/4*ones(prams.nv,1));
% Randomly place a template vesicle until the desired volume
% fraction is achieved.  Uses collision detection to decide if they
% a randomly placed vesicle intersects a solid wall or another
% vesicle

%theta = (0:prams.N-1)'*2*pi/prams.N;
%X = [cos(theta);3*sin(theta)];

addpath output;
file = 'output/couette4plyCrashed.bin';
[posx,posy,ten,wallx,wally,ea,el,time,n,nv] = loadFile(file);
X = [posx(:,1:prams.nv,1);posy(:,1:prams.nv,1)];
%X = [interpft(posx(:,:,1),prams.N);interpft(posy(:,:,1),prams.N)];
%X = [posx(:,:,659);posy(:,:,659)];
%X = [interpft(posx(:,1:prams.nv,550),prams.N);interpft(posy(:,1:prams.nv,550),prams.N)];

Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code


