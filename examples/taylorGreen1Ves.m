clear all; clc

fprintf('Four elliptical vesicles in a Taylor-Green flow.\n');
fprintf('Second-order semi-implicit time stepping.\n');
fprintf('Implicit vesicle-vesicle interactions.\n');

% Physics parameters
prams.N = 32;                           % points per vesicle
prams.nv = 1;                          % number of vesicles
prams.T = 500;                          % time horizon
prams.m = 1e5;                         % number of time steps
prams.kappa = 1e-1;                   % bending coefficient
prams.errorTol = 1e-1;
prams.viscCont = 100*ones(prams.nv,1); % viscosity contrast
prams.gmresTol = 1e-10;
options.farField = 'taylorGreen'; % background velocity
options.farFieldSpeed = 1.5;
% method of enforcing inextensibility.
% Can be 'method1' or 'method2'
options.order = 1;                % time stepping order
% Discretization of vesicle-vesicle intearctions.
% Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.antiAlias = true;
options.fmm = false;
options.fmmDLP = false;
% use FMM to compute single-layer potential
options.axis = [-4*pi 4*pi -4*pi 4*pi] + [-0.5 0.5 -0.5 0.5];
% axis for plot
options.logFile = 'output/taylorGreen1VesVC100.log';
% Name of log file for saving messages
options.dataFile = 'output/taylorGreen1VesVC100Data.bin';
% Name of binary data file for storing vesicle information

options.profile = false;

% ADD-ONS
options.alignCenterAngle = true;
options.correctShape = true;
options.reparameterization = true;
prams.maxReparamIter = 100;

options.repulsion = ~true;
prams.minDist = 0.3; %0.3
prams.minSpecRatio = 90; %30
prams.repStrength = 90; %90

options.timeAdap = true;
prams.rtolArea = 1e-3;
prams.rtolLength = 1e-3;

prams.dtMax = 1e-2;
prams.dtMin = 1e-4;
prams.betaInc = 1e-2;
prams.betaDec = 5e-2;
prams.betaUp = 1.5;
prams.betaDown = 0.5;
prams.alpha = 0.9;


% Plot on-the-fly
options.usePlot = ~true;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't

oc = curve;
% ############# 4 VESICLES ######################
cenx = [0.9*pi];
ceny = [0.8*pi];
scale = 1/7;
% x and y coordinates of centers
% Have added pertubation to centers so that flow is 
% more interesting
angle = pi/2;

X = oc.initConfig(prams.N,'nv',prams.nv,...
  'reducedArea',0.65,...
  'angle',angle,...
  'center',[cenx;ceny],...
  'scale',scale);

% Initial configuration

% PLOT INITIAL CONFIGURATION
% figure(1);clf;
% xvec = [X(1:end/2,:);X(1,:)];
% yvec = [X(end/2+1:end,:);X(end/2+1,:)];
% plot(xvec,yvec,'r','linewidth',2)
% axis equal
% max(yvec)-min(yvec)
% pause

Ves2D(X,[],[],[],prams,options,[]);
% Run vesicle code

