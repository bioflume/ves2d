function driver_rotational_vortex(runName, VC, RA, kappa, Vinf, vortexSize, Xg, IA) 

% Physics parameters
prams.N = 96;                           % points per vesicle
prams.nv = 1;                          % number of vesicles
prams.T = 20000;  % time horizon
prams.m = 2e7;                         % number of time steps
prams.kappa = kappa;                   % bending coefficient
prams.errorTol = 1e-1;
prams.viscCont = VC*ones(prams.nv,1); % viscosity contrast
prams.vortexSize = vortexSize;
prams.gmresTol = 1e-10;
options.farField = 'taylorCouette'; % background velocity
options.farFieldSpeed = Vinf;
% method of enforcing inextensibility.
% Can be 'method1' or 'method2'
options.order = 1;                % time stepping order
% Discretization of vesicle-vesicle intearctions.
% Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.antiAlias = true;
options.fmm = false;
options.fmmDLP = false;
options.track = true;
% use FMM to compute single-layer potential

options.logFile = [runName '.log'];
% Name of log file for saving messages
options.dataFile = [runName 'Data.bin'];
% Name of binary data file for storing vesicle information

options.profile = false;
options.axis = [-15 15 -15 15];

% ADD-ONS
options.alignCenterAngle = ~true;
options.correctShape = true;
options.reparameterization = ~true;
prams.maxReparamIter = 5;

options.repulsion = ~true;
prams.minDist = 0.3; %0.3
prams.minSpecRatio = 90; %30
prams.repStrength = 90; %90

options.timeAdap = true;
prams.rtolArea = 1e-3;
prams.rtolLength = 1e-3;

prams.dtMax = 1e-3;
prams.dtMin = 1e-5;
prams.betaInc = 1e-2;
prams.betaDec = 5e-2;
prams.betaUp = 1.5;
prams.betaDown = 0.5;
prams.alpha = 0.9;


% Plot on-the-fly
options.usePlot = ~false;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't

oc = curve;
% ############# 4 VESICLES ######################
cenx = Xg; % From left to 
ceny = 0;

% x and y coordinates of centers
% Have added pertubation to centers so that flow is 
% more interesting
angle = IA;

X = oc.initConfig(prams.N,'nv',prams.nv,...
  'reducedArea',RA,...
  'angle',angle,...
  'center',[cenx;ceny]);

% Scale the vesicle's radius 
[~, area, ~] = oc.geomProp(X);
rad = sqrt(area/pi);
scale = 1/rad;

X(1:end/2) = scale * (X(1:end/2)-cenx) + cenx;
X(end/2+1:end) = scale * (X(end/2+1:end)-ceny) + ceny;
disp(cenx) 
disp(ceny)
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
end
% Run vesicle code

