clear; clc;
restoredefaultpath;

% Background: u = (Umax-Gamma*y^2), Umax: maximum velocity, Gamma: Curvature
% Velocity changes direction at y = sqrt(Umax/Gamma) and is max at y = 0
options.farFieldUmax = 10; % Umax
options.farFieldCurvature = 1/5; % Gamma in the above formula
options.farFieldSpeed = 1; % This scales the velocity (Umax-Gamma*y^2)
prams.nv = 2;    % number of vesicles (no need to change)
options.farField = 'parabolic2'; % background velocity (no need to change) 

% Spatial-resolution
prams.N = 96;    % points per vesicle

% Temporal Resolution 
prams.T = 50;   % time horizon 

% physical properties of vesicles
prams.viscCont = [1; 1];   % viscosity contrast of each vesicle 
prams.kappa = 1e-1;   % bending stiffness (nondim), should be in the range [1e-2, 1], I suggest do not change this

% Save vesicle information and create a log file
% Name of log file for saving messages
options.logFile = 'output/parabolic2Ves.log';
% Name of binary data file for storing vesicle information
options.dataFile = 'output/parabolic2VesData.bin';

options.usePlot = false; % Plot on-the-fly
% Axis for the plot (if options.usePlot = true)
options.axis = [-20 20 -10 30]; 

% Initial configuration of vesicles
angles = [pi/2; pi/2]; % initial angle w.r.t. the x-axis for each vesicle
centerxs = [0 0]; % x-coordinates of the vesicles' centers (ROW VECTOR)
centerys = [-2 2]; % y-coordinates of the vesicles' centers (ROW VECTOR)
ras = [0.6; 0.6]; % reduced areas of the vesicles
effectiveRadii = [1; 1]; % this the effective radius (= sqrt(vesicle_area/pi))

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% NO NEED TO CHANGE ANYTHING BELOW THIS POINT
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
prams.gmresTol = 1e-10;  % tolerance for gmres
prams.errorTol = 1e-1; % maximum error (used for constant time stepping)
% Low-Resolution Correction Algorithms
options.antiAlias = true; % upsampling for anti-aliasing
options.reparameterization = true; % reparametrization of membrane
options.correctShape = true; % correcting area-length errors
options.alignCenterAngle = true; % aligning inc. angle and center after LRCA
options.timeAdap = true; % time adaptivity
prams.dtMax = 1e-3/options.farFieldSpeed;
prams.dtMin = 1e-5/options.farFieldSpeed;
prams.betaInc = 1e-2;
prams.betaDec = 5e-2;
prams.betaUp = 1.5;
prams.betaDown = 0.5;
prams.alpha = 0.9;
prams.m = ceil(prams.T/0.005); 
prams.rtolArea = 1e-3;
prams.rtolLength = 1e-3;

% Set options and parameters that the user doesn't
% Also add src to path
[options,prams] = initVes2D(options,prams);


oc = curve;


X = oc.initConfig(prams.N,'nv',prams.nv,'angle',angles,...
    'center',[centerxs;centerys],'reducedArea',ras);

% Get the length and velocity scale of simulation to decide on params.
[~,area,~] = oc.geomProp(X); 
for k = 1 : prams.nv
  rad = sqrt(area(k)/pi);
  scale = 1/rad*effectiveRadii(k);
  X(1:end/2,k) = scale * (X(1:end/2,k) - centerxs(k)) + centerxs(k);
  X(end/2+1:end,k) = scale * (X(end/2+1:end,k) - centerys(k)) + centerys(k);  
end


% Run vesicle code
Xfinal = Ves2D(X,[],[],[],prams,options,[],[]);


