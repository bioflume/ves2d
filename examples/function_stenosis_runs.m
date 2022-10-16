function function_stenosis_runs(runName,VC,Vinf,RA,kappa,tube_length,tube_height,min_height,contraction_width)

%runName = 'test';
%VC = 1; % viscosity contrast 
%Vinf = 1; % farfield speed
%RA = 0.65; % reduced area
%kappa = 1e-1; % bending stiffness
IA = pi/2; % initial inclination angle

% Tube dimensions
%tube_length = 50;
%tube_height = 6;
%min_height = 3;
%contraction_width = 2;
Nwall = 512;


% Physics parameters
prams.N = 96;                           % points per vesicle
prams.Nbd = Nwall;                     % points per wall
prams.nvbd = 1;                        % number of solid walls
prams.nv = 1;                          % number of vesicles
prams.T = ceil(2000/Vinf);  % time horizon
prams.m = 2e4;                         % number of time steps
prams.kappa = kappa;                   % bending coefficient
prams.errorTol = 1e-1;
prams.viscCont = VC*ones(prams.nv,1); % viscosity contrast
prams.gmresTol = 1e-10;
options.farField = 'choke'; % background velocity
options.confined = true;
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
% use FMM to compute single-layer potential

options.logFile = ['./output/' runName '.log'];
% Name of log file for saving messages
options.dataFile = ['./output/' runName 'Data.bin'];
% Name of binary data file for storing vesicle information

options.profile = false;

% ADD-ONS
options.alignCenterAngle = ~true;
options.correctShape = true;
options.reparameterization = ~true;
prams.maxReparamIter = 5;

options.repulsion = true;
prams.minDist = 0.3; %0.3
prams.minSpecRatio = 90; %30
prams.repStrength = 90; %90

options.timeAdap = true;
prams.rtolArea = 1e-3;
prams.rtolLength = 1e-3;

prams.dtMax = 1e-2/Vinf;
prams.dtMin = 1e-4/Vinf;
prams.betaInc = 1e-2;
prams.betaDec = 5e-2;
prams.betaUp = 1.5;
prams.betaDown = 0.5;
prams.alpha = 0.9;


% Plot on-the-fly
options.usePlot = false;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't

oc = curve;


%% Build the wall
a = tube_length/2; 
b = tube_height/2;
c = tube_height/min_height - 1;
order = 8;
Nsides = ceil(0.5*b/(2*a+2*b)*prams.Nbd);
Ntop = (prams.Nbd-4*Nsides)/2;

t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
t = [t1;t2;t3;t4;t5];

r = (cos(t).^order + sin(t).^order).^(-1/order);
x = a*r.*cos(t); y = b*r.*sin(t);
ind = abs(x)/contraction_width * pi < pi & y < 0;
y(ind) = y(ind).*(1-c*cos(x(ind)/contraction_width * pi))/(1+c);
Xwalls = [x;y];
%%
% ############# VESICLE ######################
% x and y coordinates of centers
% Have added pertubation to centers so that flow is 
% more interesting
angle = IA;

X = oc.initConfig(prams.N,'nv',prams.nv,...
  'reducedArea',RA,...
  'angle',angle,...
  'center',[0;0]);

% Scale the vesicle's radius 
maxL = max(X(end/2+1:end))-min(X(end/2+1:end));

X(1:end/2) = tube_height*0.8/maxL * (X(1:end/2)) ;
X(end/2+1:end) = tube_height*0.8/maxL * (X(end/2+1:end));
maxL = max(X(end/2+1:end))-min(X(end/2+1:end));
maxW = max(X(1:end/2))-min(X(1:end/2));
Xg(1) = -tube_length/2+maxW;
Xg(2) = 0;
X(1:end/2) = X(1:end/2) + Xg(1);
X(end/2+1:end) = X(end/2+1:end) + Xg(2);

[reduced_area,area,length] = oc.geomProp(X);

%%

% Print some information about the run
runInfoFile = ['./output/' runName '_runInfo.log'];
fid = fopen(runInfoFile,'w');
fprintf(fid,'%s\n',['Tube length is ' num2str(tube_length)]);
fprintf(fid,'%s\n',['Tube height is ' num2str(tube_height)]);
fprintf(fid,'%s\n',['Gap size is ' num2str(min_height)]);
fprintf(fid,'%s\n',['Gap width is ' num2str(contraction_width)]);
fprintf(fid,'%s\n',['Poiseuille flow Vinf is ' num2str(Vinf)]);
fprintf(fid,'%s\n',['Height of the vesicle is ' num2str(maxL)]);
fprintf(fid,'%s\n',['Area of the vesicle is ' num2str(area)]);
fprintf(fid,'%s\n',['Length of the vesicle is ' num2str(length)]);
fprintf(fid,'%s\n',['Reduced Area of the vesicle is ' num2str(reduced_area)]);
fprintf(fid,'%s\n',['Bending stiffness of the vesicle is ' num2str(kappa)]);
fprintf(fid,'%s\n',['Viscosity contrast is ' num2str(VC)]);

%% Initial configuration

% PLOT INITIAL CONFIGURATION
% figure(1);clf;
% xvec = [X(1:end/2,:);X(1,:)];
% yvec = [X(end/2+1:end,:);X(end/2+1,:)];
% plot(Xwalls(1:end/2),Xwalls(end/2+1:end),'k','linewidth',2)
% hold on
% plot(xvec,yvec,'r','linewidth',2)
% axis equal
% pause
%%
Ves2D_stenosis(X,Xwalls,[],[],prams,options,[]);

% Run vesicle code

