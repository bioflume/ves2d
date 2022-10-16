clear; clc;
addpath ../src/
addpath ./output/
addpath ../examples/

saveFile = 'nv10N64Truevelocity';
if 0
fileName = 'nv10N24DNN10VesSavedIC_bgFlowcouette_speed100.mat';
%fileName = 'nv5N24DNNsavedIC_bgFlowcouette_speed100.mat';
load(fileName)
X = Xhist;
else
%fileName = 'nv5N96TrueLoadIC_bgFlowcouette_speed100.mat';   
fileName = 'nv10N6410VesTrueLoadIC_bgFlowcouette_speed100.mat';
load(fileName)
X = XhistTrue;
time = timeTrue;
end


% time steps when the velocity is computed
timeSteps = [2501 5001 7501 10001 15001];


% get the walls x and y coordinates 
wallx = Xwalls(1:end/2,:); wally = Xwalls(end/2+1:end,:);

% tracers
radius = [1.1; 2.1];
Ntheta = 255;
Nradii = 96;
thetapos = [0:Ntheta-1]'/(Ntheta)*2*pi;
thetapos = [thetapos;2*pi];
rpos   = linspace(radius(1),radius(2),Nradii)';
[ttheta,rrpos]      = meshgrid(thetapos,rpos); 
xtracers = rrpos.*cos(ttheta); ytracers = rrpos.*sin(ttheta);          

%xpoints = linspace(-2.2,2.2,200);
%ypoints = linspace(-2.2,2.2,200);
%[xx,yy] = meshgrid(xpoints,ypoints);
%inPnts = find(sqrt(xx.^2+yy.^2)>=1.15 & sqrt(xx.^2+yy.^2) <=2.1);
%xin = xx(inPnts); yin = yy(inPnts);

% Options and Prams
prams.T = time(end);
options.farField = 'couette'; % background velocity
options.farFieldSpeed = 100; % scaling of background velocity
options.confined = true; % confined or unbounded geometry
options.diffDiscWalls = false; % walls are discretized by different Nbd
options.order = 1;
% Descriptive run name
prams.runName = 'couetteTracers'; 
options.usePlot = false; % Plot on-the-fly
options.track = false; % trackers on membrane

prams.nv = numel(X(1,:,1));
prams.N = 96;
prams.Nbd = 256;
prams.nvbd = 2;
prams.NbdInt = 0; prams.NbdExt = 0; prams.nvbdInt = 0; prams.nvbdExt = 0;
options.verbose = 0;
options.saveData = 0;
options.usePlot = 0;
options.track = 0;
options.quiver = 0;
options.axis = 0;

% Vesicle
prams.kappa = 1; % bending coefficient
prams.vesViscCont = 1; % Vesicle viscosity contrast
prams.viscCont = ones(prams.nv,1);

% parameters for numerics
prams.gmresTol = 1e-8;  % tolerance for gmres

options.fmm = false;  % fmm for single-layer potentials
options.fmmDLP = false; % fmm for double-layer potentials
options.matFreeWalls = false; % W2W interactions are done without a matrix
options.antiAlias = true; % upsampling for anti-aliasing
oc = curve;

[options,prams] = initVes2D(options,prams);

om = monitor(X(:,1,1),options,prams);
tt = tstep(options,prams,om);
Xwalls = [interpft(Xwalls(1:end/2,:),256);interpft(Xwalls(end/2+1:end,:),256)];
[walls,~,~] = tt.initialConfined(prams,Xwalls,[],[]);


% -------------------------------------------------------------------------
% PUT TRACERS TOGETHER
%Xtra = [xin(:);yin(:)];
%xtracers = xx; ytracers = yy;
Xtra = [xtracers(:);ytracers(:)];
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;
Ntra = tracers.N;

% Get the near structure
[~,NearW2T] = walls.getZone(tracers,2);

velxTra = zeros(numel(xtracers(:,1)),numel(xtracers(1,:)),numel(timeSteps));
velyTra = velxTra;

for it = 1:numel(timeSteps)
  
  Xves = [interpft(X(1:end/2,:,timeSteps(it)),96);...
    interpft(X(end/2+1:end,:,timeSteps(it)),96)];

  % first compute sigma and eta
  vesicle = capsules(Xves,[],[],prams.kappa,prams.viscCont,1);
  vesicle.setUpRate();
  [sigma,eta,~,~,RS,~,~] = vesicle.computeSigAndEta(...
    tt,walls,[],[]);

  % now compute velocity on the tracers
  vel = tt.tracersVelCouette(vesicle.X,sigma,[],prams.kappa,...
      prams.viscCont,walls,eta,RS,Xtra,NearW2T);
  %velxIn = zeros(numel(xtracers(:,1)),numel(xtracers(:,1)));
  %velyIn = velxIn;
  %velxIn(inPnts) = vel(1:end/2);
  %velyIn(inPnts) = vel(end/2+1:end);
  %velxTra(:,:,it) = velxIn; velyTra(:,:,it) = velyIn;
  velxTra(:,:,it) = reshape(vel(1:end/2),numel(xtracers(:,1)),numel(xtracers(1,:)));
  velyTra(:,:,it) = reshape(vel(end/2+1:end),numel(xtracers(:,1)),...
    numel(xtracers(1,:)));
  
  disp(['Step #' num2str(it) ' is done.'])
end
Xstore = X(:,:,timeSteps);
fileName = ['./output/' saveFile '.mat'];
save(fileName,'xtracers','ytracers','timeSteps','velxTra','velyTra','Xstore','Xwalls','-v7.3')
