function [options,prams] = initVes2D(options,prams)
% Set a path pointing to src directory and set options and
% prams to default values if they are left blank

P = path; ii = find(pwd == filesep); ii = ii(end);
subPath = pwd; subPath = [subPath(1:ii) 'src'];
if isempty(strfind(P, subPath)),addpath(subPath);end

PramList = {'N','nv','T','m','Nbd','NbdInt','NbdExt','nvbd','nvbdInt',...
    'nvbdExt','kappa',...
    'viscCont','gmresTol','gmresMaxIter','errorTol','areaLenTol',...
    'betaUp','betaDown','alpha','adStrength','adRange','minDist',...
    'repStrength','dtMax','dtMin','betaInc','betaDec',...
    'maxReparamIter','lev_max','lmax','etolFD','runName',...
    'deltaTtracers','totnv','xrange','yrange','scaleVes','vesViscCont',...
    'streamRate','nrow','Dpost','Dx','Dy','epsilon','ncol','xfreeze',...
    'ves2ViscCont','ves2Concent','folderName','nSeed','Dpostx','Dposty',...
    'yZigZag','freezeIfZigZag','nPutBackMax','vortexSize'};
defaultPram.N = 64;
defaultPram.nv = 1;
defaultPram.Nbd = 0;
defaultPram.NbdInt = 0;
defaultPram.NbdExt = 0;
defaultPram.nvbd = 0;
defaultPram.nvbdInt = 0;
defaultPram.nvbdExt = 0;
defaultPram.T = 1;
defaultPram.m = 100;
defaultPram.kappa = 1e-1;
defaultPram.viscCont = 1;
defaultPram.gmresTol = 1e-12;
defaultPram.gmresMaxIter = 200;
defaultPram.errorTol = 1e-1;
defaultPram.areaLenTol = 1e-2;
defaultPram.betaUp = 1.2;
defaultPram.betaDown = 0.5;
defaultPram.alpha = 0.9;
defaultPram.adRange = 8e-1;
defaultPram.adStrength = 4e-1;
defaultPram.repStrength = 900;
defaultPram.minDist = 0.4;
defaultPram.dtMax = 2;
defaultPram.dtMin = 1e-4;
defaultPram.betaInc = 1e-2;
defaultPram.betaDec = 5e-2;
defaultPram.maxReparamIter = 50;
defaultPram.lev_max = 3;
defaultPram.lmax = 700;
defaultPram.etolFD = 1e-8;
defaultPram.runName = 'defaultRun';
defaultPram.deltaTtracers = 0.01;
defaultPram.totnv = 0;
defaultPram.xrange = [];
defaultPram.yrange = [];
defaultPram.scaleVes = 0;
defaultPram.vesViscCont = 1;
defaultPram.nrow = 0;
defaultPram.Dpost = 0;
defaultPram.Dpostx = 0;
defaultPram.Dposty = 0;
defaultPram.Dy = 0;
defaultPram.Dx = 0;
defaultPram.epsilon = 0;
defaultPram.ncol = 0;
defaultPram.streamRate = 0;
defaultPram.xfreeze = 0;
defaultPram.ves2Concent = 0;
defaultPram.ves2ViscCont = 1;
defaultPram.folderName = './output/streamSim/';
defaultPram.nSeed = 0;
defaultPram.yZigZag = 0;
defaultPram.freezeIfZigZag = 0;
defaultPram.nPutBackMax = [];
defaultPram.vortexSize = 1;

for k = 1:length(PramList)
  if ~isfield(prams,PramList{k})
    eval(['prams.' PramList{k} '=defaultPram.' PramList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end

OptList = {'diffDiscWalls','order','inextens','farField',...
    'farFieldSpeed','vesves','fmm','fmmDLP','confined','usePlot',...
    'track','quiver','axis','saveData','logFile','dataFile','verbose',...
    'profile','collision','timeAdap','bending','tracers',...
    'pressure','orderGL','nsdc','adhesion','repulsion'...
    'periodic','correctShape','antiAlias','reparameterization',...
    'filterShape','fmmPrecision','alignCenterAngle','fastDirect',...
    'wallMatFile','outOfCore','haveWallMats','memsize',...
    'saveWallMat','matFreeWalls','streaming','freezing','putBackOrigin',...
    'putBackDLD','equiDistArcLength','alsoExplicit','usePreco','HODLRforW2W',...
    'randVesicles','saveVinf','saveVtotal'};
defaultOpt.diffDiscWalls = false;
defaultOpt.order = 1;
defaultOpt.inextens = 'method1';
defaultOpt.farField = 'shear';
defaultOpt.farFieldSpeed = 1;
defaultOpt.vesves = 'implicit';
defaultOpt.fmm = false;
defaultOpt.fmmDLP = false;
defaultOpt.confined = false;
defaultOpt.usePlot = true;
defaultOpt.track = false;
defaultOpt.quiver = false;
defaultOpt.axis = [];
defaultOpt.saveData = true;
defaultOpt.logFile = 'output/example.log';
defaultOpt.dataFile = 'output/exampleData.bin';
defaultOpt.verbose = true;
defaultOpt.profile = false;
defaultOpt.collision = false;
defaultOpt.timeAdap = false;
defaultOpt.bending = false;
defaultOpt.tracers = false;
defaultOpt.SDCcorrect = false;
defaultOpt.filterShape = false;
defaultOpt.pressure = false;
defaultOpt.orderGL = 2;
defaultOpt.nsdc = 0;
defaultOpt.adhesion = false;
defaultOpt.repulsion = false;
defaultOpt.periodic = false;
defaultOpt.correctShape = false;
defaultOpt.antiAlias = true;
defaultOpt.reparameterization = false;
defaultOpt.alignCenterAngle = false;
defaultOpt.fmmPrecision = 4;
defaultOpt.fastDirect = false;
defaultOpt.wallMatFile = 'output/LargeMatDLD/wallMat';
defaultOpt.outOfCore = false;
defaultOpt.haveWallMats = false;
defaultOpt.memsize = 1; % in gb
defaultOpt.saveWallMat = false;
defaultOpt.matFreeWalls = true;
defaultOpt.streaming = false;
defaultOpt.freezing = false;
defaultOpt.putBackOrigin = false;
defaultOpt.putBackDLD = false;
defaultOpt.equiDistArcLength = false;
defaultOpt.alsoExplicit = false;
defaultOpt.usePreco = true;
defaultOpt.HODLRforW2W = false;
defaultOpt.randVesicles = [];
defaultOpt.saveVinf = false;
defaultOpt.saveVtotal = false;

for k = 1:length(OptList)
  if ~isfield(options,OptList{k})
    eval(['options.' OptList{k} '=defaultOpt.' OptList{k} ';'])
    % Set any unassigned options to a default value
  end
end

% If the geometry is unbounded, make sure to set the number
% of points and number of components of the solid walls
% to 0.  Otherwise, later components will crash
if ~options.confined
  prams.Nbd = 0;
  prams.Nbdcoarse = 0;
  prams.nvbd = 0;
end

if numel(prams.viscCont) ~=prams.nv
  prams.viscCont = prams.viscCont*ones(1,prams.nv);
end


if options.nsdc > 0
  if options.order > 1
    fprintf('***************************************************\n')
    fprintf('Can only do sdc updates with first-order\n');
    fprintf('Setting sdc corrections to zero\n');
    fprintf('PUSH ANY KEY TO CONTINUE\n');
    pause
    fprintf('***************************************************\n')
    options.nsdc = 0;
  end

  if (any(prams.viscCont ~= 1) && strcmp(options.vesves,'explicit') && options.nsdc > 0)
    fprintf('***************************************************\n')
    fprintf('Not sure if this combination works\n');
    fprintf('See how rhs1 is built at SDC corrections\n');
    fprintf('PUSH ANY KEY TO CONTINUE\n');
    pause
    fprintf('***************************************************\n')
%    options.nsdc = 0;
  end
end




