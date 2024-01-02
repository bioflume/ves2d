clear; clc;
% rng('shuffle')

runName = ['interpTestRun'];
folderName = './output/';
nv = 10;
speed = 1000;

addpath ../src/
addpath ../examples/

prams.folderName = folderName;
fileName = [prams.folderName runName '.bin'];
logFile = [prams.folderName runName '.log'];

prams.farField = 'couetteTwoSep'; % 'rotation' or 'couette' (w/ solid boundaries)
prams.speed = speed; 
iplot = 1;
iCalcVel = 0;


% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
prams.Th = 1.5/(prams.speed/100); % time horizon

prams.N = 16; % num. points for true solve in DNN scheme
prams.nv = nv; 
prams.viscCont = ones(prams.nv,1);
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.nLayers = 5;
prams.dt = 1E-3/(prams.speed/100); % time step size 1E-5
dtInit = prams.dt;
tsave = 5*dtInit;

prams.repStrength = 1E+5;
prams.outWallRad = 2;
prams.inWallRad = 1;
prams.NbdExt = 256;
prams.NbdInt = 128;
prams.nvbdInt = 1;
prams.nvbdExt = 1;
prams.nvbd = prams.nvbdInt + prams.nvbdExt;
oc = curve;


% VESICLES and WALLS:
% -------------------------------------------------------------------------
vesShape = 'ellipse'; % 'ellipse' or 'curly' 
reducedArea = 0.65;
[XwallsInt,XwallsExt] = initializeWalls(prams,oc);

% Generate a sample vesicle
X = oc.ellipse(prams.N,reducedArea);
[ra,area,length] = oc.geomProp(X);
X = X/length;
Xref = X;

% Build tt to initialize vesicles
tt = buildTstep(X,prams);
[~,wallsInt,wallsExt] = tt.initialConfined(prams,[],XwallsInt,XwallsExt);

Uwall = wallsExt.u;

% Randomly place vesicles
[X,XwallsInt,XwallsExt,area0,len0] = initializeVesicles(X,...
    XwallsExt, XwallsInt, prams,oc, tt);
prams.area0 = area0; prams.len0 = len0;

[~,areaExt,~] = oc.geomProp(XwallsExt);
[~,areaInt,~] = oc.geomProp(XwallsInt);
[~,areaVes,~] = oc.geomProp(X);

if 0
figure(1); clf;
plot([XwallsExt(1:end/2); XwallsExt(1)], [XwallsExt(end/2+1:end); XwallsExt(end/2+1)],'k','linewidth',2)
hold on
axis equal
plot([XwallsInt(1:end/2,:); XwallsInt(1,:)], [XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k','linewidth',2)
hWalls = fill([XwallsInt(1:end/2,:); XwallsInt(1,:)],[XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k');
set(hWalls,'edgecolor','k')
plot([X(1:end/2,:); X(1,:)], [X(end/2+1:end,:);X(end/2+1,:)], 'r', 'linewidth',2)
hVes = fill([X(1:end/2,:); X(1,:)],[X(end/2+1:end,:);X(end/2+1,:)],'r');
set(hVes,'edgecolor','r')
title('Initial Condition')
pause
end

% Initialize matrices and counters
% -------------------------------------------------------------------------
tt.dt = prams.dt;
sig = zeros(prams.N,prams.nv); 
etaInt = zeros(2*prams.NbdInt,prams.nvbdInt); 
etaExt = zeros(2*prams.NbdExt,1); RS = zeros(3,prams.nvbdInt+1);

% initialize files
% -------------------------------------------------------------------------
fid = fopen(fileName, 'w');
fwrite(fid, [prams.dt; prams.nv; prams.N; prams.nvbdExt; prams.NbdExt; prams.nvbdInt; prams.NbdInt],'double');
fwrite(fid,XwallsExt(:),'double');
fwrite(fid,XwallsInt(:),'double')
fclose(fid);


it = 1; time = 0; nextSave = tsave;
while time < prams.Th 
  message = '********************************************';
  writeMessage(logFile,message,'%s\n')
  message = ['Time = '  num2str(time) ' out of ' num2str(prams.Th)];
  writeMessage(logFile,message,'%s\n')

  vesicle = capsules(X,[],[],prams.kappa,prams.viscCont,1);
  vesicle.setUpRate();
  [Xnew,sigNew,etaIntNew,etaExtNew,RSNew,iter,iflag] = tt.timeStepSimpleDiffDisc(X,sig,etaInt,etaExt,...
      RS,prams.viscCont,wallsInt,wallsExt,vesicle);
  message = ['GMRES took ' num2str(iter) ' iterations.'];
  writeMessage(logFile,message,'%s\n')


   % AREA-LENGTH CORRECTION
  message = 'Correcting area-length errors...';
  writeMessage(logFile,message,'%s\n')
  [Xcorr,ifail, failedIdcs] = oc.correctAreaAndLength2(Xnew,prams.area0,prams.len0);
  if ifail
    message = 'AREA-LENGTH CANNOT BE CORRECTED!!!';
    writeMessage(logFile,message,'%s\n')
  end
  Xnew = oc.alignCenterAngle(Xnew,Xcorr); 

  % Equally distribute points in arc-length
  Xiter = Xnew;
  for iter = 1 : 5
%     [Xiter,~,~] = oc.redistributeArcLength(Xiter);
    [Xiter, ~] = oc.reparametrize(Xiter,[],6,20);
  end
  
  % Fix misalignment in center and angle due to reparametrization
  Xnew = oc.alignCenterAngle(Xnew,Xiter);
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(Xnew);
  errArea = max(abs(area-prams.area0)./prams.area0); 
  errLen = max(abs(len-prams.len0)./prams.len0);
  errAreaLength = max(errArea,errLen);
  
  message = ['Max. error in area and length is ' num2str(errAreaLength)];
  writeMessage(logFile,message,'%s\n');

  % Collision check
  vesicleProv = capsules(Xnew,[],[],[],[],1);
  vesicleProv.setUpRate(tt.op);
  [NearV2V,NearV2Wint] = vesicleProv.getZone(wallsInt,3);  
  [~,NearV2Wext] = vesicleProv.getZone(wallsExt,2);
  
  [icollisionVes,collidingVesVesIdcs,~] = vesicleProv.collisionWVesOutput(NearV2V,prams.fmm,tt.op);
  [icollisionWallInt,collidingVesWallIntIdcs,~] = vesicleProv.collisionWwallOutput(wallsInt,NearV2Wint,prams.fmm,tt.op);
  [icollisionWallExt,collidingVesWallExtIdcs,~] = vesicleProv.collisionWwallOutput(wallsExt,NearV2Wext,prams.fmm,tt.op);

  problemVesIDs = [failedIdcs; collidingVesVesIdcs; collidingVesWallIntIdcs; collidingVesWallExtIdcs];
  problemVesIDs = unique(problemVesIDs);
  
  
  if icollisionWallInt
    message = 'Vesicle-wall interior collision occurs'; 
    writeMessage(logFile,message,'%s\n');
  end
  if icollisionWallExt
    message = 'Vesicle-wall exterior collision occurs'; 
    writeMessage(logFile,message,'%s\n');
  end
  if icollisionVes
    message = 'Vesicle-vesicle collision occurs'; 
    writeMessage(logFile,message,'%s\n');
  end
  
  
  if ~isempty(problemVesIDs)
    message = ['There are problematic vesicles to be replaced by randomly placed vesicles'];
    writeMessage(logFile,message,'%s\n');
    % First remove the problematic vesicles
    Xreplaced = zeros(2*prams.N,prams.nv-numel(problemVesIDs));
    XLarge = Xreplaced;
    count = 1;
    for k = 1 : prams.nv
      if all(k ~= problemVesIDs)
        Xreplaced(:,count) = Xnew(:,k);
        XLarge(1:end/2,count) = 1.1*(Xnew(1:end/2,k)-mean(X(1:end/2,k))) + mean(X(1:end/2,k));
        XLarge(end/2+1:end,count) = 1.1*(Xnew(end/2+1:end,k)-mean(X(end/2+1:end,k))) + mean(X(end/2+1:end,k));
        count = count + 1;
      end
    end

    % Now initialize new vesicles to maintain the volume fraction
    Xnew = replaceVesicles(Xref,XLarge,Xreplaced,numel(problemVesIDs),XwallsExt, XwallsInt, prams, tt);
  end
  
  X = Xnew; sig = sigNew; etaExt = etaExtNew; etaInt = etaIntNew; RS = RSNew;  
  time = time + prams.dt;
  it = it + 1;
  if time >= nextSave 
    writeDataWithEta(fileName,X,sig,etaExt,etaInt,RS,time);
    nextSave = nextSave + tsave;
    message = ['Saving Data'];
    writeMessage(logFile,message,'%s\n');
  end
  
  

  if iplot
    figure(1); clf;
    plot([XwallsExt(1:end/2); XwallsExt(1)], [XwallsExt(end/2+1:end); XwallsExt(end/2+1)],'k','linewidth',2)
    hold on
    axis equal
    plot([XwallsInt(1:end/2,:); XwallsInt(1,:)], [XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k','linewidth',2)
    hWalls = fill([XwallsInt(1:end/2,:); XwallsInt(1,:)],[XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k');
    set(hWalls,'edgecolor','k')
    plot([X(1:end/2,:); X(1,:)], [X(end/2+1:end,:);X(end/2+1,:)], 'r', 'linewidth',2)
    hVes = fill([X(1:end/2,:); X(1,:)],[X(end/2+1:end,:);X(end/2+1,:)],'r');
    set(hVes,'edgecolor','r')
    title(['Time = ' num2str(time)])
    pause(0.1)
  end


end % end for it = 1 : ntime
% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeDataWithEta(fileName,X,sigma,etaExt,etaInt,RS,time)
% writeDatawithEta(X,sigma,eta,ea,el,time,res) writes the position,
% tension, density function on solid walls, errors, and time to a
% binary file.  Matlab can later read this file to postprocess the
% data.  Note that sigma, eta, and RS can be post processed given X
% using computeSigAndEta.  This gives a more accurate velocity field
% for studing problems like an advection-diffusion solver in a couette
% apparatus (see Gokberk's work)
 
output = [time;X(:);sigma(:);etaExt(:);etaInt(:);RS(:)];
% format that postProcess/loadfile2.m reads the output
fid = fopen(fileName,'a');
fwrite(fid,output,'double');
fclose(fid);

end % writeData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeMessage(logFile,message,format)
% function writeMessage(message,format) appends message to o.fileName
% with format


fid = fopen(logFile,'a');
fprintf(fid,format,message);
fclose(fid);

disp(message)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vel, tracers] = calculateVelocityField(wallsInt, wallsExt, tt, prams)

[xx, yy] = meshgrid(linspace(-prams.outWallRad, prams.outWallRad,20)', linspace(-prams.outWallRad, prams.outWallRad,20)');

Xtra = [xx(:); yy(:)];

tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;

[~,NearWint2T] = wallsInt.getZone(tracers,2);
[~,NearWext2T] = wallsExt.getZone(tracers,2);

rhs3Ext = wallsExt.u;
rhs3Int = wallsInt.u;
rhs = [rhs3Ext(:); rhs3Int(:); zeros(3*(prams.nvbd-1),1)];

Xn = tt.bdiagWall * rhs;
etaInt = zeros(2*prams.NbdInt,prams.nvbdInt);
RS = zeros(3,prams.nvbd); % rotlets and stokeslets

% unstack the density function
etaExt = Xn(1:2*prams.NbdExt);  % assuming nvbdExt = 1
for k = 1:prams.nvbdInt
  etaInt(:,k) = Xn(2*prams.NbdExt+(k-1)*2*prams.NbdInt+1:2*prams.NbdExt+2*k*prams.NbdInt);
end  

% unstack the rotlets and stokeslets
otlets = Xn(2*prams.NbdExt+2*prams.nvbdInt*prams.NbdInt+1:end);  
for k = 2:prams.nvbd
  istart = (k-2)*3+1;
  iend = 3*(k-1);
  RS(:,k) = otlets(istart:iend);
end

% Now calculate the velocity at the tracer points
potWallInt = tt.opWallInt;
potWallExt = tt.opWallExt;

kernelInt = @potWallInt.exactStokesDL;
kernelExt = @potWallExt.exactStokesDL;
kernelDirectInt = @potWallInt.exactStokesDL;
kernelDirectExt = @potWallExt.exactStokesDL;

jump = -1/2;
DLP = tt.wallDLPext + jump*eye(2*prams.NbdExt);
DLPfun = @(X) potWallExt.exactStokesDLdiag(wallsExt,DLP,X);
FwallExt2Tra = potWallExt.nearSingInt(wallsExt,etaExt,DLPfun,[],...
    NearWext2T,kernelExt,kernelDirectExt,tracers,false,false);

DLP = tt.wallDLPint;
for k = 1:prams.nvbd-1
  DLP(:,:,k) = DLP(:,:,k) + jump*eye(2*prams.NbdInt);
end
DLPfun = @(X) potWallInt.exactStokesDLdiag(wallsInt,DLP,X);
FwallInt2Tra = potWallInt.nearSingInt(wallsInt,etaInt,DLPfun,[],...
    NearWint2T,kernelInt,kernelDirectInt,tracers,false,false);

% compute the velocity on the tracers due to the
% solid walls
Fwall2Tra = FwallExt2Tra + FwallInt2Tra;

% velocity due to the stokeslets and rotlets    
FLets2Tra = zeros(2*tracers.N,1);
for k = 2:prams.nvbd
  stokeslet = RS(1:2,k);
  rotlet = RS(3,k);
  FLets2Tra = FLets2Tra + tt.RSlets(Xtra,wallsInt.center(:,k-1),...
    stokeslet,rotlet);
end  

vel = Fwall2Tra + FLets2Tra;

% Remove points inside pillars
InOut = wallsInt.sortPts(Xtra, 0, NearWint2T, potWallInt);

velx = vel(1:end/2); vely = vel(end/2+1:end);

velx(InOut>0) = NaN; vely(InOut>0) = NaN; 

% Remove points outside the domain
drs = sqrt(Xtra(1:end/2).^2 + Xtra(end/2+1:end).^2);
velx(drs>prams.outWallRad) = NaN;
vely(drs>prams.outWallRad) = NaN;

vel = [velx;vely];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,XwallsInt,XwallsExt,area0,len0] = initializeVesicles(Xref,...
    XwallsExt, XwallsInt, prams,oc, tt)

% Wall structures
wallsExt = capsules(XwallsExt,[],[],0,0,true);
% wallsInt = capsules(XwallsInt,[],[],0,0,true);
XwallsIntLarge = zeros(size(XwallsInt));
for k = 1 : numel(XwallsInt(1,:))
  XwallsIntLarge(1:end/2,k) =1.1*(XwallsInt(1:end/2,k)-mean(XwallsInt(1:end/2,k))) + mean(XwallsInt(1:end/2,k));
  XwallsIntLarge(end/2+1:end,k) =1.1*(XwallsInt(end/2+1:end,k)-mean(XwallsInt(end/2+1:end,k))) + mean(XwallsInt(end/2+1:end,k));
end
wallsInt = capsules(XwallsIntLarge,[],[],0,0,true);

% Coordinates of the reference vesicle
x0 = Xref(1:end/2); y0 = Xref(end/2+1:end);

X = zeros(2*prams.N, prams.nv);
XLarge = zeros(2*prams.N, prams.nv);

k = 1;

while k <= prams.nv
  cx = -0.9*prams.outWallRad + 1.8*prams.outWallRad*rand;
  cy = -0.9*prams.outWallRad + 1.8*prams.outWallRad*rand;

  phi = -pi/4+pi/2*rand;

  xpot = cx + x0*cos(phi) + y0*sin(phi);
  ypot = cy - x0*sin(phi) + y0*cos(phi);

  xpotLarge = cx + 1.1*(x0*cos(phi) + y0*sin(phi));
  ypotLarge = cy - 1.1*(x0*sin(phi) - y0*cos(phi));

  accept = true; % tentatively accept the vesicle

  % create capsule with the potential vesicle
  vesicle = capsules([xpotLarge;ypotLarge],[],[],[],[],true);

  [~,NearV2WInt] = vesicle.getZone(wallsInt,3);
  [~,NearV2WExt] = vesicle.getZone(wallsExt,3);


  [~,icollisionWallExt] = vesicle.collision(wallsExt,[],NearV2WExt,...
    tt.fmm,tt.op);
  [~,icollisionWallInt] = vesicle.collision(wallsInt,[],NearV2WInt,...
    tt.fmm,tt.op);
  if icollisionWallExt || icollisionWallInt
    accept = false;

    message = ['Vesicle crossed wall.'];
    disp(message)
    % at least one of the vesicles's points is outside of one
    % of the solid wall components
  end

  if sqrt(cx.^2 + cy.^2) > prams.outWallRad; accept = false; end;


  if accept 
    % if vesicle is not outside of physical walls, accept it as
    % a potential new vesicle.  It will be kept as long as it intersects
    % no other vesicles  
    X(:,k) = [xpot;ypot];
    XLarge(:,k) = [xpotLarge;ypotLarge];

    % create an object with the current configuration of vesicles
    vesicle = capsules(XLarge(:,1:k),[],[],0,0,true);  


    % see if vesicles have crossed using collision detection code
    % Can be used with or without the fmm    
    NearV2V = vesicle.getZone([],1);
    icollisionVes = vesicle.collision(...
        [],NearV2V,[],tt.fmm,tt.op);

    if icollisionVes
      X(:,k) = 0;
      XLarge(:,k) = 0;
      message = ['Vesicles crossed.'];
      disp(message)
      % if they've crossed, reject it
    else
      k = k + 1;
      message = [num2str(prams.nv-k+1,'%d') ' vesicles left to fill the domain'];
      disp(message)
      % if they haven't crossed, increase the total number of vesicles
    end
  end % if accept

end % end while
[~,area0,len0] = oc.geomProp(X);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = replaceVesicles(Xref,XLarge,X,num2create,...
    XwallsExt, XwallsInt, prams, tt)

% Wall structures
wallsExt = capsules(XwallsExt,[],[],0,0,true);
% wallsInt = capsules(XwallsInt,[],[],0,0,true);
XwallsIntLarge = zeros(size(XwallsInt));
for k = 1 : numel(XwallsInt(1,:))
  XwallsIntLarge(1:end/2,k) =1.1*(XwallsInt(1:end/2,k)-mean(XwallsInt(1:end/2,k))) + mean(XwallsInt(1:end/2,k));
  XwallsIntLarge(end/2+1:end,k) =1.1*(XwallsInt(end/2+1:end,k)-mean(XwallsInt(end/2+1:end,k))) + mean(XwallsInt(end/2+1:end,k));
end
wallsInt = capsules(XwallsIntLarge,[],[],0,0,true);

% Coordinates of the reference vesicle
x0 = Xref(1:end/2); y0 = Xref(end/2+1:end);


k = 1;
while k <= num2create
  cx = -0.9*prams.outWallRad + 1.8*prams.outWallRad*rand;
  cy = -0.9*prams.outWallRad + 1.8*prams.outWallRad*rand;

  phi = -pi/4+pi/2*rand;

  xpot = cx + x0*cos(phi) + y0*sin(phi);
  ypot = cy - x0*sin(phi) + y0*cos(phi);

  xpotLarge = cx + 1.1*(x0*cos(phi) + y0*sin(phi));
  ypotLarge = cy - 1.1*(x0*sin(phi) - y0*cos(phi));

  accept = true; % tentatively accept the vesicle

  % create capsule with the potential vesicle
  vesicle = capsules([xpotLarge;ypotLarge],[],[],[],[],true);

  [~,NearV2WInt] = vesicle.getZone(wallsInt,3);
  [~,NearV2WExt] = vesicle.getZone(wallsExt,3);


  [~,icollisionWallExt] = vesicle.collision(wallsExt,[],NearV2WExt,...
    tt.fmm,tt.op);
  [~,icollisionWallInt] = vesicle.collision(wallsInt,[],NearV2WInt,...
    tt.fmm,tt.op);
  if icollisionWallExt || icollisionWallInt
    accept = false;

    message = ['Vesicle crossed wall.'];
    disp(message)
    % at least one of the vesicles's points is outside of one
    % of the solid wall components
  end

  if sqrt(cx.^2 + cy.^2) > prams.outWallRad; accept = false; end;


  if accept 
    % if vesicle is not outside of physical walls, accept it as
    % a potential new vesicle.  It will be kept as long as it intersects
    % no other vesicles  
    XPotLarge = [XLarge [xpotLarge;ypotLarge]];

    % create an object with the current configuration of vesicles
    vesicle = capsules(XPotLarge,[],[],0,0,true);  


    % see if vesicles have crossed using collision detection code
    % Can be used with or without the fmm    
    NearV2V = vesicle.getZone([],1);
    icollisionVes = vesicle.collision(...
        [],NearV2V,[],tt.fmm,tt.op);

    if icollisionVes
      message = ['Vesicles crossed.'];
      disp(message)
      % if they've crossed, reject it
    else
      X = [X [xpot;ypot]];
      k = k + 1;
      message = [num2str(num2create-k+1,'%d') ' vesicles left to fill the domain'];
      disp(message)
      % if they haven't crossed, increase the total number of vesicles
    end
  end % if accept

end % end while

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tt = buildTstep(X,pramsIn)
N = pramsIn.N; nv = pramsIn.nv;
prams.Nbd = 0;
prams.N = N;
prams.nv = nv;
prams.NbdInt = pramsIn.NbdInt;
prams.NbdExt = pramsIn.NbdExt;
prams.nvbd = pramsIn.nvbdInt+pramsIn.nvbdExt;
prams.nvbdInt = pramsIn.nvbdInt; prams.nvbdExt = pramsIn.nvbdExt;
prams.kappa = pramsIn.kappa;


options.diffDiscWalls = 1;
options.verbose = 0;
options.saveData = 0;
options.usePlot = 0;
options.track = 0;
options.quiver = 0;
options.axis = 0;
prams.T = pramsIn.Th;
prams.gmresTol = 1e-10;
prams.m = pramsIn.Th/pramsIn.dt;
prams.errTol = 1e-1;
options.tracers = 0;
options.timeAdap = 0;
options.order = 1;
prams.areaLenTol = 1e-2;
options.fmm = pramsIn.fmm;
options.fmmDLP = pramsIn.fmmDLP;
options.antiAlias = 1;
options.correctShape = 1;
options.adhesion = 0;
options.repulsion = 1;

options.confined = true;
options.farField = pramsIn.farField;
options.farFieldSpeed = pramsIn.speed;

[options,prams] = initVes2D(options,prams);
    
om = monitor(X,options,prams);
tt = tstep(options,prams,om);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XwallsInt,XwallsExt] = initializeWalls(prams,oc)

% outer boundary
thet = (0:prams.NbdExt-1)'*2*pi/prams.NbdExt;
XwallsExt = [prams.outWallRad * cos(thet); prams.outWallRad*sin(thet)];

% inner boundary 
thet = (0:prams.NbdInt-1)'*2*pi/prams.NbdInt;

XwallsInt = [prams.inWallRad*cos(-thet); prams.inWallRad*sin(-thet)];;

end % initializeWalls
