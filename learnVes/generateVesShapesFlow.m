% clear; clc;
% runName = ['dataGenRun'];
% folderName = './output/';
% nv = 5;
% speed = 10000;
function generateVesShapesFlow(runName, folderName, nv, speed)
addpath ../src/
addpath ../examples/
rng('shuffle')

prams.folderName = folderName;
fileName = [prams.folderName runName '.bin'];
logFile = [prams.folderName runName '.log'];

prams.farField = 'rotateDataGen'; % 'rotation' or 'couette' (w/ solid boundaries)
prams.speed = speed; 
iplot = 0;
iCalcVel = 0;


% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
prams.Th = 1.5/(prams.speed/100); % time horizon

prams.N = 128; % num. points for true solve in DNN scheme
prams.nv = nv; 
prams.viscCont = ones(prams.nv,1);
prams.fmm = ~false; % use FMM for ves2ves
prams.fmmDLP = ~false; % use FMM for ves2walls
prams.kappa = 1;

prams.dt = 1E-5/(prams.speed/100); % time step size
dtInit = prams.dt;
tsave = 50*dtInit;

prams.outWallRad = 2;
prams.inWallScale = 0.45;
prams.NbdExt = 1024;
prams.NbdInt = 512;
prams.nvbdInt = 4;
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


fid = fopen(logFile, 'w');
fclose(fid);
message = ['Speed is ' num2str(prams.speed)];
writeMessage(logFile,message, '%s\n');
message = ['Viscosity Contrast is ' num2str(prams.viscCont(1))];
writeMessage(logFile,message, '%s\n');
message = ['Bending stiffness is ' num2str(prams.kappa)];
writeMessage(logFile,message, '%s\n');
areaFraction = sum(areaVes)/(sum(areaExt)-sum(areaInt));
message = ['Area fraction is ' num2str(areaFraction)];
writeMessage(logFile,message, '%s\n');

if iCalcVel
% Calculate velocity field
[vel, tracers] = calculateVelocityField(wallsInt, wallsExt, tt, prams);
ng = sqrt(tracers.N);
xx = tracers.X(1:end/2); xx = reshape(xx,[ng ng]);
yy = tracers.X(end/2+1:end); yy = reshape(yy,[ng ng]);
velx = vel(1:end/2); vely = vel(end/2+1:end);

Vx = reshape(velx, [ng ng]);
Vy = reshape(vely, [ng ng]);

figure(1); clf;
plot([XwallsExt(1:end/2); XwallsExt(1)], [XwallsExt(end/2+1:end); XwallsExt(end/2+1)],'k','linewidth',2)
hold on
axis equal
plot([XwallsInt(1:end/2,:); XwallsInt(1,:)], [XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k','linewidth',2)
hWalls = fill([XwallsInt(1:end/2,:); XwallsInt(1,:)],[XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k');
set(hWalls,'edgecolor','k')

quiver(xx, yy, Vx, Vy)
pause
end

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
  [Xcorr,ifail] = oc.correctAreaAndLength2(Xnew,prams.area0,prams.len0);
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
  vesicleProv = capsules(X,[],[],[],[],1);
  vesicleProv.setUpRate(tt.op);
  [NearV2V,NearV2Wint] = vesicleProv.getZone(wallsInt,3);  
  [~,NearV2Wext] = vesicleProv.getZone(wallsExt,2);
  [~,icollisionWallExt] = vesicleProv.collision(wallsExt,...
      NearV2V,NearV2Wext,prams.fmm,tt.op);
  [icollisionVes,icollisionWallInt] = vesicleProv.collision(wallsInt,...
    NearV2V,NearV2Wint,prams.fmm,tt.op);
  icollisionWall = icollisionWallInt || icollisionWallExt;
  
  if icollisionWall
    message = 'Vesicle-wall collision occurs'; 
    writeMessage(logFile,message,'%s\n');
  end
  if icollisionVes
    message = 'Vesicle-vesicle collision occurs'; 
    writeMessage(logFile,message,'%s\n');
  end

  if errAreaLength > 1e-2 
    writeMessage(logFile,message,'%s\n');  
    prams.dt = prams.dt/2;
    message = ['Time step rejected, taking it with a smaller step: ' num2str(prams.dt)];
    writeMessage(logFile,message,'%s\n');
    if prams.dt < 1e-9
      break
    end
  else
    X = Xnew; sig = sigNew; etaExt = etaExtNew; etaInt = etaIntNew; RS = RSNew;  
    time = time + prams.dt;
    prams.dt = min(dtInit, prams.dt*1.2);
    it = it + 1;
    if time >= nextSave 
      writeDataWithEta(fileName,X,sig,etaExt,etaInt,RS,time);
      nextSave = nextSave + tsave;
      message = ['Saving Data'];
      writeMessage(logFile,message,'%s\n');
    end
  end
  tt.dt = prams.dt;

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
end % end function
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
angle = pi * rand(prams.nvbdInt,1);
% inner boundary 
thet = (0:prams.NbdInt-1)'*2*pi/prams.NbdInt;
% sample shapes
% 1: circle, 2: optim-shape1, 3: cShape
% 4: star3, 5: star4, 6: optim-shape2, 7: optim-shape3
% 8: optim-shape4, 9: optim-shape5

ishapes = randi(9,[1,prams.nvbdInt]);
ishapes(1) = 1; % circle for the first one

XwallsInt = zeros(prams.NbdInt*2,prams.nvbdInt);
Xps = zeros(prams.NbdInt*2,prams.nvbdInt);
sizes = zeros(prams.nvbdInt,1);
for ik = 1 : prams.nvbdInt
  if ishapes(ik) == 1 % circle
    Xp = [cos(-thet); sin(-thet)];
  elseif ishapes(ik) == 2 % optim-shape1
    points = [[6.6 3.9]; [15.1 6.3]; [5.6 6.2]; [-10.6 5.6]; [-2.0 4.1]; ...
        [-2.8 -6.5]; [6.0 -13.7]; [1.8 -10.8]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);
  elseif ishapes(ik) == 3 % cShape
    Xp = [-(1.5+sin(-thet)).*(cos(0.85*pi*cos(-thet))); (1.5+sin(-thet)).*(sin(0.85*pi*cos(-thet)))];
  elseif ishapes(ik) == 4 % star3
    folds = 3;
    radius = 1 + 0.3*cos(-folds*thet);
    Xp = [radius.*cos(-thet);radius.*sin(-thet)];
  elseif ishapes(ik) == 5 % star4
    folds = 4;
    radius = 1 + 0.3*cos(-folds*thet);
    Xp = [radius.*cos(-thet);radius.*sin(-thet)];  
  elseif ishapes(ik) == 6 % optim-shape2
    points = [[6.3 2.7]; [11.6 6.0]; [4.3 7.7]; [-12.8 5.1]; [-6.3 4.9]; ...
        [-4.7 -6.1]; [5.3 -13.6]; [1.0 -8.6]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);
  elseif ishapes(ik) == 7 % optim-shape3
    points = [[0.3 4.2]; [13.8 11.6]; [-2.0 10.6]; [-12.0 10.4]; [-6.9 -4.3]; ...
        [-6.3 -0.9]; [2.5 -9.2]; [1.0 -6.1]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);
  elseif ishapes(ik) == 8 % optim-shape4
    points = [[0.6 4.4]; [14.5 10.1]; [-1.8 10.2]; [-12.1 7.9]; [-6.9 1.4]; ...
        [-10.3 -6.9]; [3.6 -7.5]; [2.3 -9.9]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);
  elseif ishapes(ik) == 9 % optim-shape5
    points = [[10.0 9.1]; [0.0 9.1]; [-10.0 9.1]; [-9.1 7.5]; [-0.9 -7.5]; ...
        [0 -9.1]; [0.9 -7.5]; [9.1 7.5]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);  
  end
  % Scale down the arc-length
  Xp = [Xp(1:end/2)-mean(Xp(1:end/2)); Xp(end/2+1:end)-mean(Xp(end/2+1:end))];
  Lx = max(Xp(1:end/2))-min(Xp(1:end/2));
  Ly = max(Xp(end/2+1:end))-min(Xp(end/2+1:end));
  sizes(ik) = max(Lx,Ly);
  Xp = prams.inWallScale * Xp/sizes(ik);
  Xps(1:end/2,ik) = Xp(1:end/2)*cos(angle(ik)) - Xp(end/2+1:end)*sin(angle(ik));
  Xps(end/2+1:end,ik) = Xp(1:end/2)*sin(angle(ik)) + Xp(end/2+1:end)*cos(angle(ik));
  Lx = max(Xps(1:end/2,ik))-min(Xps(1:end/2,ik));
  Ly = max(Xps(end/2+1:end,ik))-min(Xps(end/2+1:end,ik));
  sizes(ik) = max(Lx,Ly);
end

% Place the interior walls
xcs = zeros(prams.nvbdInt,1); ycs = zeros(prams.nvbdInt,1); 
for ik = 1 : prams.nvbdInt
  K = [(1:ik-1) (ik+1:prams.nvbdInt)];
  iplace = false; iter = 0;
  while ~iplace
    
    rad = 0.8*prams.outWallRad * rand();
    if ik == 1; rad = 0; end; % put in the center the first one

    thet = 2*pi*rand();
    xc = rad*cos(thet); yc = rad*sin(thet);

    dx = xc-xcs(K); dy = yc-ycs(K);
    dists = sqrt(dx.^2 + dy.^2)-(sizes(ik)/2+sizes(K)/2);

    % Check with the wall
    
    if sqrt(xc^2 + yc^2) + sizes(ik)/2 >= 0.9*prams.outWallRad
        iplace = false; 
    else
        iplace = true;
    end
    if iplace == true
        if ik>1 && any(dists < 0.2*prams.inWallScale)
            iplace = false; 
        end
    end
    iter = iter + 1;
    disp([num2str(ik) 'th post, iteration: ' num2str(iter)])
  end
  if iplace
    
    disp('Placed.') 
    xcs(ik) = xc; ycs(ik) = yc;
    XwallsInt(:,ik) = [Xps(1:end/2,ik)+xc; Xps(end/2+1:end,ik)+yc];
  end
end



end % initializeWalls
