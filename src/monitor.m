classdef monitor
% Used for doing input/output of data.
% Can do plotting, write to files, compute error in area
% and length, and write to console

properties
areaAll         % keeps areas of all vesicles (already streamed, to be streamed)
lengthAll       % keeps lengths of all vesicles (already streamed, to be streamed)
area            % area of initial vesicles 
length          % length of initial vesicles
reducedArea     % reduced area
verbose         % write data to console
saveData        % save the vesicle positions, tension, etc
usePlot         % plot the vesicles
track           % place tracker points on the vesicles
quiver          % plot the velocity field on the vesicles
axis            % axis of the plot
dataFile        % name of the file to write the data
logFile         % name of the file to write the log
T               % time horizion
m               % number of time steps
errorTol        % error tolerance for errors in area and length
tracers         % use tracers to monitor passive particles
timeAdap        % use time adaptivity
reparameterization % use of reparameterization
repulsion       % use of repulsion
N               % points per vesicle
nv              % number of vesicles
Nbd             % points per solid wall (if any)
nvbd            % number of solid wall components (if any)
NbdInt          % points per interior solid wall
NbdExt          % points on outer solid wall
nvbdInt         % number of interior walls
nvbdExt         % number of exterior walls
order           % time stepping order
orderGL         % order of Gauss-Lobatto quadrature
pressure        % calculate pressure
areaLenTol      % allowable relative area and length change
adhesion        % flag to tell if adhesion is in model
adStrength      % strength of adhesion (W_0)
adRange         % range of adhesion (d_0)
confined        % bounded flow
diffDiscWalls   % whether we have two sets of walls with different disc.
farField        % farField velocity
fmm             % flag for FMM
correctShape    % whether error in area-length is corrected
istreaming      % whether vesicles are streamed or not
end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = monitor(X,options,prams)
% monitor(X,Xwalls,options,prams) saves options and parameters needed
% by the class and also computes the initial error in area and length
% so that the errors in area and length can be monitored throughout the
% simulation.
% This is the constructor

o.N = prams.N;                      % points per vesicle
o.nv = prams.nv;                    % number of vesicles
o.Nbd = prams.Nbd;                  % points per solid wall
o.nvbd = prams.nvbd;                % number of solid walls
o.NbdInt = prams.NbdInt;            % points per inner solid wall
o.NbdExt = prams.NbdExt;            % points on outer solid wall
o.nvbdInt = prams.nvbdInt;          % number of inner solid walls
o.nvbdExt = prams.nvbdExt;          % number of exterior walls
o.verbose = options.verbose;        % write data to console
o.saveData = options.saveData;      % save the data
o.usePlot = options.usePlot;        % plot the data
o.track = options.track;            % include tracker points
o.quiver = options.quiver;          % include velocity vectors
o.axis = options.axis;              % axis of plot
o.dataFile = options.dataFile;      % data file name
o.logFile = options.logFile;        % log file name
o.T = prams.T;                      % time horizon
o.m = prams.m;                      % number of time steps
o.errorTol = prams.errorTol;        % error tolerance
o.tracers = options.tracers;        % tracers
o.timeAdap = options.timeAdap;      % time adpativity
o.order = options.order;            % time stepping order
o.orderGL = options.orderGL;        % order of Gauss-Lobatto 
                                    % quadtrature
o.pressure = options.pressure;      % calculate pressure
o.areaLenTol = prams.areaLenTol;    % allowable error in area and length

o.reparameterization = options.reparameterization;
o.repulsion = options.repulsion;
o.fmm = options.fmm;
o.correctShape = options.correctShape;

o.adhesion = options.adhesion;
o.adStrength = prams.adStrength;
o.adRange = prams.adRange;

o.confined = options.confined;
o.diffDiscWalls = options.diffDiscWalls;
o.farField = options.farField;

if ~isempty(X)
  oc = curve;
  % area, length, and reduced area of initial shape
  [o.reducedArea,o.area,o.length] = oc.geomProp(X);
  o.areaAll = o.area; o.lengthAll = o.length;
end

o.istreaming = options.streaming;

end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initializeFiles(o,X,sig,Xwalls,XwallsInt,XwallsExt,Xtra,pressTar)
% initializeFiles(o,X,sig,eta,etaInt,etaExt,RS,Xwalls,...
%         XwallsInt,XwallsExt,Xtra,pressTar) does the initial
% writing of data to files and the console.  It first deletes any
% previous data and then writes the number of points, tracer initial
% conditions, pressure targets X and Xwalls are the vesicles and solid
% walls, Xtra and pressTar are initial conditions for the tracers and
% pressure/stress target locations

N = size(X,1)/2; % points per vesicle
nv = size(X,2); % number of vesicles

NbdInt = o.NbdInt; % points per inner wall
NbdExt = o.NbdExt; % points on exterior wall
nvbdInt = o.nvbdInt; % number of inner walls
nvbdExt = o.nvbdExt; % number of exterior walls

Nbd = o.Nbd; % points per solid wall
nvbd = o.nvbd; % number of solid walls
oc = curve;

% seperate x and y coordinates of solid walls
[xx,yy] = oc.getXY(Xwalls);
[xxInt,yyInt] = oc.getXY(XwallsInt);
[xxExt,yyExt] = oc.getXY(XwallsExt);


[xtra,ytra] = oc.getXY(Xtra);
% seperate x and y coordinates of tracers
if o.saveData
  fid = fopen(o.dataFile,'w');
  % write number of points and vesicles to data file
  fwrite(fid,[N;nv;Nbd;nvbd;NbdExt;NbdInt;nvbdExt;nvbdInt],'double');
  
  % write the solid walls to the data file
  if o.confined
    if~o.diffDiscWalls
      fwrite(fid,[xx(:),yy(:)],'double');
    else
      fwrite(fid,[xxExt(:);yyExt(:)],'double');
      fwrite(fid,[xxInt(:);yyInt(:)],'double');
    end
  
    fclose(fid);
  end
end

fid = fopen(o.logFile,'w');
fclose(fid);

o.writeMessage(' ','%s\n')
message = ['************* PHYSICAL PARAMETERS *************'];
o.writeMessage(message,'%s\n')

% write number of points
message = [num2str(N) ' points per vesicle'];
o.writeMessage(message,'%s\n')

% write number of vesicles
message = [num2str(nv) ' total vesicles'];
o.writeMessage(message,'%s\n')

message = ['Far-field velocity is ' o.farField];
o.writeMessage(message,'%s\n')

% write time stepping information
message = ['Time horizon is ' num2str(o.T)];
o.writeMessage(message,'%s\n')

if ~o.timeAdap
  message = ['Time step size is ' num2str(o.T/o.m)];
else
  message = ['Area-length tolerance is ' num2str(o.areaLenTol)];
end
o.writeMessage(message,'%s\n')

if o.confined
  if o.diffDiscWalls
    message = ['Number of points per inner wall is ' num2str(NbdInt,'%d')];
    o.writeMessage(message,'%s\n')
    message = ['Number of points per outer wall is  ' num2str(NbdExt,'%d')];
    o.writeMessage(message,'%s\n')
    message = ['Number of inner walls is ' num2str(nvbdInt,'%d')];
    o.writeMessage(message,'%s\n')
  else
    message = ['Number of points per wall is  ' num2str(Nbd,'%d')];
    o.writeMessage(message,'%s\n')
    message = ['Number of walls is ' num2str(nvbd,'%d')];
    o.writeMessage(message,'%s\n')  
  end
end

message = ['FMM is on? ' num2str(o.fmm)];
o.writeMessage(message,'%s\n')
message = ['Repulsion is on? ' num2str(o.repulsion)];
o.writeMessage(message,'%s\n')
message = ['Area-length correction is on? ' num2str(o.correctShape)];
o.writeMessage(message,'%s\n')


if o.saveData
  if o.timeAdap
    fileName = [o.logFile(1:end-4) 'Res.dat'];
    fid = fopen(fileName,'w');
    fclose(fid);
    % initiate the residual file
  end

  if o.tracers
    fileName = [o.logFile(1:end-4) 'Tracers.bin'];
    fid = fopen(fileName,'w');
    fwrite(fid,[xtra;ytra],'double');
    fclose(fid);
  end
  % write tracer intial conditions

  if o.pressure
    fileName = [o.logFile(1:end-4) 'Pressure.bin'];
    fid = fopen(fileName,'w');
    fclose(fid);
    o.writePressure(pressTar.X);
    % write the pressure's target locations
  end

  if o.pressure
    fileName = [o.logFile(1:end-4) 'Stress11.bin'];
    fid = fopen(fileName,'w');
    fclose(fid);
    fileName = [o.logFile(1:end-4) 'Stress12.bin'];
    fid = fopen(fileName,'w');
    fclose(fid);
    fileName = [o.logFile(1:end-4) 'Stress21.bin'];
    fid = fopen(fileName,'w');
    fclose(fid);
    fileName = [o.logFile(1:end-4) 'Stress22.bin'];
    fid = fopen(fileName,'w');
    fclose(fid);
    o.writeStress(pressTar.X,'pts');
    % write the stress's target locations
  end
  % Erase anything left over in log and binary files
  % Write the parameters to log file.  Also write the target
  % locations of the tracers, pressure, and stress if these
  % are being calculated

  % save initial configuartion
  o.writeData(X,sig,0,0);
  
 
  message = ['Initial Areas are:            ' ...
      num2str(o.area(1),'%10.2e')];
  o.writeMessage(message,'%s\n')
  message = ['Initial Lengths are:          ' ...
      num2str(o.length(1),'%10.2e')];
  o.writeMessage(message,'%s\n')
  message = ['Initial Reduced Areas are:    ' ...
      num2str(o.reducedArea(1),'%10.2e')];
  o.writeMessage(message,'%s\n\n')
end
% write initial reduced area, area, and length to log file

o.writeStars
o.writeMessage(' ','%s\n')

end % initializeFiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function terminate = outputInfo(o,X,sigma,u,Xwalls,...
    XwallsInt,XwallsExt,Xtra,time,iter,dtScale,res,rem_ids,dt,allVesViscConts,iflag)
% terminate = outputInfo(o,X,sigma,u,Xwalls,...
%     XwallsInt,XwallsExt,Xtra,time,iter,dtScale,res,rem_ids,iflag)
% computes the error in area and length and write messages to the data
% file, the log file, and the console.  Tells the simulation to stop if
% error in area or length is too large.  X, sigma, u are the position,
% tension, and velocity of vesicles, Xwalls is the parameterization of
% the solid walls, time is the current time, iter is the number of
% GMRES iterations, dtScale is the amount the time step is scaled by,
% res is the residual incurred from the previous time step to the
% current, and iflag is the success of gmres

% check if error in area or length is greater than threshold
errorTol = o.errorTol;
terminate = false;
if ~isempty(rem_ids)
  [ea,el] = o.errors(X(:,rem_ids));
  if max(ea,el) > errorTol
    message = 'ERROR IS TOO LARGE!  I AM STOPPING!';
    o.writeMessage(message,'%s\n');
    message = ['Max error in area is   ',num2str(ea,'%4.2e')];
    o.writeMessage(message,'%s\n');
    message = ['Max error in length is ',num2str(el,'%4.2e')];
    o.writeMessage(message,'%s\n');
    terminate = true;
    return
  end
   % if time step shrinks too much, then terminate the run
  if dt < 1e-6
    message = ['Time step shrinks too much, it is not feasible to complete it'];
    o.writeMessage(message,'%s\n')
    message = ['I AM STOPPING'];
    o.writeMessage(message,'%s\n')
    terminate = true;
  end
  
  if iter == 0 
    message = 'LAST GMRES REQUIRED 0 ITERATIONS!';    
    o.writeMessage(message,'%s\n');
    message = 'I AM STOPPING!';    
    o.writeMessage(message,'%s\n');
    terminate = true;
  end
else
  oc = curve;
  [~,area,length] = oc.geomProp(X);  
  ea = max(abs(area./o.areaAll-1)); 
  el = max(abs(length./o.lengthAll-1));  
  message = 'All vesicles have been streamed and completed simulation';
  o.writeMessage(message,'%s\n');  
  message = 'SIMULATION SUCCESSFULLY COMPLETED';
  o.writeMessage(message,'%s\n');  
  terminate = true; 
end

% Begin plotting
if o.usePlot
  o.plotData(X,u,Xwalls,XwallsInt,XwallsExt,Xtra,time,rem_ids,allVesViscConts);
  pause(0.01)    
end
% End plotting

% Begin saving data
if o.saveData
  % don't want to save initial small time steps, but still
  % want to check the error in area and length so that the
  % simulation is killed early on if need be
  o.writeData(X,sigma,time,res);
end
% End saving data


% Begin sending messages to log file and console
message = ['GMRES required ' num2str(iter) ...
    ' iterations to couple vesicles and solid walls'];
o.writeMessage(message,'%s\n')

if iflag == 1
  message = 'GMRES DID NOT CONVERGE: didn''t achieve tolerance';
elseif iflag == 2
  message = 'GMRES DID NOT CONVERGE: preconditioner ill-conditioned';
elseif iflag == 3
  message = 'GMRES DID NOT CONVERGE: successive iterates were the same';
end
if iflag ~=0
  o.writeStars
  o.writeMessage(message,'%s\n')
  o.writeStars
end

message = ['t = ' num2str(time,'%4.2e') ...
      ' of T = ' num2str(o.T,'%4.2e')]; 
o.writeMessage(message,'%s\n')
% write the new time step size if doing time adaptivity

message1 = ['Max error in area is   ' num2str(ea,'%4.2e')];
message2 = ['Max error in length is ' num2str(el,'%4.2e')];
o.writeMessage(message1,'%s\n')
o.writeMessage(message2,'%s\n')
%if o.timeAdap
%  message3 = ['Residual is            ' num2str(res,'%4.2e')];
%  o.writeMessage(message3,'%s\n')
%end
o.writeMessage(' ','%s\n')
% End sending data to files and console


if (o.saveData && o.tracers)
  fileName = [o.logFile(1:end-4) 'Tracers.bin'];
  fid = fopen(fileName,'a');
  fwrite(fid,Xtra,'double');
  fclose(fid);
  % save the tracer locations
end

end % outputInfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writePressure(o,press)
% writePressure(press) writes the pressure to a .bin file

output = press;

fileName = [o.logFile(1:end-4) 'Pressure.bin'];
fid = fopen(fileName,'a');
fwrite(fid,output,'double');
fclose(fid);
% save the pressure values

end % writePressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeStress(o,output,flag)
% writeStress(X,output,flag) writes the stress to a .bin file.  flag is
% used to decide if it needs to write the locations or one of the four
% stress components

if (strcmp(flag,'pts') || strcmp(flag,'11'))
  fileName = [o.logFile(1:end-4) 'Stress11.bin'];
  fid = fopen(fileName,'a');
  fwrite(fid,output,'double');
  fclose(fid);
  % save the (1,1) component of the stress
end

if (strcmp(flag,'pts') || strcmp(flag,'12'))
  fileName = [o.logFile(1:end-4) 'Stress12.bin'];
  fid = fopen(fileName,'a');
  fwrite(fid,output,'double');
  fclose(fid);
  % save the (1,2) component of the stress
end

if (strcmp(flag,'pts') || strcmp(flag,'21'))
  fileName = [o.logFile(1:end-4) 'Stress21.bin'];
  fid = fopen(fileName,'a');
  fwrite(fid,output,'double');
  fclose(fid);
  % save the (2,1) component of the stress
end

if (strcmp(flag,'pts') || strcmp(flag,'22'))
  fileName = [o.logFile(1:end-4) 'Stress22.bin'];
  fid = fopen(fileName,'a');
  fwrite(fid,output,'double');
  fclose(fid);
  % save the (2,2) component of the stress
end

end % writeStress

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ea,el] = errors(o,X)
% function [ea,el] = errors(X) computes the errors in area and length
% of the new vesicle position

oc = curve;
[~,a,l] = oc.geomProp(X);
% compute new areas and length

ea = max(abs(a./o.area - 1));
el = max(abs(l./o.length - 1));
% compute error in area and length

end % errors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeStars(o)
% writeStars writes a message of stars to the console and the log file
% depending on verbose and saveData

messageStars = '*********************************************';
o.writeMessage(messageStars,'%s\n')

end % writeStars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeMessage(o,message,format)
% function writeMessage(message,format) appends message to o.fileName
% with format


fid = fopen(o.logFile,'a');
fprintf(fid,format,message);
fclose(fid);

% save to log file
if o.verbose
  disp(message)
end
% write to console


end % writeMessage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotData(o,X,u,Xwalls,XwallsInt,XwallsExt,Xtra,time,rem_ids,allVesViscConts)
    
% plotData(X,u,Xwalls,XwallsInt,XwallsExt,Xtra,time,ea,el,rem_ids)
% plots the current configuration
% with title X is the vesicle position, u is the vesicle velocity,
% Xwalls is the solid wall positions, Xtra is the tracer locations,
% time is the current time, ea and el are the errors in area and length


oc = curve;
% seperate x and y coordinates
[x,y] = oc.getXY(X);

% Plot all vesicles and walls (if there is any)
if ~o.istreaming
figure(1); clf; hold on
plot([x;x(1,:)],[y;y(1,:)],'Color',[1 0.6 0.78],'linewidth',2)
plot([x(:,rem_ids);x(1,rem_ids)],[y(:,rem_ids);y(1,rem_ids)],'r','linewidth',2)
if ~isempty(Xwalls)
  [xwall,ywall] = oc.getXY(Xwalls);
  plot([xwall;xwall(1,:)],[ywall;ywall(1,:)],'k','linewidth',2)
elseif ~isempty(XwallsInt)
  [xwallInt,ywallInt] = oc.getXY(XwallsInt);
  [xwallExt,ywallExt] = oc.getXY(XwallsExt);
  plot([xwallExt;xwallExt(1,:)],[ywallExt;ywallExt(1,:)],'k','linewidth',2)
  plot([xwallInt;xwallInt(1,:)],[ywallInt;ywallInt(1,:)],'k','linewidth',2)
end
else
figure(1); clf; hold on;    
% if streaming occurs, we do not want to plot the frozen ones, plot only 
% active ones, but show the ones with different VCs with different colors
xrem = x(:,rem_ids); yrem = y(:,rem_ids);
VCsOfActVes = allVesViscConts(rem_ids); % VCs of active ones
VCs = unique(allVesViscConts); % get the viscosity contrasts
idVC1s = find(VCsOfActVes==VCs(1)); % the smaller viscosity contrast (supposed to displace)
idVC2s = find(VCsOfActVes==VCs(2));
% first plot the one with smaller VC, which is supposed to displace
plot([xrem(:,idVC1s);xrem(1,idVC1s)],[yrem(:,idVC1s);yrem(1,idVC1s)],'b','linewidth',2)
% then plot the one with higher VC, which is supposed to zig-zag
plot([xrem(:,idVC2s);xrem(1,idVC2s)],[yrem(:,idVC2s);yrem(1,idVC2s)],'r','linewidth',2)
% plot the walls
[xwallInt,ywallInt] = oc.getXY(XwallsInt);
[xwallExt,ywallExt] = oc.getXY(XwallsExt);
plot([xwallExt;xwallExt(1,:)],[ywallExt;ywallExt(1,:)],'k','linewidth',2)
plot([xwallInt;xwallInt(1,:)],[ywallInt;ywallInt(1,:)],'k','linewidth',2)
end
% Plot a tracker point at the first and middle point
% so that we can observe features like trank-treading
if o.track
  plot(x,y,'bo','markersize',20)
%  plot(x(1,:),y(1,:),'b.','markersize',20)
%  plot(x(N/2+1,:),y(N/2+1,:),'b.','markersize',20)
end

% plot velocity field on surface of vesicles
if o.quiver
  quiver(X(1:N,:),X(N+1:2*N,:),u(1:N,:,end),u(N+1:2*N,:,end))
end

% plot the tracers
if o.tracers
  [xtra,ytra] = oc.getXY(Xtra);
  plot(xtra,ytra,'bo')
  
end

titleStr = ['t = ' num2str(time,'%4.2e')];
title(titleStr)
% Title
axis equal
if ~isempty(o.axis)
  axis(o.axis)
end
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ycolor','w')
set(gca,'ycolor','w')

end % plotData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotnoVesTracersData(o,XwallsInt,XwallsExt,Xtra,epsilon)
% plotData(X,u,Xwalls,Xtra,time,ea,el) plots the current configuration
% with title X is the vesicle position, u is the vesicle velocity,
% Xwalls is the solid wall positions, Xtra is the tracer locations,
% time is the current time, ea and el are the errors in area and length
oc = curve;

[extWallx,extWally] = oc.getXY(XwallsExt);
[intWallx,intWally] = oc.getXY(XwallsInt);

figure(1);clf;
axes('position',[0 0.4 1 .5])
box on
plot([extWallx;extWallx(1)],[extWally;extWally(1)],'k','linewidth',1)
hold on
xwall = interpft(intWallx,96);
ywall = interpft(intWally,96);
vec1 = [xwall;xwall(1,:)];
vec2 = [ywall;ywall(1,:)];
plot(vec1,vec2,'k','linewidth',1)
h1 = fill(vec1,vec2,'k');
set(h1,'edgecolor','k')

nTra = size(Xtra,2)/2;
for traIdx = 1:nTra
  plot(Xtra(1:20:end,traIdx),Xtra(1:20:end,traIdx+nTra),'b','linewidth',1.5);
end

titleStr = ['\epsilon = ' num2str(epsilon)];
title(titleStr,'FontSize',28,'FontName','Palatino')
% Title

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
axis equal
if ~isempty(o.axis)
  axis(o.axis)
end

end % plotData


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(o,X,sigma,time,res)
% writeData(X,sigma,ea,el,time,res) writes the position, tension,
% errors, and time to a binary file.  Matlab can later read this file
% to postprocess the data
 
oc = curve;
[x,y] = oc.getXY(X);
output = [x(:);y(:);sigma(:);time];
% format that postProcess/loadfile.m reads the output
fid = fopen(o.dataFile,'a');
fwrite(fid,output,'double');
fclose(fid);

if o.timeAdap
  fileName = [o.logFile(1:end-4) 'Res.dat'];
  fid = fopen(fileName,'a');
  fprintf(fid,'%10.5e\n',res);
  fclose(fid);
end
%% write the reisdual to a seperate dat file.

end % writeData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeDataWithEta(o,X,sigma,eta,RS,ea,el,time,res)
% writeDatawithEta(X,sigma,eta,ea,el,time,res) writes the position,
% tension, density function on solid walls, errors, and time to a
% binary file.  Matlab can later read this file to postprocess the
% data.  Note that sigma, eta, and RS can be post processed given X
% using computeSigAndEta.  This gives a more accurate velocity field
% for studing problems like an advection-diffusion solver in a couette
% apparatus (see Gokberk's work)
 
oc = curve;
[x,y] = oc.getXY(X);
[etax,etay] = oc.getXY(eta);
%output = [x(:);y(:);sigma(:);etax(:);etay(:);ea;el;time];
% format that postProcess/loadfile2.m reads the output

file = [o.dataFile(1:end-4) '2.bin'];
output = [x(:);y(:);sigma(:);etax(:);etay(:);RS(:);ea;el;time];
% format that postProcess/loadfile2.m reads the output
fid = fopen(file,'a');
fwrite(fid,output,'double');
fclose(fid);

end % writeData


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function summary(o,matvecs,derivs,fmms,repuls,ncollVes,...
        ncollWal,naccept,nreject,countGMRES,totTime)
% summary(o,matvecs,derivs,fmms,repuls,nitrepar,ncollVes,...
%         ncollWal,naccept,nreject,countGMRES,totTime)

% print the total number of matvecs
message = ['Number of times in TimeMatVec is ' ...
    num2str(matvecs)];
o.writeStars
o.writeMessage(message,'%s\n')
o.writeStars

% print the total number of times differential operators are
% formed in matrix form
message = ['Number of times in computeDeriv is ' ...
    num2str(derivs)];
o.writeStars
o.writeMessage(message,'%s\n')
o.writeStars

% print the total number of fmm calls
message = ['Number of times calling Stokes FMM is ' ...
    num2str(fmms)];
o.writeStars
o.writeMessage(message,'%s\n')
o.writeStars

% print the number of times repulsion is actively used (includes SDC)
if o.repulsion
  message = ['Number of times having non-zero repulsion is ' ...
  num2str(repuls)];
  o.writeStars
  o.writeMessage(message,'%s\n')
  o.writeStars   
end


% print the number of vesicle-vesicle and vesicle-wall collisions
if o.timeAdap
  message = ['Number of rejected time steps when vesicle-vesicle collision occurs is ' ...
      num2str(ncollVes)];
  o.writeStars
  o.writeMessage(message,'%s\n')
  o.writeStars
  message = ['Number of rejected time steps when vesicle-wall collision occurs is ' ...
      num2str(ncollWal)];
  o.writeStars
  o.writeMessage(message,'%s\n')
  o.writeStars
end

%print the total number of GMRES iterations
message = ['Number of total GMRES iterations is ' ...
    num2str(countGMRES)];
o.writeStars
o.writeMessage(message,'%s\n')
o.writeStars

% print the total number of accepted and rejected time steps 
if o.timeAdap  
  message = ['Number of accepted time steps is ' num2str(naccept)];
  o.writeStars
  o.writeMessage(message,'%s\n')
  message = ['Number of rejected time steps is ' ...
      num2str(nreject)];
  o.writeStars
  o.writeMessage(message,'%s\n')
  o.writeStars
end

% Save total time spent
message = ['Total time is ' num2str(totTime,'%5.2e') ];
o.writeMessage(message,'%s\n')
o.writeStars


end % summary


end % methods


end % classdef

