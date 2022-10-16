function getTrajectories
%% ------------------------------------------------------------------------
% PARAMETERS and OPTIONS
% -------------------------------------------------------------------------
% Spatial Discretization
params.Ntheta = 512;
params.Nradii = 127;
params.radii = [10; 20];

% Temporal Discretization
params.Th = 90;
params.deltaT_diff = 0.04;
params.deltaT_adv = params.deltaT_diff/2;

% Trajectories
% xcoords = [10:1.25:20]';
% params.nTrajPoints = size(xcoords(2:end-1));

params.nTrajPoints = 100;

options.velInterpScheme = 'spline';
options.saveData = 1;

% files.velocity     = '/scratch/gokberk/Case211/Case211TA_Vel128x512_p01.bin';
% files.trajectories = '/scratch/gokberk/Case211/Case211_UniformTrajectories.bin';
% files.LogFile      = '/scratch/gokberk/Case211/TrajectoriesLog.log';

% files.velocity     = '/scratch/gokberk/LargestSim/couette150_Vel256x1024_intp.bin';
% files.trajectories = '/scratch/gokberk/LargestSim/LargSim_UniformTrajectories.bin';
% files.LogFile      = '/scratch/gokberk/LargestSim/TrajectoriesLog.log';

% files.velocity     = '/scratch/gokberk/Case411/Case411TA_Vel256x1024_Quick_p01.bin';
% files.trajectories = '/scratch/gokberk/Case411/Case411_UniformTrajectories.bin';
% files.LogFile      = '/scratch/gokberk/Case411/TrajectoriesLog.log';

% files.velocity     = '/scratch/gokberk/Case111/Case111TA_48N_Vel128x512_intp.bin';
% files.trajectories = '/scratch/gokberk/Case111/Case111_UniformTrajectories.bin';
% files.LogFile      = '/scratch/gokberk/Case111/TrajectoriesLog.log';

% files.velocity     = '/scratch/gokberk/Case112/Case112TA_Vel128x512NewRec_intp.bin';
% files.trajectories = '/scratch/gokberk/Case112/Case112_UniformTrajectories.bin';
% files.LogFile      = '/scratch/gokberk/Case112/TrajectoriesLog.log';
% params.readRemove  = 0;

files.velocity     = '/scratch/gokberk/Case115/Case115TA_Vel128x512_intp.bin';
files.trajectories = '/scratch/gokberk/Case115/Case115_UniformTrajectories.bin';
files.LogFile      = '/scratch/gokberk/Case115/TrajectoriesLog.log';
params.readRemove  = 1000;

%params.readRemove = 10000;

%% ------------------------------------------------------------------------
% SPACE DISCRETIZATION
% -------------------------------------------------------------------------
[x,y,r1D,theta,s,R1D,l] = generateGrid(params.Ntheta,params.Nradii,params.radii);
[THETA,R] = meshgrid([theta;2*pi],r1D);

thetaEx = [theta(end-9:end)-2*pi ;theta; theta(1:10)+2*pi]; % extend theta due to periodicity
[thetg,sg] = meshgrid(thetaEx,s); % for new Cheb grid (s does not contain 0 and pi)

% trajPoints = floor((params.Nradii-2).*rand(params.nTrajPoints,1) + 2);

xUniform = linspace(params.radii(1),params.radii(2),params.nTrajPoints+2);
xTrajPoints = xUniform(2:end-1)';
yTrajPoints = zeros(params.nTrajPoints,1);

% xTrajPoints = xcoords(2:end-1);
% yTrajPoints = zeros(7,1);

%% ------------------------------------------------------------------------
% TIME DISCRETIZATION
% -------------------------------------------------------------------------
ntime = params.Th/params.deltaT_adv; % Number of time intervals for diffusion
time  = linspace(0,params.Th,ntime+1);

%% ------------------------------------------------------------------------
% LOAD VELOCITY FILE
% -------------------------------------------------------------------------
fidVel = fopen(files.velocity,'r');
valVel = fread(fidVel,2,'double');
% Read and remove velocity 
 
for i = 1 : params.readRemove
    vel = fread(fidVel,valVel(1),'double');
end
clear vel

% ---------------------------------------------------------------------
%% Initialize Data Files
% ---------------------------------------------------------------------
if options.saveData
    fid = fopen(files.trajectories,'w');
    fwrite(fid,[params.Nradii+1;params.Ntheta+1;ntime+1],'double');
    fwrite(fid,params.Th,'double');
    fwrite(fid,params.radii,'double');
    fclose(fid);

    fid = fopen(files.LogFile,'w');
    dateTime = datestr(clock);
    fprintf(fid,'%s\n',['Finding trajectories is started on: ' dateTime]);
    fprintf(fid,'%s\n',['Velocity file: ' files.velocity]);
    fprintf(fid,'%s\n',['Nradial x Ntheta = ' num2str(params.Nradii+1) 'x' num2str(params.Ntheta) ]);
    fprintf(fid,'%s\n',['Time horizon = ' num2str(params.Th) '. Time step for advection = ' num2str(params.deltaT_adv)]);
    fclose(fid);


    fid = fopen(files.trajectories,'w');
    fwrite(fid,params.nTrajPoints,'double');
    fwrite(fid,ntime+1,'double');
    fclose(fid);
    writeData(files.trajectories,[xTrajPoints;yTrajPoints]);
end
% ---------------------------------------------------------------------

%% Trace the particles
N_radii = params.Nradii;
N_theta = params.Ntheta;
radius = params.radii;
Th = params.Th;

% Read velocity at every time step (deltaT_diff/deltaT_adv/2 steps)
vel0 = fread(fidVel,valVel(1),'double');

% Velocity components at t_n and on the regular grid
Vx0 = reshape(vel0(1:end/2),(N_radii+1),(N_theta+1));
Vy0 = reshape(vel0(end/2+1:end),(N_radii+1),(N_theta+1));
Vx0 = Vx0(:,1:end-1);
Vy0 = Vy0(:,1:end-1);
% Enforce periodicity in theta
Vx0Ex  = [Vx0(:,end-9:end) Vx0 Vx0(:,1:10)];
Vy0Ex  = [Vy0(:,end-9:end) Vy0 Vy0(:,1:10)];

remove = 1;

% PLOT THE APPARATUS
%--------------------------------------
xwalls1 = 20*cos((0:360)*2*pi/360);
ywalls1 = 20*sin((0:360)*2*pi/360);
xwalls2 = 10*cos((0:360)*2*pi/360);
ywalls2 = 10*sin((0:360)*2*pi/360);
plot(xwalls1,ywalls1,'k','linewidth',2)
hold on
plot(xwalls2,ywalls2,'k','linewidth',2)
axis square
%--------------------------------------

plot(xTrajPoints,yTrajPoints,'.k','markersize',6)

for tt = 1 : ntime
    
 disp([num2str(time(tt+1)) ' unit of ' num2str(Th) ' :'])
 disp('--------------------------------------------------------------')

 tic
 if remove
    % Remove deltaT = 0.01;
    velRemove = fread(fidVel,valVel(1),'double');
    clear velRemove
 end

 % Map to the uniform theta-s grid to interpolate velocity at initial points and t_n 
 [s0,thet0] = arbiTOreg(xTrajPoints,yTrajPoints,radius,R1D);

 Vx_tPoints = interp2(thetg,sg,Vx0Ex,thet0,s0,options.velInterpScheme);
 Vy_tPoints = interp2(thetg,sg,Vy0Ex,thet0,s0,options.velInterpScheme);
 
 % Find midpoints
 x_mid = xTrajPoints + params.deltaT_adv * Vx_tPoints;
 y_mid = yTrajPoints + params.deltaT_adv * Vy_tPoints;
 
 % Map to the uniform theta-s grid to interpolate velocity at midPoint and t_n+1 
 [s_mid,thet_mid] = arbiTOreg(x_mid,y_mid,radius,R1D);
 
 % Read velocity at every time step (deltaT_diff/deltaT_adv/2 steps)
 vel1 = fread(fidVel,valVel(1),'double');
 
 % Velocity components at t_n+1 and on the regular grid
 Vx1 = reshape(vel1(1:end/2),(N_radii+1),(N_theta+1));
 Vy1 = reshape(vel1(end/2+1:end),(N_radii+1),(N_theta+1));
 Vx1 = Vx1(:,1:end-1);
 Vy1 = Vy1(:,1:end-1);
 
 % Enforce periodicity in theta
 Vx1Ex  = [Vx1(:,end-9:end) Vx1 Vx1(:,1:10)];
 Vy1Ex  = [Vy1(:,end-9:end) Vy1 Vy1(:,1:10)];
 
 Vx_MidtPoints = interp2(thetg,sg,Vx1Ex,thet_mid,s_mid,options.velInterpScheme);
 Vy_MidtPoints = interp2(thetg,sg,Vy1Ex,thet_mid,s_mid,options.velInterpScheme);
 
 xTrajPointsNew = xTrajPoints + 0.5 * params.deltaT_adv * (Vx_tPoints + Vx_MidtPoints);
 yTrajPointsNew = yTrajPoints + 0.5 * params.deltaT_adv * (Vy_tPoints + Vy_MidtPoints);

 Vx0Ex  = Vx1Ex;
 Vy0Ex  = Vy1Ex;
 xTrajPoints = xTrajPointsNew;
 yTrajPoints = yTrajPointsNew;
 
 plot(xTrajPoints,yTrajPoints,'.k','markersize',6)
 
 tim = toc;
 disp(['Elapsed time is ' num2str(tim) 'sec'])
  
 if options.saveData
     keepDiary(files.LogFile,ntime,tt,tim);
     writeData(files.trajectories,[xTrajPoints;yTrajPoints]);
 end
   
end
 
end
% #########################################################################
%% MAPPING FROM ARBITRARY DOMAIN (r,theta) TO UNIFORM DOMAIN (s,theta)
% #########################################################################
function [sp,thetp] = arbiTOreg(xp,yp,radius,R1D)
rinr = radius(1);
rout = radius(2);

% Find theta at the arbitrary point
thetp = atan2(yp,xp);
thetp(:,ceil(end/2)+1:end) = thetp(:,ceil(end/2)+1:end) + 2*pi;
thetp(thetp>=2*pi) = mod(thetp(thetp>=2*pi),2*pi);
thetp(thetp<0)    = mod(thetp(thetp<0),2*pi);

% Find r at the arbitrary point [r0, r1]
rp = (xp.^2+yp.^2).^(0.5);

% Find R at the arbitrary point [-1, 1]
% Rp = 2*(rp-rinr)/(rout-rinr) - 1;
Rp = (R1D(end)-R1D(1))*(rp - rinr)/(rout - rinr) + R1D(1);

Rp(Rp>1) = 1;
Rp(Rp<-1) = -1;

% Find s at the arbitrary point
sp = acos(-Rp);
end

% #########################################################################
%% WRITE DATA to BIN FILE
% #########################################################################
function writeData(fileName,dataArr)
fid = fopen(fileName,'a');
fwrite(fid,dataArr,'double');
fclose(fid);
end

% #########################################################################
%% KEEP DIARY
% #########################################################################
function keepDiary(fileName,ntime,tst,tim)
fid = fopen(fileName,'a');
fprintf(fid,'%s\n',['Elapsed Time at step = ' num2str(tst) ' of ' num2str(ntime) ' = ' num2str(tim) ' sec']);

fclose(fid);
end 
 