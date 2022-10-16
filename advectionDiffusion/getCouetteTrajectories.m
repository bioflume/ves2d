function getCouetteTrajectories

V1 = @(x,y)  y/3 .* (1-400./(x.^2+y.^2)); % x-velocity
V2 = @(x,y) -x/3 .* (1-400./(x.^2+y.^2)); % y-velocity

%% ------------------------------------------------------------------------
% PARAMETERS and OPTIONS
% -------------------------------------------------------------------------
% Spatial Discretization
params.Ntheta = 512;
params.Nradii = 127;
params.radii = [10; 20];

% Temporal Discretization
params.Th = 100;
params.deltaT_diff = 0.04;
params.deltaT_adv = params.deltaT_diff/2;

% Trajectories
% xcoords = [10:1.25:20]';
% params.nTrajPoints = size(xcoords(2:end-1));

params.nTrajPoints = 100;

options.velInterpScheme = 'spline';
options.saveData = 1;

files.trajectories = '/scratch/gokberk/Case111/Couette_Trajectories.bin';
files.LogFile      = '/scratch/gokberk/Case111/TrajectoriesLog.log';

%% ------------------------------------------------------------------------
% SPACE DISCRETIZATION
% -------------------------------------------------------------------------
[x,y,r1D,theta,s,R1D,l] = generateGrid(params.Ntheta,params.Nradii,params.radii);
[THETA,R] = meshgrid([theta;2*pi],r1D);

thetaEx = [theta(end-9:end)-2*pi ;theta; theta(1:10)+2*pi]; % extend theta due to periodicity
[thetg,sg] = meshgrid(thetaEx,s); % for new Cheb grid (s does not contain 0 and pi)

trajPoints = floor((params.Nradii-2).*rand(params.nTrajPoints,1) + 2);
xTrajPoints = x(trajPoints,1);
yTrajPoints = y(trajPoints,1);

%% ------------------------------------------------------------------------
% TIME DISCRETIZATION
% -------------------------------------------------------------------------
ntime = params.Th/params.deltaT_adv; % Number of time intervals for diffusion
time  = linspace(0,params.Th,ntime+1);

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
Th = params.Th;


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

 Vx_tPoints = V1(xTrajPoints,yTrajPoints);
 Vy_tPoints = V2(xTrajPoints,yTrajPoints);
  
 % Find midpoints
 x_mid = xTrajPoints + params.deltaT_adv * Vx_tPoints;
 y_mid = yTrajPoints + params.deltaT_adv * Vy_tPoints;
 
 Vx_MidtPoints = V1(x_mid,y_mid);
 Vy_MidtPoints = V2(x_mid,y_mid);
 
 xTrajPointsNew = xTrajPoints + 0.5 * params.deltaT_adv * (Vx_tPoints + Vx_MidtPoints);
 yTrajPointsNew = yTrajPoints + 0.5 * params.deltaT_adv * (Vy_tPoints + Vy_MidtPoints);

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
