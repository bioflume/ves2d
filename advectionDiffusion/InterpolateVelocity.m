load Ves40Vel_254x1024
vel = vel(:,31:151);
pause
N_theta = 1024; 
N_radii = 255;  
% Fourier Grid
theta = [0:N_theta-1]'/N_theta*2*pi;
% Chebyshev Grid
s = [0:N_radii]*pi/N_radii;
R1D = -cos(s);

geometry = [10 20];
% Linear Map
r1D = (R1D+1)*(geometry(2)-geometry(1))/2 + geometry(1);
% MESH 
[theta, R] = meshgrid(theta,R1D); % in matrix form
r = (R+1)*(geometry(2)-geometry(1))/2 + geometry(1); % [r1,r2] in matrix form
x = r.*cos(theta); y = r.*sin(theta); % in cartesian coordinates

% Time Interpolation
Th = 12;
delta_t_old = 0.1;
delta_t_new = 0.01;

nTimeOld = Th/delta_t_old;
nTimeNew = Th/delta_t_new;

NoP = N_theta * (N_radii-1);

Told = linspace(0,Th,nTimeOld+1);
Tnew = linspace(0,Th,nTimeNew+1);

VELOCITY = zeros(2*NoP,nTimeNew+1);

% log file
logFile = 'Ves40InterpLogFile.log';

if 1 
    fid = fopen(logFile,'w');
    fprintf(fid,'%s\n',['# of vesicles = ' num2str(40) ', Grid Size = ', num2str(N_radii+1) 'x' num2str(N_theta)]);
    fprintf(fid,'%s\n',['DeltaT_old = ' num2str(delta_t_old) ', DeltaT_new = ' num2str(delta_t_new)]);
    fclose(fid);
end


for pp = 1 : 2*NoP
    tic
    VELOCITY(pp,:) = interp1(Told,vel(pp,:),Tnew,'spline');
    tim = toc;

    fid = fopen(logFile,'a');
    fprintf(fid,'%s\n',['Elapsed Time at point = ' num2str(pp) ' is ' num2str(tim) ' sec']);
    fclose(fid);
end

vel_new = zeros(2*(N_radii+1)*(N_theta+1),nTimeNew+1);

for tt = 1 : nTimeNew+1
    
    Vt = VELOCITY(:,tt);
    Vx = Vt(1:(N_radii-1)*(N_theta));
    Vy = Vt((N_radii-1)*(N_theta)+1:end);
    Vx = reshape(Vx,N_radii-1,N_theta);
    Vy = reshape(Vy,N_radii-1,N_theta);
    
    Vxt = zeros(N_radii+1,N_theta+1); Vyt = Vxt;
    
    Vxt(2:end-1,1:N_theta) = Vx;
    Vxt(1,1:N_theta) = -10*sin(theta(1,:));
    Vxt(end,1:N_theta) = zeros(1,N_theta);
    Vxt(:,end) = Vxt(:,1);
    
    Vyt(2:end-1,1:N_theta) = Vy;
    Vyt(1,1:N_theta) = 10*cos(theta(1,:));
    Vyt(end,1:N_theta) = zeros(1,N_theta);
    Vyt(:,end) = Vyt(:,1);
    
    vel_new(:,tt) = [Vxt(:);Vyt(:)];
    
end


% save('Vel2Vesic_256x1025_P01.mat','vel_new','-v7.3')
