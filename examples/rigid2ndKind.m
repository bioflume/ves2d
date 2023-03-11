clear; clc;

addpath ../src/

%% Problem setup
% Fluid properties
viscosity = 1; % Pa.s -- 1000x water's viscosity

% Body geometry
radius = 5; % micro meters
Npoints = 64;

% Body shape discretization assuming a sphere
theta = [0:Npoints-1]'/Npoints * 2 * pi;
body_x = radius*cos(theta);
body_y = radius*sin(theta);

% Force and torque under body moves
ext_force = [1; 0]; % pico Newton
ext_torque = [0];

% Time scale -- Stokes' law
% Fd = 6 * pi * mu * R* v
exp_vel = ext_force ./ (6 * pi * viscosity * radius); % in um/sec
% time horizon is reached when 10x of body's diameter is traveled
time_horizon = (10*2*radius)/norm(exp_vel); 
time_scale = 2*radius/norm(exp_vel); 
dt = 1E-4 * time_scale; 

% Call kernels
ker = kernels(Npoints);

%% Take time steps
time = 0;
nsteps = 0;
X = [body_x; body_y];
center = [0; 0];
orientation = 0; % angle
traction = zeros(2*Npoints,1);
U = zeros(3,1);
while time < time_horizon
  % update time
  time = time + dt; 
  nsteps = nsteps + 1;
  
  % build body object
  bb = body(X, radius, center, orientation,U);
  bb.traction = traction; 

  % calculate the K matrix
  K = bb.calc_K_matrix();
  KT = bb.calc_KT_matrix(bb.sa);

  % calculate the single and double layer matrices
  SLP = ker.stokesSLmatrixAlpert(bb,viscosity);  
  DLP = ker.stokesDLmatrix(bb);
  DLPT = ker.stokesDLTmatrix(bb);

  % form RHS
  RHS = zeros(2*Npoints+3,1);
  RHS(2*Npoints+1:end) = -[ext_force;ext_torque];
    
  % form the LHS matrix
  MAT = [SLP              -0.5*K-DLP*K;...
        -0.5*KT-KT*DLPT   zeros(3)];

  % solve the system
  sol = gmres(MAT,RHS);

  % If we want to check the N matrix:
  Ktil = 0.5*K+DLP*K;
  KtilT = 0.5*KT+KT*DLPT;
  Nmat = (KtilT*(SLP\eye(size(SLP)))*Ktil)\eye(3);

  % dissect the solution
  traction = sol(1:2*Npoints);
  ui = sol(2*Npoints+1:2*Npoints+2);
  wi = sol(2*Npoints+3);
  U = [ui;wi];
  

  figure(2);clf;
  plot([X(1:end/2);X(1)],[X(end/2+1:end);X(end/2+1)],'linewidth',2)
  hold on
  quiver(X(1:end/2),X(end/2+1:end),traction(1:end/2),traction(end/2+1:end))
  axis equal

  % update the position and orientation
  center = bb.center;
    
  X0 = [X(1:Npoints)-center(1); X(Npoints+1:end)-center(2)];
  center = center + ui*dt;

  xnew = center(1)+cos(wi*dt)*X0(1:Npoints)+sin(wi*dt)*X0(Npoints+1:end);
  ynew = center(2)-sin(wi*dt)*X0(1:Npoints)+cos(wi*dt)*X0(Npoints+1:end);
  X = [xnew; ynew];
  
  orientation = orientation + wi*dt;
  
  % Display the results  
  display_results(ui,wi,exp_vel,time,X,radius)

  pause

end % end while time < time_horizon


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_results(ui, wi, exp_vel, time, X, radius)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Time = ' num2str(time) 's'])
disp(['Expected velocity is ' num2str(reshape(exp_vel,1,[])) 'um/s'])
disp(['Velocity is found: ' num2str(reshape(ui,1,[])) 'um/s'])
disp(['Angular velocity is found: ' num2str(wi)])

figure(1); clf;
vecx = [interpft(X(1:end/2),128); X(1)];
vecy = [interpft(X(end/2+1:end),128); X(end/2+1)];

plot(vecx, vecy, 'k');
fill(vecx, vecy,'k');
axis equal
xlim([-2*radius 20*radius])
ylim([-2*radius 2*radius])

grid off
box on
title(['Time = ' num2str(time) 's'])
pause(0.1)

end



