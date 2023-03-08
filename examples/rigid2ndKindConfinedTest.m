clear; clc;

addpath ../src/
iflag = 1;

%% Problem setup
% Fluid properties
viscosity = 1; % Pa.s -- 1000x water's viscosity

% Body geometry
radius = 1; % micro meters
Npoints = 64;

% Wall geometry
NpointsWall = 512;
if 1
Xwalls = generateWall(NpointsWall, 20*radius, 10*radius);
else
theta = [0:NpointsWall-1]'/NpointsWall * 2 * pi;
Xwalls = [10*radius*cos(theta); 10*radius*sin(theta)];
end

kerWall = kernels(NpointsWall);
wall = walls(Xwalls, [0 0], zeros(size(Xwalls)));
wallDLP = kerWall.stokesDLmatrixWall(wall);
wallN0 = kerWall.stokesN0matrix(wall);
wallPrecond = wallsPrecond(wall, wallDLP, wallN0); % inverse of wall matrices

% Body shape discretization assuming a sphere
theta = [0:Npoints-1]'/Npoints * 2 * pi;
body_x = radius*cos(theta);
body_y = radius*sin(theta);

% Force and torque under body moves
ext_force = [5; 0]; % pico Newton
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
surf_vel = zeros(2*Npoints,1);

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

    
  % FIRST CALCULATE TRACTION DUE TO FORCE
  % traction = solve K adjoint * traction = Force using pseudoinverse
  if nsteps == 1
%   traction = KT\([ext_force; ext_torque]);
  traction = (pinv(K*KT))*(K*[ext_force; ext_torque]);
  end
  
  if iflag
  figure(1); clf;
  plot(bb.X(1:end/2),bb.X(end/2+1:end))
  hold on
  quiver(bb.X(1:end/2),bb.X(end/2+1:end),traction(1:end/2),traction(end/2+1:end))
  axis equal
  title('Traction on body--1')
  pause
  end
  % CALCULATE VELOCITY DUE TO TRACTION ON THE WALL
  % uwall = -exactStokesDL*surf_vel + exactStokesSLP*traction  
  [~,uwall] = ker.exactStokesSL(bb, viscosity, traction, [], Xwalls, 1);
  if iflag
  uMag = mean(sqrt(uwall(1:end/2).^2 + uwall(end/2+1:end).^2));    
  disp(['Velocity due to traction: ' num2str(uMag)]);
  end

  if nsteps > 1
  [~,uwall2] = ker.exactStokesDL(bb,K*U,[],Xwalls,1);
  uwall = uwall + uwall2;
  if iflag
  uMag = mean(sqrt(uwall2(1:end/2).^2 + uwall2(end/2+1:end).^2));
  disp(['Velocity due to u: ' num2str(uMag)]);
  end
  end
  if iflag
  figure(2); clf;
  plot(Xwalls(1:end/2),Xwalls(end/2+1:end))
  hold on
  quiver(Xwalls(1:end/2),Xwalls(end/2+1:end),uwall(1:end/2),uwall(end/2+1:end))
  axis equal
  title('Velocity on the wall')
  pause
  end

  % SOLVE FOR WALL DENSITY
  wallDensity = wallPrecond*(-uwall);
  
  % CALCULATE VELOCITY ON THE BODY DUE TO WALL DENSITY
  [~,ubody] = ker.exactStokesDL(wall,wallDensity,[],bb.X,1);
  
  if iflag
  figure(1);  hold on;
  quiver(bb.X(1:end/2),bb.X(end/2+1:end),ubody(1:end/2),ubody(end/2+1:end))
  title('Velocity on the body')
  uMag = mean(sqrt(ubody(1:end/2).^2 + ubody(end/2+1:end).^2));
  disp(['Velocity mag on the body due to wall: ' num2str(uMag)]);
  pause
  end

  % THEN SOLVE FOR BODY VELOCITY AND TRACTION

  % form RHS
  RHS = zeros(2*Npoints+3,1);
  RHS(1:2*Npoints) = -ubody;
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

  % update the position and orientation
  center = bb.center;
    
  X0 = [X(1:Npoints)-center(1); X(Npoints+1:end)-center(2)];
  center = center + ui*dt;

  xnew = center(1)+cos(wi*dt)*X0(1:Npoints)+sin(wi*dt)*X0(Npoints+1:end);
  ynew = center(2)-sin(wi*dt)*X0(1:Npoints)+cos(wi*dt)*X0(Npoints+1:end);
  X = [xnew; ynew];
  
  orientation = orientation + wi*dt;
  
  % Display the results  
  display_results(ui,wi,exp_vel,time,X,Xwalls,radius)

%   pause

end % end while time < time_horizon


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_results(ui, wi, exp_vel, time, X, Xwalls,radius)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Time = ' num2str(time) 's'])
disp(['Expected velocity is ' num2str(reshape(exp_vel,1,[])) 'um/s'])
disp(['Velocity is found: ' num2str(reshape(ui,1,[])) 'um/s'])
disp(['Angular velocity is found: ' num2str(wi)])

figure(1); clf;
vecx = [interpft(X(1:end/2),128); X(1)];
vecy = [interpft(X(end/2+1:end),128); X(end/2+1)];

plot(vecx, vecy, 'k');
hold on;
fill(vecx, vecy,'k');
axis equal


grid off
box on
title(['Time = ' num2str(time) 's'])

if ~isempty(Xwalls)
vecx = [Xwalls(1:end/2); Xwalls(1)];
vecy = [Xwalls(end/2+1:end); Xwalls(end/2+1)];
plot(vecx, vecy, 'k','linewidth',2)
else
xlim([-2*radius 20*radius])
ylim([-2*radius 2*radius])
end
pause(0.1)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = generateWall(N,a ,b)
order = 10;
  % parameters for the boundary
  
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
%   t = (0:N-1)'*2*pi/N;
  % Parameterize t so that geometry is closer to 
  % equispaced in arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  X = [x;y];
  % rounded off cylinder.  a and b control the length and height 
  % and order controls the regularity

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mat = wallsPrecond(walls, wallDLP, wallN0)
% wallsPrecond(walls) computes the matrix which is the 
% exact inverse of
% the double-layer potential for stokes flow in a bounded domain.  Used
% in the preconditioner for vesicle simulations and capsules.m/computeEta
% which computes eta and RS when there is no vesicle.


Nbd = walls.N;

% Allocate space for blocks of matrix that carries the double- layer
% potential, rotlets, and stokeslets to the velocity and the conditions
% in (A4) and (A5) in Rahimian et al.
M11 = zeros(2*Nbd,2*Nbd);


% Self interaction terms with the jump coming from the double layer
% potential
M11(1:2*Nbd,1:2*Nbd) = M11(1:2*Nbd,1:2*Nbd) + wallN0(:,:,1);
jump = - 1/2; 
istart = 1;
iend = 2*Nbd;
M11(istart:iend,istart:iend) = M11(istart:iend,istart:iend) + ...
    jump*eye(2*Nbd) + wallDLP(:,:,1);


% invert the matrix
Mat = M11 \ eye(2*Nbd);

    
end % wallsPrecond