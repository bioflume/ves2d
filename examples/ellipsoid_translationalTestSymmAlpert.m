function error = ellipsoid_translationalTestSymmAlpert(N,wallRadScale)

addpath ../src/
iflag = 1;

oc = curve;

%% Problem setup
% Fluid properties
viscosity = 1; % Pa.s -- 1000x water's viscosity

% Body geometry
radius_x = 1; % micro meters
radius_y = 2;
Npoints = N;
radius = max([radius_x;radius_y]);

% Wall geometry
NpointsWall = (N/radius)*wallRadScale;
wall_radius = wallRadScale*radius;
theta = [0:NpointsWall-1]'/NpointsWall * 2 * pi;
Xwalls = [wall_radius*cos(theta); wall_radius*sin(theta)];


% Wall kernels
kerWall = kernels(NpointsWall);
wall = walls(Xwalls, [0 0], zeros(size(Xwalls)));
wallSLP = kerWall.stokesSLmatrixAlpertWeightless(wall,viscosity);
% Symmetric Alpert:
wallSLP = 0.5 * (wallSLP + wallSLP');

Hwall = [wall.sa; wall.sa] * 2 * pi / wall.N;

% Body shape discretization assuming a sphere
theta = [0:Npoints-1]'/Npoints * 2 * pi;
body_x = radius_x*cos(theta);
body_y = radius_y*sin(theta);
% X = [body_x; body_y];
% for it = 1 : 100
% [X,~,~] = oc.redistributeArcLength(X);
% end
% body_x = X(1:end/2);
% body_y = X(end/2+1:end);

%% work out the analytical solution
given_force = 1;
R1 = radius;
R2 = wall_radius;

h = (R1^2 + R2^2)/(R1^2*(log(R2)+1) + R2^2*log(R2)-R2^2 - (R1^2+R2^2)*log(R1));
am2 = R1^2 * R2^2 / (2 * (R1^2*(log(R2)+1)+R2^2*log(R2)-R2^2-(R1^2+R2^2)*log(R1)));
a0 = (R1^2*(2*log(R2)+1) + 2*R2^2*log(R2)-R2^2)/(2*R1^2*(log(R2)+1)+R2^2*log(R2)-R2^2-(R1^2+R2^2)*log(R1));
a2 = 1/(2*(R1^2*(log(R2)+1)+R2^2*log(R2)-R2^2-(R1^2+R2^2)*log(R1)));


expected_velocity = given_force/(4*pi*viscosity*h);
% Force and torque under body moves
ext_force = [given_force; 0]; % pico Newton
ext_torque = [0];

% Time scale
time_horizon = (10*2*radius)/norm(expected_velocity); 
time_scale = 2*radius/norm(expected_velocity); 
dt = 1E-3 * time_scale; 
time_horizon = dt;
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
  KT = bb.calc_KT_matrix_weightless(bb.sa);
  H = bb.calc_jac_matrix();

  % calculate the single and double layer matrices
  SLP = ker.stokesSLmatrixAlpertWeightless(bb,viscosity);  
  SLP = 0.5 * (SLP + SLP');
  
  DLP = ker.stokesDLmatrixWeightless(bb);
  DLPT = ker.stokesDLTmatrixWeightless(bb);

  % calculate the interaction matrices
  SLP_bo = ker.stokesSLmatrixInteractionWeightless(wall, bb, viscosity);
  SLP_ob = ker.stokesSLmatrixInteractionWeightless(bb, wall, viscosity);

  DLP_ob = ker.stokesDLmatrixInteractionWeightless(bb,wall);
  DLPT_bo = ker.stokesDLTmatrixInteractionWeightless(wall,bb);
  
   
  % THEN SOLVE FOR BODY VELOCITY AND TRACTION, WALL TRACTION

  % form RHS
  RHS = zeros(2*NpointsWall + 2*Npoints + 3,1);
  RHS(1:2*Npoints) = 0;
  RHS(2*Npoints+1:2*Npoints+3) = -[ext_force;ext_torque];
  RHS(2*Npoints+4:end) = 0;
    
  % form the LHS matrix
  
  Hmat = diag(H);
  

  mat11 = SLP;
  mat12 =  -0.5*K-(DLP*Hmat)*K;
  mat13 = SLP_bo;

  mat21 = -0.5.*KT-(KT*Hmat)*DLPT;
  mat22 = zeros(3);
  mat23 = -(KT*Hmat)*DLPT_bo;
  
  mat31 = SLP_ob;
  mat32 = -(DLP_ob*Hmat)*K;
  mat33 = wallSLP;
   

  MAT = [mat11 mat12 mat13;...
         mat21 mat22 mat23;...
         mat31 mat32 mat33];
  
  if iflag
  % Check if SLP_bo is transpose(SLP_ob)
  check_norm = norm(MAT-MAT');
  disp(['Symmetric MAT check: ' num2str(check_norm)])
  end

  % solve the system
  sol = gmres(MAT,RHS,[],1E-12,50);

  % dissect the solution
  

  traction = sol(1:2*Npoints)./H;
  ui = sol(2*Npoints+1:2*Npoints+2);
  wi = sol(2*Npoints+3);
  U = [ui;wi];
  ubody_found = (K*U)./H;
  wall_traction = sol(2*Npoints+4:end)./Hwall;

    
%     wall_vel = wallSLP * sol(2*Npoints+4:end);
%     wallVelInt = sum((wall_vel(1:end/2).*wall.normal(1:end/2) + wall_vel(end/2+1:end).*wall.normal(end/2+1:end))*2*pi/wall.N.*wall.sa);
%     disp(['Net velocity on wall: ' num2str(wallVelInt)])

%   if iflag
%   figure(2); clf; hold on;
%   
%   quiver(bb.X(1:end/2),bb.X(end/2+1:end),ubody_found(1:end/2),ubody_found(end/2+1:end));
%   axis equal
%   quiver(bb.X(1:end/2),bb.X(end/2+1:end),traction(1:end/2),traction(end/2+1:end))
%   quiver(Xwalls(1:end/2), Xwalls(end/2+1:end),wall_traction(1:end/2),wall_traction(end/2+1:end))
%   plot(Xwalls(1:end/2), Xwalls(end/2+1:end),'linewidth',2)
%   plot(bb.X(1:end/2), bb.X(end/2+1:end),'linewidth',2);
%   legend('Body velocity', 'Body traction', 'Wall traction')
%   pause
%   end

  % update the position and orientation
  center = bb.center;
    
  X0 = [X(1:Npoints)-center(1); X(Npoints+1:end)-center(2)];
  center = center + ui*dt;

  xnew = center(1)+cos(wi*dt)*X0(1:Npoints)+sin(wi*dt)*X0(Npoints+1:end);
  ynew = center(2)-sin(wi*dt)*X0(1:Npoints)+cos(wi*dt)*X0(Npoints+1:end);
  X = [xnew; ynew];
  
  orientation = orientation + wi*dt;
  
  % Display the results  
  error = display_results(ui,wi,expected_velocity,time,X,Xwalls,radius);

end % end while time < time_horizon
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function error = display_results(ui, wi, expected_velocity, time, X, Xwalls,radius)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Time = ' num2str(time) 's'])
disp(['Velocity is found: ' num2str(reshape(ui,1,[])) 'um/s'])
disp(['Expected velocity is: ' num2str(expected_velocity) 'um/s'])
error = abs(ui(1)-expected_velocity)./expected_velocity;
disp(['Error is: ' num2str(error)])
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