function [error,ave_iter,orientation,theta_theo] = ellipsoid_shearTestSDDj(shearStrength,N)

addpath ../src/
iflag = 1;

oc = curve;

%% Problem setup
% Fluid properties
viscosity = 1; % Pa.s -- 1000x water's viscosity

% Body geometry
radius_x = 1; % micro meters
radius_y = 10;
Npoints = N;
radius = max([radius_x;radius_y]);
bgFlow = @(y) [shearStrength*y;zeros(size(y))];

% Body shape discretization assuming a sphere
theta = [0:Npoints-1]'/Npoints * 2 * pi;
body_x = radius_x*cos(theta);
body_y = radius_y*sin(theta);

% X = [body_x;body_y];
% figure(1); clf;
% vecx = [interpft(X(1:end/2),128); X(1)];
% vecy = [interpft(X(end/2+1:end),128); X(end/2+1)];
% 
% plot(vecx, vecy, 'k');
% hold on;
% fill(vecx, vecy,'k');
% axis equal
% 
% 
% grid off
% box on
% pause

xlim([-2*radius 2*radius])
ylim([-2*radius 2*radius])

expected_angVelocity = @(theta) shearStrength/(radius_x^2 + radius_y^2) ...
    * (radius_x^2*sin(theta)^2 + radius_y^2*cos(theta)^2);  


% Time scale
time_horizon = 20/shearStrength;
dt = 1E-1/shearStrength; 
time_horizon = 5*dt;
% time_horizon = dt;

% Call kernels
ker = kernels(Npoints);

%% Take time steps
time = 0;
nsteps = 0;
X = [body_x; body_y];
center = [0; 0];
orientation = 0; % angle
traction = zeros(2*Npoints,1);

tot_iter = 0;
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
  Hmat = diag(H);
  HmatInv = diag(1./H);
  
  GR = ker.stokeRotmatrix(bb,bb,viscosity);

  % THEN SOLVE FOR BODY VELOCITY AND TRACTION, WALL TRACTION

  % form RHS
  RHS = zeros(2*Npoints + 3,1);
  RHS(1:2*Npoints) = -bgFlow(X(end/2+1:end));
    
  % form the LHS matrix
%   What we are suggesting
  if 0
  mat11 = SLP;
  mat12 =  -0.5*K-(DLP*Hmat)*K;

  mat21 = -0.5.*KT-(KT*Hmat)*DLPT;
  mat22 = zeros(3);
  end

  % Shravan
  if 0
  mat11 = HmatInv/2 - DLP + SLP;
  mat12 = -K;
  
  mat21 = -KT*Hmat*DLPT;
  mat22 = zeros(3);
  end

  % Power and Miranda
  if 1
  mat11 = HmatInv/2 - DLP + GR*KT;
  mat12 = -K;

  mat21 = -KT;
  mat22 = zeros(3);
  end

  MAT = [mat11 mat12;...
         mat21 mat22];
  
  if iflag
  % Check if SLP_bo is transpose(SLP_ob)
  check_norm = norm(MAT-MAT');
  disp(['Symmetric MAT check: ' num2str(check_norm)])
  end

  % solve the system
  [sol,~,~,iter] = gmres(MAT,RHS,[],1E-10,size(MAT,1));
  tot_iter = tot_iter + iter(2);
  disp(['Number of GMRES iterations is ' num2str(iter(2))])

  % dissect the solution
  

  traction = sol(1:2*Npoints)./H;
  ui = sol(2*Npoints+1:2*Npoints+2);
  wi = sol(2*Npoints+3);
  U = [ui;wi];

  % update the position and orientation
  cx = mean(X(1:Npoints));
  cy = mean(X(Npoints+1:2*Npoints));
    
  X0 = [X(1:Npoints)-cx; X(Npoints+1:end)-cy];
  cx = cx + ui(1)*dt;
  cy = cy + ui(2)*dt;

  xnew = cx+cos(wi*dt)*X0(1:Npoints)+sin(wi*dt)*X0(Npoints+1:end);
  ynew = cy-sin(wi*dt)*X0(1:Npoints)+cos(wi*dt)*X0(Npoints+1:end);
  X = [xnew; ynew];
  
  expected_velocity = expected_angVelocity(orientation);

  orientation = orientation + wi*dt;
  
  
  % Display the results  
  error(nsteps,1) = display_results(ui,wi,expected_velocity,time,X,radius);
end % end while time < time_horizon
ave_iter = tot_iter / nsteps;

theta_theo = atan(radius_y/radius_x * tan(-radius_x*radius_y*shearStrength*time/(radius_x^2 + radius_y^2)));
orientation = oc.getIncAngle(X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function error = display_results(ui, wi, expected_velocity, time, X,radius)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Time = ' num2str(time) 's'])
disp(['Angular velocity is found: ' num2str(wi)])
disp(['Expected velocity is: ' num2str(expected_velocity) 'um/s'])
error = abs(wi-expected_velocity)./abs(expected_velocity);
disp(['Error is: ' num2str(error)])
disp(['Velocity is found: ' num2str(reshape(ui,1,[])) 'um/s'])

if 1
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

xlim([-2*radius 2*radius])
ylim([-2*radius 2*radius])

pause(0.1)
end
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