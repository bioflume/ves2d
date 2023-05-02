function [ave_iter,max_iter,trajectories] = starSuspensionConfined_TestSymmAlpert(shearStrength, dt, N)
addpath ../src/

%% Problem setup
% Fluid properties
viscosity = 1; % Pa.s -- 1000x water's viscosity

% Confinement
wall_radius = 10;
Nwall = 256;

theta = [0:Nwall-1]'/Nwall * 2 * pi;
Xwalls = [wall_radius*cos(theta); wall_radius*sin(theta)];

% Body geometry
Npoints = N;
xcenters = [-4 -3 -1 0];
ycenters = [-1 -2 -4 -5];

% xcenters = [-5 -3 0 2];
% ycenters = [-1 -3 -6 -8];

angles = [pi/4 0 pi/4 0;
    0 pi/4 0 pi/4;
    pi/4 0 pi/4 0;
    0 pi/4 0 pi/4];

[xc,yc] = meshgrid(xcenters,ycenters);
folds = 4;
X = zeros(2*N,numel(xc(:)));
t = (0:N-1)'*2*pi/N;
radius = 1 + 0.3*cos(folds*t);
Xt = [radius.*cos(t);radius.*sin(t)];
dy = max(Xt(end/2+1:end))-min(Xt(end/2+1:end));
Xt = Xt/dy;
for ik = 1 : numel(xc(:))
% rotate Xt and locate the center
X(1:end/2,ik) = xc(ik) + cos(angles(ik))*Xt(1:N)+sin(angles(ik))*Xt(N+1:end);
X(end/2+1:end,ik) = yc(ik) - sin(angles(ik))*Xt(1:N)+cos(angles(ik))*Xt(N+1:end);
end

% Time scale
time_horizon = 17;

% Call kernels
ker = kernels(Npoints);


% Wall kernels
kerWall = kernels(Nwall);
wall = walls(Xwalls, [0 0], zeros(size(Xwalls)));
wallSLP = kerWall.stokesSLmatrixAlpertWeightless(wall,viscosity);
% Symmetric Alpert:
wallSLP = 0.5 * (wallSLP + wallSLP');
wall_velocity = [-Xwalls(end/2+1:end);Xwalls(1:end/2)];
%% Take time steps
time = 0;
nsteps = 0;
tot_iter = 0;
velocities = [];
trajectories = [];
max_iter = -inf;
while time < time_horizon
  % update time
  time = time + dt; 
  nsteps = nsteps + 1;
  
  % build body object
  bb = body(X, [], [], [], []);

  % calculate the K matrix
  K = bb.calc_K_matrix();
  KT = bb.calc_KT_matrix_weightless([]);
  Hmat = bb.calc_jacDiag_matrix();
  H = bb.calc_jac_matrix();

  % calculate the single and double layer matrices
  SLP = ker.stokesSLmatrixAlpertWeightless(bb,viscosity);  
  DLP = ker.stokesDLmatrixWeightless(bb);
  DLPT = ker.stokesDLTmatrixWeightless(bb);
  
   
  % THEN SOLVE FOR BODY VELOCITY AND TRACTION, WALL TRACTION
  RHS = zeros(2*bb.N*bb.nv + 3*bb.nv + 2*wall.N,1);
  RHS(2*bb.N*bb.nv+3*bb.nv+1:end) = wall_velocity;
  % form RHS
  for k = 1 : bb.nv
    SLP(:,:,k) = 0.5 * (SLP(:,:,k) + SLP(:,:,k)');
  end


  % solve the system
  sys_size = 2*bb.N*bb.nv+3*bb.nv+2*wall.N;
  [sol,iflag,~,iter,~] = gmres(@(X) ...
      TimeMatVec(X,bb,wall,ker,kerWall,wallSLP,SLP,DLP,DLPT,K,KT,Hmat,viscosity),...
      RHS,[],1E-4,sys_size,[]);
  
  tot_iter = tot_iter + iter(2);
  if iter(2) > max_iter; max_iter = iter(2); end;
  disp(['Number of GMRES iterations is ' num2str(iter(2))])
  if nsteps == 1; pause; end;

  % dissect the solution
  velocity = zeros(3,bb.nv);
  centers = zeros(2,bb.nv);
  traction = zeros(2*bb.N,bb.nv);

  for k = 1 : bb.nv
    istart = (k-1)*(2*bb.N+3)+1;
    iend = istart+2*bb.N-1;
    traction(:,k) = sol(istart:iend)./H(:,k);
    istart = iend+1;
    iend = istart+2;
    velocity(:,k) = sol(istart:iend);
    
    cx = mean(bb.X(1:bb.N,k));
    cy = mean(bb.X(bb.N+1:2*bb.N,k));

    X0 = [bb.X(1:bb.N,k)-cx; bb.X(bb.N+1:2*bb.N,k)-cy];
    cx = cx + velocity(1,k)*dt;
    cy = cy + velocity(2,k)*dt;
    centers(:,k) = [cx;cy];
    
    xnew = cx+cos(velocity(3,k)*dt)*X0(1:bb.N)+sin(velocity(3,k)*dt)*X0(bb.N+1:end);
    ynew = cy-sin(velocity(3,k)*dt)*X0(1:bb.N)+cos(velocity(3,k)*dt)*X0(bb.N+1:end);
    X(:,k) = [xnew;ynew];

  end
  velocities = [velocities velocity(:)];
  trajectories = [trajectories centers(:)];

  
  % Display the results  
  display_results(time,X,Xwalls,traction,trajectories,velocity);
end % end while time < time_horizon
ave_iter = tot_iter / nsteps;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = TimeMatVec(sol, bb, wall, ker, kerWall, wallSLP, SLP, DLP, DLPT, K, KT, H, viscosity)

val = zeros(size(sol));
N = bb.N;
nv = bb.nv;

valTraction = zeros(2*N,nv);
valVelocity = zeros(3,nv);

traction = zeros(2*N,nv);
velocity = zeros(3,nv);

for k = 1 : nv
  istart = (k-1)*(2*N+3)+1;
  iend = istart+2*N-1;
  traction(:,k) = sol(istart:iend);
  istart = iend+1;
  iend = istart+2;
  velocity(:,k) = sol(istart:iend);
end

% traction on the wall
wallTraction = sol(iend+1:end);

% Calculate the self-interactions
SLvals = zeros(2*N,nv);
KDLvals = zeros(2*N,nv); % includes -1/2 * K - DHK
KTDTvals = zeros(3,nv); % includes -1/2*KT - KT *H *DT
surfVels = zeros(2*N,nv);
for k = 1 : nv
  SLvals(:,k) = SLP(:,:,k) * traction(:,k);
  KDLvals(:,k) = -1/2*K(:,:,k)*velocity(:,k)-DLP(:,:,k)*(H(:,:,k)*(K(:,:,k)*velocity(:,k)));
  KTDTvals(:,k) = -1/2*KT(:,:,k)*traction(:,k) - ((KT(:,:,k)*H(:,:,k))*DLPT(:,:,k))*traction(:,k);
  surfVels(:,k) = H(:,:,k)* (K(:,:,k) * velocity(:,k));
end


% Calculate the hydrodynamic interactions
SLinteract = ker.stokesSL_times_density(bb,bb,viscosity,traction,1);
DLinteract = ker.stokesDL_times_density(bb,bb,surfVels,1);
DLTinteract = ker.stokesDLT_times_density(bb,bb,traction,1);

% Calculate interactions on the wall
SLonWall = ker.stokesSL_times_density(bb,wall,viscosity,traction,0);
DLonWall = ker.stokesDL_times_density(bb,wall,surfVels, 0);

% Solve the wall equation
SLwallVal = wallSLP * wallTraction;
valWall = SLwallVal + SLonWall - DLonWall;

% Calculate interactions from the wall on the bodies
SLwall2body = kerWall.stokesSL_times_density(wall,bb,viscosity,wallTraction,0);
DLTwall2body = kerWall.stokesDLT_times_density(wall,bb,wallTraction,0);
SLinteract = SLinteract + SLwall2body;
DLTinteract = DLTinteract + DLTwall2body;

valTraction = valTraction + SLvals + KDLvals + SLinteract - DLinteract;
valVelocity = valVelocity + KTDTvals;
for k = 1 : nv
  valVelocity(:,k) = valVelocity(:,k) - (KT(:,:,k)*H(:,:,k))*DLTinteract(:,k);
  istart = (k-1)*(2*N+3)+1;
  iend = istart - 1 + 2*N + 3;
  val(istart:iend) = [valTraction(:,k);valVelocity(:,k)];
end
val(iend+1:end) = valWall;

end % val = TimeMatVec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_results(time,X,Xwalls,traction,trajectories,velocity)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Time = ' num2str(time) 's'])

if 1
figure(1); clf;
vecx = [interpft(X(1:end/2,:),128); X(1,:)];
vecy = [interpft(X(end/2+1:end,:),128); X(end/2+1,:)];

plot(vecx, vecy, 'k');
hold on;
fill(vecx, vecy,'k');
axis equal

vecx = [Xwalls(1:end/2); Xwalls(1)];
vecy = [Xwalls(end/2+1:end); Xwalls(end/2+1)];

plot(vecx,vecy,'k');

% for k = 1 : size(X,2)
%   istart = (k-1)*2+1;
%   plot(trajectories(istart,:),trajectories(istart+1,:),'linewidth',2)
%   quiver(trajectories(istart,end),trajectories(istart+1,end),1*velocity(1,k),1*velocity(2,k),'c','LineWidth',2)
%   quiver(X(1:end/2,k),X(end/2+1:end,k),traction(1:end/2,k),traction(end/2+1:end,k),'r')
% end

grid off
box on
title(['Time = ' num2str(time) 's'])

% xlim([-7 3])
% ylim([-8 2])

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