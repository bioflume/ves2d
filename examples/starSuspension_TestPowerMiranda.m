function [ave_iter,max_iter,trajectories] = starSuspension_TestPowerMiranda(shearStrength, dt, N)


addpath ../src/

%% Problem setup
% Fluid properties
viscosity = 1; % Pa.s -- 1000x water's viscosity

% Body geometry
Npoints = N;
bgFlow = @(y) [shearStrength*y;zeros(size(y))];
xcenters = [-4 -3 -1 0];
ycenters = [-1 -2 -4 -5];

% xcenters = [-5 -3 0 2];
% ycenters = [-1 -3 -6 -8];

angles = [pi/4 0 pi/4 0;
    0 pi/4 0 pi/4;
    pi/4 0 pi/4 0;
    0 pi/4 0 pi/4];

% angles = [pi/4 pi/4];

[xc,yc] = meshgrid(xcenters,ycenters);
folds = 4;
X = zeros(2*N,numel(xc(:)));
t = (0:N-1)'*2*pi/N;
radius = 1 + 0.3*cos(folds*t);
Xt = [radius.*cos(t);radius.*sin(t)];
dy = max(Xt(end/2+1:end))-min(Xt(end/2+1:end));
Xt = Xt/dy;
nrigid = numel(xc);
for ik = 1 : numel(xc(:))
% rotate Xt and locate the center
X(1:end/2,ik) = xc(ik) + cos(angles(ik))*Xt(1:N)+sin(angles(ik))*Xt(N+1:end);
X(end/2+1:end,ik) = yc(ik) - sin(angles(ik))*Xt(1:N)+cos(angles(ik))*Xt(N+1:end);
end

% Time scale
time_horizon = 5*dt;

% Call kernels
ker = kernels(Npoints);

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
  DLP = ker.stokesDLmatrixWeightless(bb);  
   
  % THEN SOLVE FOR BODY VELOCITY AND TRACTION, WALL TRACTION
  RHS = zeros(2*bb.N*bb.nv + 3*bb.nv,1);
  % form RHS
  for k = 1 : bb.nv
    istart = (k-1)*(2*bb.N+3)+1;
    iend = istart + 2*bb.N - 1;
    RHS(istart:iend) = -bgFlow(X(end/2+1:end,k));
  end

  % build the preconditioner
  bdiagL = zeros(2*Npoints+3,2*Npoints+3,nrigid);
  bdiagU = zeros(2*Npoints+3,2*Npoints+3,nrigid);
  for k=1:nrigid
    invHmat = diag(1./H(:,k));
    [bdiagL(:,:,k),bdiagU(:,:,k)] = lu([0.5*invHmat-DLP(:,:,k) K(:,:,k); KT(:,:,k) zeros(3)]);
  end
  % solve the system
  sys_size = 2*bb.N*bb.nv+3*bb.nv;

  if 0
  [sol,iflag,~,iter,~] = gmres(@(X) ...
      TimeMatVec(X,bb,ker,DLP,K,KT,H,viscosity),...
      RHS,[],1E-10,sys_size,[]);
  else
  [sol,iflag,~,iter,~] = gmres(@(X) ...
      TimeMatVec(X,bb,ker,DLP,K,KT,H,viscosity),...
      RHS,[],1E-10,sys_size, @(z) precondBD(z,bdiagL,bdiagU));
  end
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
  display_results(time,X,traction,trajectories,velocity);
end % end while time < time_horizon
ave_iter = tot_iter / nsteps;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = precondBD(z,bdiagL,bdiagU)
Nrd = (size(bdiagL,1)-3)/2;
nvrd = size(bdiagL,3);

valRigid = zeros(2*Nrd+3,nvrd);
for k = 1 : nvrd
  valRigid(:,k) = bdiagU(:,:,k)\(bdiagL(:,:,k)\z((k-1)*(2*Nrd+3)+1:k*(2*Nrd+3)));
end
val = valRigid(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = TimeMatVec(sol, bb, ker, DLP, K, KT, H, viscosity)

val = zeros(size(sol));
N = bb.N;
nv = bb.nv;

valDensity = zeros(2*N,nv);
valVelocity = zeros(3,nv);

density = zeros(2*N,nv);
velocity = zeros(3,nv);

for k = 1 : nv
  istart = (k-1)*(2*N+3)+1;
  iend = istart+2*N-1;
  density(:,k) = sol(istart:iend);
  istart = iend+1;
  iend = istart+2;
  velocity(:,k) = sol(istart:iend);
end

% Calculate the self-interactions
DLvals = zeros(2*N,nv); % includes 1/2 - D*density
GRKvals = zeros(2*N,nv);
KTvals = zeros(3,nv);
stokeslet = zeros(2,nv);
rotlet = zeros(1,nv);
surfVels = zeros(2*N,nv);
for k = 1 : nv
  
  DLvals(:,k) = +(1./(2*H(:,k))).*density(:,k)-DLP(:,:,k)*density(:,k);
  
  stokeRot = KT(:,:,k)*density(:,k);
  KTvals(:,k) = stokeRot;
  stokeslet(:,k) = stokeRot(1:2);
  rotlet(:,k) = stokeRot(3);
  
  GRKvals(:,k) = ker.stokeRotTimes_SelfForceTorque(bb.X(:,k),stokeslet(:,k),rotlet(:,k),viscosity);
  surfVels(:,k) = (K(:,:,k) * velocity(:,k));
end


% Calculate the hydrodynamic interactions

DLinteract = ker.stokesDL_times_density(bb,bb,density,1);
GRinteract = ker.stokeRotTimes_forceTorque(bb,bb,stokeslet,rotlet,viscosity,1);
valDensity = valDensity + DLvals + 0*GRKvals - DLinteract + 0*GRinteract - surfVels;
valVelocity = valVelocity + KTvals;
for k = 1 : nv
  istart = (k-1)*(2*N+3)+1;
  iend = istart - 1 + 2*N + 3;
  val(istart:iend) = [valDensity(:,k);valVelocity(:,k)];
  
end

end % val = TimeMatVec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_results(time,X,traction,trajectories,velocity)

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

% for k = 1 : size(X,2)
%   istart = (k-1)*2+1;
%   plot(trajectories(istart,:),trajectories(istart+1,:),'linewidth',2)
%   quiver(trajectories(istart,end),trajectories(istart+1,end),1*velocity(1,k),1*velocity(2,k),'c','LineWidth',2)
%   quiver(X(1:end/2,k),X(end/2+1:end,k),traction(1:end/2,k),traction(end/2+1:end,k),'r')
% end

grid off
box on
title(['Time = ' num2str(time) 's'])

xlim([-7 3])
ylim([-8 2])

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