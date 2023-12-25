classdef capsules < handle
% This class implements standard calculations that need to
% be done to a vesicle, solid wall, or a collection of arbitrary
% target points (such as tracers or pressure/stress targets)
% Given a vesicle, the main tasks that can be performed are
% computing the required derivatives (bending, tension, surface
% divergence), the traction jump, the pressure and stress, 
% and constructing structures required for near-singluar
% integration

properties
N;          % number of points per vesicle
nv;         % number of vesicles
X;          % positions of vesicles
sig;        % tension of vesicles
u;          % velocity field of vesicles
kappa;      % bending modulus
viscCont;   % viscosity contrast
xt;         % tangent unit vector
sa;         % Jacobian
isa;        % inverse of Jacobian
length;     % minimum length over allc vesicles
cur;        % curvature
center;     % center of the point required for stokeslets
            % and rotlets
IK;         % index of Fourier modes for fft and ifft
            % that are needed repetatively
antiAlias;  % flag to use anti-aliasing
uprate;     % the rate of upsampling if doing anti-aliasing

end %properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = capsules(X,sigma,u,kappa,viscCont,antiAlias)
% capsules(X,sigma,u,kappa,viscCont) sets parameters and options for
% the class; no computation takes place here.  
%
% sigma and u are not needed and typically unknown, so just set them to
% an empty array.

o.N = size(X,1)/2;              % points per vesicle
o.nv = size(X,2);               % number of vesicles
o.X = X;                        % position of vesicle
oc = curve;
% Jacobian, tangent, and curvature
[o.sa,o.xt,o.cur] = oc.diffProp(o.X);
o.isa = 1./o.sa;
o.sig = sigma;          % Tension of vesicle
o.u = u;                % Velocity of vesicle
o.kappa = kappa;        % Bending modulus
o.viscCont = viscCont;  % Viscosity contrast
% center of vesicle.  Required for center of rotlets and
% stokeslets in confined flows
o.center = [mean(X(1:o.N,:));mean(X(o.N+1:2*o.N,:))];

% minimum arclength needed for near-singular integration
[~,~,len] = oc.geomProp(X);
o.length = min(len);

% ordering of the fourier modes.  It is faster to compute once here and
% pass it around to the fft differentitation routine
o.IK = fft1.modes(o.N,o.nv);

% default is to do no upsampling.  setUpRate can be called to find an
% upsampling rate so that the aliasing error is controlled
o.antiAlias = antiAlias;
o.uprate = 1;


end % capsules

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setUpRate(vesicle,op)
% setUpRate(vesicle,op) picks the upsampling rate for each vesicle.
% This is done by looking at a maximum sized upsampled geometry, and
% then choosing a minimal half integer so that the high frequencies and
% low freqeuncies of the single-layer potential times bending above a
% desired ratio.

maxUp = 16;
N = vesicle.N;
nv = vesicle.nv;

Nup = maxUp*N;
modesUp = (-Nup/2:Nup/2-1)';

oc = curve;
[x,y] = oc.getXY(vesicle.X);
x = interpft(x,Nup);
y = interpft(y,Nup);

vesicleUp = capsules([x;y],[],[],1,1,false);

Gf = vesicleUp.tracJump(vesicleUp.X,zeros(Nup,vesicle.nv));
%Gf = op.exactStokesSLdiag(vesicleUp,[],Gf);

uprates = vesicle.uprate*ones(vesicle.nv,1);
for k = 1:nv
  z = Gf(1:Nup,k) + 1i*Gf(Nup+1:end,k);
  zh = fftshift(fft(z))/Nup;

  alpha = 1.5;
  uprates(k) = alpha;
  lowEnergy = sqrt(sum(abs(zh(abs(modesUp) < alpha/4*N)).^2));
  highEnergy = sqrt(sum(abs(zh(abs(modesUp) < alpha/2*N & ...
      abs(modesUp) >= alpha/4*N)).^2));


  % no more than upsampling by maxUp
  desiredRatio = 1.0e0;
  while (highEnergy/lowEnergy > desiredRatio & alpha < maxUp)
    alpha = alpha + 5e-1;
    % upsample by the halves
    uprates(k) = alpha;
    lowEnergy = sqrt(sum(abs(zh(abs(modesUp) < alpha/4*N)).^2));
    highEnergy = sqrt(sum(abs(zh(abs(modesUp) < alpha/2*N & ...
        abs(modesUp) >= alpha/4*N)).^2));
  end
end
vesicle.uprate = max(uprates);

% disp(['Vesicle uprate: ' num2str(vesicle.uprate)])
%highEnergy/lowEnergy

%figure(3); clf
%s = find(abs(modesUp) < alpha/4*N);
%semilogy(modesUp(s),abs(zh(s)),'r.')
%hold on
%s = find(abs(modesUp) >= alpha/4*N);
%semilogy(modesUp(s),abs(zh(s)),'g.')
%xlim([-alpha/2*N alpha/2*N])

end % setUpRate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = tracJump(o,f,sigma)
% tracJump(f,sigma) computes the traction jump where the derivatives
% are taken with respect to a linear combiation of previous time steps
% which is stored in object o Xm is 2*N x nv and sigma is N x nv

if ~o.antiAlias
  f = o.bendingTerm(f) + o.tensionTerm(sigma);
else
  oc = curve;
  [x,y] = oc.getXY(o.X);
  [fx,fy] = oc.getXY(f);
  Nup = o.uprate*o.N;
  xup = interpft(x,Nup); yup = interpft(y,Nup);
  fxup = interpft(fx,Nup); fyup = interpft(fy,Nup);
  sigma = interpft(sigma,Nup);
  
  % compute upsampled vesicle
  vesicleUp = capsules([xup;yup],[],[],o.kappa,1,true);
  
  fup = vesicleUp.bendingTerm([fxup;fyup]) + ...
      vesicleUp.tensionTerm(sigma);
  
  % Get low frequencies of forces
  [fupx,fupy] = oc.getXY(fup); 
  fxh = fft(fupx)/o.uprate;
  fyh = fft(fupy)/o.uprate;
  fx = real(ifft([fxh(1:o.N/2,:);fxh(Nup - o.N/2+1:Nup,:)]));
  fy = real(ifft([fyh(1:o.N/2,:);fyh(Nup - o.N/2+1:Nup,:)]));
  f = oc.setXY(fx,fy);
end


end % tracJump

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ben = bendingTerm(o,f)
% ben = bendingTerm(f) computes the term due to bending
% -kappa*fourth-order derivative

ben = [-o.kappa*curve.arcDeriv(f(1:o.N,:),4,o.isa,o.IK);...
  -o.kappa*curve.arcDeriv(f(o.N+1:2*o.N,:),4,o.isa,o.IK)];

end % bendingTerm
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ten = tensionTerm(o,sig)
% ten = tensionTerm(o,sig) computes the term due to tension (\sigma *
% x_{s})_{s}

ten = [curve.arcDeriv(sig.*o.xt(1:o.N,:),1,o.isa,o.IK);...
    curve.arcDeriv(sig.*o.xt(o.N+1:2*o.N,:),1,o.isa,o.IK)];

end % tensionTerm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = surfaceDiv(o,f)
% divf = surfaceDiv(f) computes the surface divergence of f with respect
% to the vesicle stored in object o.  f has size N x nv

oc = curve; 
[fx,fy] = oc.getXY(f);

if ~o.antiAlias
  [tangx,tangy] = oc.getXY(o.xt);
  f = curve.arcDeriv(fx,1,o.isa,o.IK).*tangx + ...
    curve.arcDeriv(fy,1,o.isa,o.IK).*tangy;
else
  Nup = o.uprate*o.N;
  [x,y] = oc.getXY(o.X);
  xup = interpft(x,Nup); yup = interpft(y,Nup);
  fxup = interpft(fx,Nup); fyup = interpft(fy,Nup);
  vesicleUp = capsules([xup;yup],[],[],o.kappa,1,true);
  [tangx,tangy] = oc.getXY(vesicleUp.xt);
  f = curve.arcDeriv(fxup,1,vesicleUp.isa,vesicleUp.IK).*tangx + ...
    curve.arcDeriv(fyup,1,vesicleUp.isa,vesicleUp.IK).*tangy;

  f = fft(f)/o.uprate;
  f = real(ifft([f(1:o.N/2,:);f(Nup - o.N/2+1:Nup,:)]));
end

end % surfaceDiv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function adh = adhesionTerm(o,Wconst,d0)
% adh = adhsionTerm(Wconst,d0) computes the term due to adhesion.
% Wconst is the strength of the adhesion force and d0 corresponds to
% the distance where the adhesive force is the strongest

adh = zeros(2*o.N,o.nv);
for k = 1:o.nv
  notk = [(1:k-1) (k+1:o.nv)];
  adhesive_x = o.X(1:o.N,notk);
  adhesive_y = o.X(o.N+1:2*o.N,notk);
  ds = o.sa(:,notk);
  
  for j = 1:o.N
    dist2 = (o.X(j,k) - adhesive_x).^2 + ...
        (o.X(j+o.N,k) - adhesive_y).^2;
    Dw = -Wconst * 4*d0^2.*(d0^2 - dist2)./dist2.^(2.5);
    adh(j,k) = sum(sum(-Dw .* (o.X(j,k) - adhesive_x).*ds));
    adh(j+o.N,k) = sum(sum(-Dw .* (o.X(j+o.N,k) - adhesive_y).*ds));
  end
end
end % adhesionTerm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = repulStrengthScale(o,minDist,tt,uWall)
% W = repulStrengthScale(o,minDist,tt,uWall) automatically 
% determines the strength of the repulsion. It considers a single 
% point on the first vesicle.
% Then assumes that there is another point in the close proximity 
% of that point. It computes the repulsion force and SLP of that 
% force. The strength is the ratio of the maximum velocity on 
% the inlet to the maximum SLP.

op = tt.op;

oc = curve;
[x,y] = oc.getXY(o.X);

[~,arcLen] = oc.arcLengthParameter(x,y);
hmax = max(diff(arcLen));

eta = minDist*hmax;
L = 1;
dist = eta/1.2;

dF = -L*(2*L+1)*(L+1)/3 + L*(L+1)*1.2;

repx = dF .* dist;

rep = zeros(o.N*2,1);
rep(1) = repx;

if~isempty(tt.Galpert)
  uRep = tt.Galpert(:,:,1)*rep;
else
  Galpert = op.stokesSLmatrix(o);
  uRep = Galpert*rep;
end

W = max(1*uWall./(sqrt(uRep(1:end/2).^2+uRep(end/2+1:end).^2)));
end % repulStrengthScale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,eta] = repulStrengthScaleSimple(o,minDist,tt,uWall)
% W = repulStrengthScale(o,minDist,tt,uWall) automatically 
% determines the strength of the repulsion. It considers a single 
% point on the first vesicle.
% Then assumes that there is another point in the close proximity 
% of that point. It computes the repulsion force and SLP of that 
% force. The strength is the ratio of the maximum velocity on 
% the inlet to the maximum SLP.

op = tt.op;

oc = curve;
[x,y] = oc.getXY(o.X);

[~,arcLen] = oc.arcLengthParameter(x,y);
hmax = max(diff(arcLen));

eta = minDist*hmax;
L = 1;
dist = eta/1.2;

dF = -L*(2*L+1)*(L+1)/3 + L*(L+1)*1.2;

repx = dF .* dist;

rep = zeros(o.N*2,1);
rep(1) = repx;

if~isempty(tt.Galpert)
  uRep = tt.Galpert(:,:,1)*rep;
else
  Galpert = op.stokesSLmatrix(o);
  uRep = Galpert*rep;
end

W = 1*max(1*uWall./(sqrt(uRep(1:end/2).^2+uRep(end/2+1:end).^2)));
end % repulStrengthScaleSimple
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rep = repulsionScheme(o,X,W,minDist,walls,wallsInt,wallsExt)
% rep = repulsionScheme(o,W,minDist,walls,wallsInt,wallsExt) computes the
% repulsion force. if the distance between points is greater than
% minDist/2*hmax then there is no need for repulsion. 
% Repulsion range is minDist of the maximum grid size (h).
% Repulsion is computed using the discrete penalty layers given in Grinspun
% et al. (2009), Asynchronuous Contact Mechanics.

% flag whether we have walls or not
confined = 0;
if ~isempty(walls) || ~isempty(wallsInt)
  confined = 1;
end
% flag for whether there are two sets of walls with different
% discretization
diffDiscWalls = 0;
if ~isempty(wallsInt)
  diffDiscWalls = 1;
end

oc = curve;

% Compute x,y coordinates of net repulsive force on each point of each
% vesicle due to all other vesicles and walls
rep = zeros(2*o.N,o.nv);

% Upsample if necessary
Nup = max(2*o.N,96);

[ox,oy] = oc.getXY(o.X);
oxup = interpft(ox,Nup); oyup = interpft(oy,Nup);

% X might be o.X or deltaX (when we do SDC correction)
[x,y] = oc.getXY(X);
xup = interpft(x,Nup); yup = interpft(y,Nup);

if confined
  if ~diffDiscWalls 
    [xwalls,ywalls] = oc.getXY(walls.X);
    Nbd = max(2*walls.N,256);
    xwalls = interpft(xwalls,Nbd);
    ywalls = interpft(ywalls,Nbd);
  else
    [xwallsInt,ywallsInt] = oc.getXY(wallsInt.X);
    NbdInt = max(wallsInt.N,256);
    xwallsInt = interpft(xwallsInt,NbdInt);
    ywallsInt = interpft(ywallsInt,NbdInt);
    
    [xwallsExt,ywallsExt] = oc.getXY(wallsExt.X);
    NbdExt = max(wallsExt.N,256);
    xwallsExt = interpft(xwallsExt,NbdExt);
    ywallsExt = interpft(ywallsExt,NbdExt);
  end
end

for k = 1:o.nv  
  repx = zeros(o.N,1); repy = zeros(o.N,1);  
  [~,arcLen] = oc.arcLengthParameter(ox(:,k),oy(:,k));
  hmax = max(diff(arcLen));
  
  notk = [(1:k-1) (k+1:o.nv)];
  notk_ox = oxup(:,notk);
  notk_oy = oyup(:,notk);
  
  notk_x = xup(:,notk);
  notk_y = yup(:,notk);
  
  eta = minDist * hmax;
  
  for j = 1 : o.N
    if o.nv ~= 1  
      % Find the distances to each point on a vesicle  
      dist = ((ox(j,k) - notk_ox).^2 + ...
          (oy(j,k) - notk_oy).^2).^0.5;
      % Find out the maximum amount of layers necessary for each distance
      L = floor(eta./dist);

      % Stiffness l^2
      dF = -L.*(2*L+1).*(L+1)/3 + L.*(L+1).*eta./dist;
      
      repx(j) = sum(sum(dF .* (x(j,k) - notk_x)))*2*pi/Nup;
      repy(j) = sum(sum(dF .* (y(j,k) - notk_y)))*2*pi/Nup;
    end
    if confined
      etaWalls = hmax*minDist;  
      if ~diffDiscWalls  
        dist = ((ox(j,k) - xwalls).^2 + (oy(j,k) - ywalls).^2).^0.5;
        L = floor(etaWalls./dist);
      
        % Stiffness l^2
        dF = -L.*(2*L+1).*(L+1)/3 + L.*(L+1).*etaWalls./dist;
      
        wallRepx = sum(sum(dF .* (x(j,k)-xwalls)))*2*pi/Nbd;
        wallRepy = sum(sum(dF .* (y(j,k)-ywalls)))*2*pi/Nbd;
      
        % Repulsion due to wall has to beat repulsion due to vesicles if a
        % vesicle is a sandwich between a vesicle and a wall
        repx(j) = repx(j) + wallRepx;
        repy(j) = repy(j) + wallRepy;
      else
        % Interior walls    
        dist = ((ox(j,k) - xwallsInt).^2 + (oy(j,k) - ywallsInt).^2).^0.5;
        L = floor(etaWalls./dist);
      
        % Stiffness l^2
        dF = -L.*(2*L+1).*(L+1)/3 + L.*(L+1).*etaWalls./dist;
      
        wallIntRepx = sum(sum(dF .* (x(j,k)-xwallsInt)));
        wallIntRepy = sum(sum(dF .* (y(j,k)-ywallsInt)));
        
        % Exterior wall
        dist = ((ox(j,k) - xwallsExt).^2 + (oy(j,k) - ywallsExt).^2).^0.5;
        L = floor(etaWalls./dist);
      
        % Stiffness l^2
        dF = -L.*(2*L+1).*(L+1)/3 + L.*(L+1).*etaWalls./dist;
      
        wallExtRepx = sum(sum(dF .* (x(j,k)-xwallsExt)));
        wallExtRepy = sum(sum(dF .* (y(j,k)-ywallsExt)));
      
        % Repulsion due to wall has to beat repulsion due to vesicles if a
        % vesicle is a sandwich between a vesicle and a wall
        repx(j) = repx(j) + wallIntRepx + wallExtRepx;
        repy(j) = repy(j) + wallIntRepy + wallExtRepy; 
          
      end % diffDiscWalls
    end % confined
  end % o.N
  % repulsion on the kth vesicle multiplied with strength W
  rep(:,k) = W*[repx;repy];
end % o.nv
end % repulsionScheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rep = repulsionSchemeSimple(o,X,W,eta,walls,wallsInt,wallsExt)
% rep = repulsionScheme(o,W,minDist,walls,wallsInt,wallsExt) computes the
% repulsion force. if the distance between points is greater than
% minDist/2*hmax then there is no need for repulsion. 
% Repulsion range is minDist of the maximum grid size (h).
% Repulsion is computed using the discrete penalty layers given in Grinspun
% et al. (2009), Asynchronuous Contact Mechanics.

% flag whether we have walls or not
confined = 0;
if ~isempty(walls) || ~isempty(wallsInt)
  confined = 1;
end
% flag for whether there are two sets of walls with different
% discretization
diffDiscWalls = 0;
if ~isempty(wallsInt)
  diffDiscWalls = 1;
end

oc = curve;
nv = numel(X(1,:));
N = numel(X(:,1))/2;

% Compute x,y coordinates of net repulsive force on each point of each
% vesicle due to all other vesicles and walls
rep = zeros(2*N,nv);

[ox,oy] = oc.getXY(X);

if confined
  if ~diffDiscWalls 
    [xwalls,ywalls] = oc.getXY(walls.X);
  else
    [xwallsInt,ywallsInt] = oc.getXY(wallsInt.X);
    [xwallsExt,ywallsExt] = oc.getXY(wallsExt.X);
  end
end

for k = 1:nv
  repx = zeros(N,1); repy = zeros(N,1);  
  
  notk = [(1:k-1) (k+1:nv)];
  notk_ox = ox(:,notk);
  notk_oy = oy(:,notk);
  
  for j = 1 : N
    if nv ~= 1  
      % Find the distances to each point on a vesicle  
      dist = ((ox(j,k) - notk_ox).^2 + ...
          (oy(j,k) - notk_oy).^2).^0.5;
      % Find out the maximum amount of layers necessary for each distance
      L = floor(eta./dist);

      % Stiffness l^2
      dF = -L.*(2*L+1).*(L+1)/3 + L.*(L+1).*eta./dist;
      
      repx(j) = sum(sum(dF .* (ox(j,k) - notk_ox)));
      repy(j) = sum(sum(dF .* (oy(j,k) - notk_oy)));
    end
    if confined
      etaWalls = 2*eta;  
      if ~diffDiscWalls  
        dist = ((ox(j,k) - xwalls).^2 + (oy(j,k) - ywalls).^2).^0.5;
        L = floor(etaWalls./dist);
      
        % Stiffness l^2
        dF = -L.*(2*L+1).*(L+1)/3 + L.*(L+1).*etaWalls./dist;
      
        wallRepx = sum(sum(dF .* (ox(j,k)-xwalls)));
        wallRepy = sum(sum(dF .* (oy(j,k)-ywalls)));
      
        % Repulsion due to wall has to beat repulsion due to vesicles if a
        % vesicle is a sandwich between a vesicle and a wall
        repx(j) = repx(j) + wallRepx;
        repy(j) = repy(j) + wallRepy;
      else
        % Interior walls    
        dist = ((ox(j,k) - xwallsInt).^2 + (oy(j,k) - ywallsInt).^2).^0.5;
        L = floor(etaWalls./dist);
      
        % Stiffness l^2
        dF = -L.*(2*L+1).*(L+1)/3 + L.*(L+1).*etaWalls./dist;
      
        wallIntRepx = sum(sum(dF .* (ox(j,k)-xwallsInt)));
        wallIntRepy = sum(sum(dF .* (oy(j,k)-ywallsInt)));
        
        % Exterior wall
        dist = ((ox(j,k) - xwallsExt).^2 + (oy(j,k) - ywallsExt).^2).^0.5;
        L = floor(etaWalls./dist);
      
        % Stiffness l^2
        dF = -L.*(2*L+1).*(L+1)/3 + L.*(L+1).*etaWalls./dist;
      
        wallExtRepx = sum(sum(dF .* (ox(j,k)-xwallsExt)));
        wallExtRepy = sum(sum(dF .* (oy(j,k)-ywallsExt)));
      
        % Repulsion due to wall has to beat repulsion due to vesicles if a
        % vesicle is a sandwich between a vesicle and a wall
        repx(j) = repx(j) + wallIntRepx + wallExtRepx;
        repy(j) = repy(j) + wallIntRepy + wallExtRepy; 
          
      end % diffDiscWalls
    end % confined
  end % o.N
  % repulsion on the kth vesicle multiplied with strength W
  rep(:,k) = W*[repx;repy];
end % o.nv
end % repulsionSchemeSimple
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = newBending(o)
% f = newBending computes the additional two terms
% that come from defining the energy as 
% \int(\kappa - \tilde{\kappa})^2 where \tilde{\kappa} is some
% preferred configuration of the vesicle
global err

persistent curRef DcurRef D2curRef
% Only need to compute the first two derivatives of the
% desired convergence one

if isempty(DcurRef)
  a = 1;
  b = 2;
  % TODO: THIS NEEDS TO BE IMPLEMENTED BETTER
  % width and height of the ellipse that the vesicles
  % are going to try to approach
  %oc = curve;
  Xref = zeros(2*o.N,o.nv);
  % position of desired configuration
  DXref = zeros(2*o.N,o.nv);
  D2Xref = zeros(2*o.N,o.nv);
  % derivatives of the desired configuration
  theta = (0:o.N-1)'*2*pi/o.N;
  for k = 1:o.nv
    Xref(:,k) = [a*cos(theta);b*sin(theta)];
  end
  % something arbitrary for now
  DXref(1:o.N,:) = curve.arcDeriv(Xref(1:o.N,:),1,ones(o.N,1),o.IK);
  D2Xref(1:o.N,:) = curve.arcDeriv(DXref(1:o.N,:),1,ones(o.N,1),o.IK);
  DXref(o.N+1:2*o.N,:) = curve.arcDeriv(Xref(o.N+1:2*o.N,:),...
      1,ones(o.N,1),o.IK);
  D2Xref(o.N+1:2*o.N,:) = curve.arcDeriv(DXref(o.N+1:2*o.N,:),...
      1,ones(o.N,1),o.IK);
  % compute the derivatives of the desired configuration

  curRef = (DXref(1:o.N,:).*D2Xref(o.N+1:2*o.N,:) - ...
      DXref(o.N+1:2*o.N,:).*D2Xref(1:o.N,:))./...
      (DXref(1:o.N,:).^2 + DXref(o.N+1:2*o.N,:).^2).^1.5;
  % compute the curvature of the desired configuration

  DcurRef = zeros(o.N,o.nv);
  D2curRef = zeros(o.N,o.nv);
  for k=1:o.nv
    isa = 1./(DXref(1:o.N,:).^2 + DXref(o.N+1:2*o.N,:).^2).^0.5;
    DcurRef(:,k) = curve.arcDeriv(curRef(:,k),1,isa(:,k),o.IK(:,k));
    D2curRef(:,k) = curve.arcDeriv(DcurRef(:,k),1,isa(:,k),o.IK(:,k));
  end
  % compute the first and second derivatives of the curvature
end

f = zeros(2*o.N,o.nv);
% modificiation to traction jump
for k = 1:o.nv
  f(:,k) = [D2curRef(:,k);D2curRef(:,k)].*...
      [o.xt(o.N+1:2*o.N,k);-o.xt(1:o.N,k)];
  % first new term in the bending which acts in the normal direction
  f(:,k) = f(:,k) + [o.cur(:,k).*DcurRef(:,k);...
      o.cur(:,k).*DcurRef(:,k)].*o.xt(:,k);
  % second new term in the bending which acts in the tangential
  % direction
end

f = -o.kappa*f;
% constant in front is negative 1 times the bending coefficient

end % newbending


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ben,Ten,Div] = computeDerivs(o)
% [Ben,Ten,Div] = computeDerivs computes the matricies that takes a
% periodic function and maps it to the fourth derivative, tension, and
% surface divergence all with respect to arclength.  Everything in this
% routine is matrix free at the expense of having repmat calls

% In this implementation (compared to the old one below), we upsample these
% matrices as we form the block-diagonal preconditioner. We realized that
% since we upsample in TimeMatVec, we'd better upsample in the
% preconditioner to reduce the number of GMRES iterations. If not upsample,
% it takes more GMRES iterations

persistent D1 Restrict Prolong
% persistent variable for matrix that does Fourier differentiation (not
% with respect to arclength)
if o.antiAlias
  if isempty(D1) || o.N*o.uprate ~= size(D1,1)
    D1 = fft1.fourierDiff(o.N*o.uprate);
  end
  % Build  Fourier differentiation matrix if it is not constructed
  % yet.  It can be used throught the whole simulation as it only
  % depends on the spatial discretization N

  if isempty(Restrict) || o.uprate*o.N ~= size(Restrict,2)
    [Restrict,Prolong] = fft1.fourierRandP(o.N,o.N*o.uprate);
  end
  % build restriction and prolongation matrices
else
  if isempty(D1) || o.N ~= size(D1,1)
    D1 = fft1.fourierDiff(o.N);
  end
  % Build  Fourier differentiation matrix if it is not constructed yet.
  % It can be used throught the whole simulation as it only depends on
  % the spatial discretization N
end

Ben = zeros(2*o.N,2*o.N,o.nv);
Ten = zeros(2*o.N,o.N,o.nv);
Div = zeros(o.N,2*o.N,o.nv);

if o.antiAlias
  oc = curve;
  Nup = o.uprate*o.N;
  [x,y] = oc.getXY(o.X);
  xup = interpft(x,Nup); yup = interpft(y,Nup);
  vesicleUp = capsules([xup;yup],[],[],o.kappa,1,true);
end

if o.antiAlias
  for k = 1:o.nv
    isa = vesicleUp.isa(:,k);
    arcDeriv = isa(:,ones(Nup,1)).*D1;
    % This line is equivalent to repmat(o.isa(:,k),1,o.N).*D1 but much
    % faster.

    D4 = arcDeriv*arcDeriv; D4 = D4*D4;
    BenUp = [D4 zeros(Nup); zeros(Nup) D4];
    Ben(:,:,k) = [Restrict zeros(o.N,Nup); zeros(o.N,Nup) Restrict] * BenUp * ...
        [Prolong zeros(Nup,o.N);zeros(Nup,o.N) Prolong];

    TenUp = [arcDeriv*diag(vesicleUp.xt(1:Nup,k));...
                 arcDeriv*diag(vesicleUp.xt(Nup+1:end,k))];
    Ten(:,:,k) = [Restrict zeros(o.N,Nup); zeros(o.N,Nup) Restrict] * TenUp * ...
        Prolong;

    DivUp = [diag(vesicleUp.xt(1:Nup,k))*arcDeriv ...
                  diag(vesicleUp.xt(Nup+1:end,k))*arcDeriv];
    Div(:,:,k) = Restrict * DivUp *  ...
        [Prolong zeros(Nup,o.N);zeros(Nup,o.N) Prolong];
  end

else
  for k = 1:o.nv
    % compute single arclength derivative matrix
    isa = o.isa(:,k);
    arcDeriv = isa(:,ones(o.N,1)).*D1;
    % This line is equivalent to repmat(o.isa(:,k),1,o.N).*D1 but much
    % faster.

    D4 = arcDeriv*arcDeriv; D4 = D4*D4;
    Ben(:,:,k) = [D4 zeros(o.N); zeros(o.N) D4];

    Ten(:,:,k) = [arcDeriv*diag(o.xt(1:o.N,k));...
                 arcDeriv*diag(o.xt(o.N+1:end,k))];

    Div(:,:,k) = [diag(o.xt(1:o.N,k))*arcDeriv ...
                  diag(o.xt(o.N+1:end,k))*arcDeriv];
  end
end
Ben = real(Ben);
Ten = real(Ten);
Div = real(Div);
% Imaginary part should be 0 since we are preforming a real operation

end % computeDerivs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ben,Ten,Div] = computeDerivsOld(o)
% [Ben,Ten,Div] = computeDerivs computes the matricies that takes a
% periodic function and maps it to the fourth derivative, tension, and
% surface divergence all with respect to arclength.  Everything in this
% routine is matrix free at the expense of having repmat calls

% see above for the newest version that we use currently.

persistent D1 

% persistent variable for matrix that does Fourier differentiation (not with respect to arclength)
if isempty(D1);
  D1 = fft1.fourierDiff(o.N);
end
% Build  Fourier differentiation matrix if it is not constructed yet.
% It can be used throught the whole simulation as it only depends on
% the spatial discretization N

Ben = zeros(2*o.N,2*o.N,o.nv);
Ten = zeros(2*o.N,o.N,o.nv);
Div = zeros(o.N,2*o.N,o.nv);

for k = 1:o.nv
  % compute single arclength derivative matrix
  isa = o.isa(:,k);
  arcDeriv = isa(:,ones(o.N,1)).*D1;
  % This line is equivalent to repmat(o.isa(:,k),1,o.N).*D1 but much
  % faster.

  D4 = arcDeriv*arcDeriv; D4 = D4*D4;
  Ben(:,:,k) = [D4 zeros(o.N); zeros(o.N) D4];

  Ten(:,:,k) = [arcDeriv*diag(o.xt(1:o.N,k));...
               arcDeriv*diag(o.xt(o.N+1:end,k))];

  Div(:,:,k) = [diag(o.xt(1:o.N,k))*arcDeriv ...
                diag(o.xt(o.N+1:end,k))*arcDeriv];
end

Ben = real(Ben);
Ten = real(Ten);
Div = real(Div);
% Imaginary part should be 0 since we are preforming a real operation

end % computeDerivsOld

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma,eta,etaInt,etaExt,RS,vesVel,iter] = computeSigAndEta(...
    vesicle,tt,walls,wallsInt,wallsExt)
% computes the tension, density function (if flow is confined), and
% velocity of the vesicle for the vesicle configuration given by
% vesicle.X and solid wall configuration given by walls.X, if there are two
% sets of walls with different disc. then we have wallsInt and wallsExt, 
% instead of walls.

% flag for whether there are two sets of walls with different
% discretizations
diffDiscWalls = tt.diffDiscWalls;

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles
if tt.confined
  if ~diffDiscWalls  
    Nbd = walls.N; % points per solid wall
    nvbd = walls.nv; % number of solid walls
    NbdInt = 0;
    NbdExt = 0;
    nvbdInt = 0;
    nvbdExt = 0;
  else
    NbdInt = wallsInt.N; % points per inner wall
    NbdExt = wallsExt.N; % points on exterior wall
    nvbdInt = wallsInt.nv; % number of inner walls
    nvbdExt = wallsExt.nv; % number of outer walls (assuming 1 now)
    Nbd = 0;
    nvbd = nvbdInt+nvbdExt;
  end
else
  Nbd = 0;
  nvbd = 0;
  NbdInt = 0;
  NbdExt = 0;
  nvbdInt = 0;
  nvbdExt = 0;
end

op = tt.op;

% compute self bending, tension, and divegence terms
[Ben,Ten,Div] = vesicle.computeDerivs;

% single-layer potential
tt.Galpert = op.stokesSLmatrix(vesicle);
%tt.Gbarnett = op.laplaceSLcomplexMatrix(vesicle);
tt.Gbarnett = [];

% double-layer potential
tt.D = op.stokesDLmatrix(vesicle);
% also build the DLP without correction if fmm is on
if tt.fmmDLP
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
    interpft(vesicle.X(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0);
  tt.DLPnoCorr = op.stokesDLmatrixNoCorr(vesicleUp);
  
else
  tt.DLPnoCorr = [];
end

if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
    interpft(vesicle.X(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0);    
  tt.SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);  
else
  tt.SLPnoCorr = [];  
end

% near-singular integration structures
if tt.confined
 if ~diffDiscWalls
   [tt.NearV2V,tt.NearV2W] = vesicle.getZone(walls,3);
   [tt.NearW2W,tt.NearW2V] = walls.getZone(vesicle,3);
   
 else
   [tt.NearV2V,tt.NearV2Wint] = vesicle.getZone(wallsInt,3);
   [tt.NearWint2Wint,tt.NearWint2V] = wallsInt.getZone(vesicle,3);
   

   [~,tt.NearV2Wext] = vesicle.getZone(wallsExt,2);
   [~,tt.NearWext2V] = wallsExt.getZone(vesicle,2);

   [~,tt.NearWint2Wext] = wallsInt.getZone(wallsExt,2);
   
   [~,tt.NearWext2Wint] = wallsExt.getZone(wallsInt,2); 
 end % diffDiscWalls
else
  tt.NearV2V = vesicle.getZone([],1);
  tt.NearW2V = [];
  tt.NearV2W = [];
  tt.NearW2W = [];
  tt.NearV2Wint = [];
  tt.NearWint2Wint = [];
  tt.NearWint2V = [];
  tt.NearV2Wext = [];
  tt.NearWext2V = [];
  tt.NearWint2Wext = [];
  tt.NearWext2Wint = [];
end % confined

% initalize the right-hand side for [vesicle velocity;tension] stacked
% over each vesicle followed by the density function followed by the
% rotlets and stokeslets
if tt.confined
  if ~diffDiscWalls  
    rhs = zeros(3*N*nv+ 2*Nbd*nvbd + 3*(nvbd-1),1);
  else
    rhs = zeros(3*N*nv + 2*NbdExt*nvbdExt + 2*NbdInt*nvbdInt + 3*(nvbd-1),1);
  end
else
  rhs = zeros(3*N*nv,1);
end

% kernel for single-layer potential.  Only difference is if the FMM is
% used or not
if ~tt.fmm
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
  kernelDirect = @op.exactStokesSL;
end

% traction jump due only to the position
f = vesicle.tracJump(vesicle.X,zeros(N,nv));

% single-layer potential due to all vesicles except the diagonal one
SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);
Fslp = op.nearSingInt(vesicle,f,SLP,tt.SLPnoCorr,...
  tt.NearV2V,kernel,kernelDirect,vesicle,true,false);

% Evaluate single-layer potential due to all vesicles on
% the solid walls WITH near-singular integration
if tt.confined
  if ~diffDiscWalls
    FSLPwall = op.nearSingInt(vesicle,f,SLP,[],...
      tt.NearV2W,kernel,kernelDirect,walls,false,false);
  else
    FSLPwallInt = op.nearSingInt(vesicle,f,SLP,[],...
      tt.NearV2Wint,kernel,kernelDirect,wallsInt,false,false);
    FSLPwallExt = op.nearSingInt(vesicle,f,SLP,[],...
      tt.NearV2Wext,kernel,kernelDirect,wallsExt,false,false);    
  end
end


if ~tt.confined
  % If doing unbounded flow, add in the background velocity  
  Ffar = tt.farField(vesicle.X,[]);
else
  % no slip boundary condition for velocity on solid walls  
  Ffar = zeros(2*N,nv);
  if ~diffDiscWalls
    [U,~] = tt.farField(walls.X,[]);
  else
    [Uext,Uint] = tt.farField(wallsExt.X,wallsInt.X);  
  end  
end

velBen = Fslp + op.exactStokesSLdiag(vesicle,tt.Galpert,f) + Ffar;

for k = 1:nv
  istart = (k-1)*3*N + 1;
  iend = istart + 2*N - 1;
  rhs(istart:iend) = velBen(:,k);
end


if ~diffDiscWalls
  for k = 1:nvbd
    istart = 3*nv*N+1 + (k-1)*2*Nbd;
    iend = istart + 2*Nbd - 1;
    rhs(istart:iend) = U(:,k) - FSLPwall(:,k);
  end
else 
  % assuming one exterior wall  
  rhs(3*nv*N+1:3*nv*N+2*NbdExt) = Uext - FSLPwallExt; 
  for k = 1:nvbdInt
    istart = 3*nv*N+2*NbdExt + 1 + (k-1)*2*NbdInt;
    iend = istart + 2*NbdInt - 1;
    rhs(istart:iend) = Uint(:,k) - FSLPwallInt(:,k);
  end
end

% build vesicle part of the block-diagonal preconditioner
% GMRES
alpha = (1+vesicle.viscCont)/2;
tt.bdiagVes = [];
for k = 1:nv
  if any(vesicle.viscCont ~= 1)
    [tt.bdiagVes.L(:,:,k),tt.bdiagVes.U(:,:,k)] = lu( ...
      [(alpha(k)*eye(2*N) - tt.D(:,:,k)) ...
          -tt.Galpert(:,:,k)*Ten(:,:,k); ...
      Div(:,:,k) zeros(N)]);
  else
    [tt.bdiagVes.L(:,:,k),tt.bdiagVes.U(:,:,k)] = lu(...
      [eye(2*N) -tt.Galpert(:,:,k)*Ten(:,:,k); ...
      Div(:,:,k) zeros(N)]);
  end
end

% solve for tension and density function with block-diagonal 
% preconditioning
maxGMRESit = min(tt.gmresMaxIter,N*nv+2*Nbd*nvbd + 2*NbdExt*nvbdExt + ...
    2*NbdInt*nvbdInt + 3*(nvbd-1));
[sigDen,F,R,I] = gmres(@(X) tt.sigDenMatVec(X,vesicle,walls,...
    wallsInt,wallsExt),rhs,[],tt.gmresTol,maxGMRESit,@tt.preconditionerBD);

iter = I(2);

vesVel = zeros(2*N,nv);
sigma = zeros(N,nv); % tension
eta = zeros(2*Nbd,nvbd); % solid wall density
etaInt = zeros(2*NbdInt,nvbdInt); % inner wall density
etaExt = zeros(2*NbdExt,nvbdExt); % outer wall density
RS = zeros(3,nvbd); % rotlets and stokeslets

% unstack the tension and density function from the GMRES solver
for k = 1:nv
  vesVel(:,k) = sigDen((k-1)*3*N+1:(3*k-1)*N);
  sigma(:,k) = sigDen((3*k-1)*N+1:3*k*N);
end

if ~diffDiscWalls
  for k = 1:nvbd
    eta(:,k) = sigDen(3*nv*N+(k-1)*2*Nbd+1:3*nv*N+k*2*Nbd);
  end
else
  % assuming one outer wall  
  etaExt = sigDen(3*nv*N+1:3*nv*N+2*NbdExt);
  for k = 1:nvbdInt
    istart = 3*nv*N+2*NbdExt+(k-1)*2*NbdInt+1;
    iend = istart+2*NbdInt-1;
    etaInt(:,k) = sigDen(istart:iend);
  end  
end

if ~diffDiscWalls
  sizeSoFar = 3*nv*N+2*nvbd*Nbd;
else
  sizeSoFar = 3*nv*N+2*NbdExt*nvbdExt+2*NbdInt*nvbdInt;  
end

for k = 2:nvbd  
  istart = sizeSoFar + 3*(k-2) + 1;
  iend = istart + 2;
  RS(:,k) = sigDen(istart:iend);
end


end % computeSigAndEta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma,vesVel,iter,tt,sigDen] = ...
        computeSigmaWRegul(vesicle,tt,sigDenOld)
% Computes tension using regularized kernel for singular integrals,
% near-singular integrals are not corrected and computed by FMM

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles
op = tt.op;

% compute self bending, tension, and divegence terms
[Ben,Ten,Div] = vesicle.computeDerivs;

% single-layer matrix w/ zero assigned to the diagonals
tt.SLPnoCorr = op.stokesSLmatrixNoCorr(vesicle);
SLPnoCorrRegul = op.stokesSLmatrixNoCorrRegul(vesicle);
% initalize the right-hand side for [vesicle velocity;tension] stacked
% over each vesicle followed by the density function followed by the
% rotlets and stokeslets
rhs = zeros(3*N*nv,1);

% traction jump due only to the position
f = vesicle.tracJump(vesicle.X,zeros(N,nv));

% single-layer potential due to all vesicles except the diagonal one
% here the self-interaction is subtracted
if tt.fmm
  Fslp = op.exactStokesSLfmm(vesicle,f,tt.SLPnoCorr);
else
  Fslp = op.exactStokesSL(vesicle,f,tt.SLPnoCorr);  
end


% If doing unbounded flow, add in the background velocity  
Ffar = tt.farField(vesicle.X,[]);

% Here, the self-interaction is computed with regularized kernel and added
velBen = Fslp + op.exactStokesSLdiag(vesicle,tt.SLPnoCorr,f) + Ffar;

for k = 1:nv
  istart = (k-1)*3*N + 1;
  iend = istart + 2*N - 1;
  rhs(istart:iend) = velBen(:,k);
end

% build vesicle part of the block-diagonal preconditioner
% GMRES
alpha = (1+vesicle.viscCont)/2;
for k = 1:nv
  [tt.bdiagVes.L(:,:,k),tt.bdiagVes.U(:,:,k)] = lu(...
    [eye(2*N) -0.5*SLPnoCorrRegul(:,:,k)*Ten(:,:,k); ...
    Div(:,:,k) zeros(N)]);
end

% solve for tension and density function with block-diagonal 
% preconditioning
% maxGMRESit = min(tt.gmresMaxIter,3*N*nv);
maxGMRESit = 2;
% gmresTol = tt.gmresTol;
gmresTol = 1e-4;

[sigDen,F,R,I] = gmres(@(X) tt.sigDenMatVecLowAcc(X,vesicle),...
    rhs,[],gmresTol,maxGMRESit,@tt.preconditionerBD,[],sigDenOld);
iter = I(2);

vesVel = zeros(2*N,nv);
sigma = zeros(N,nv); % tension

% unstack the tension and density function from the GMRES solver
for k = 1:nv
  vesVel(:,k) = sigDen((k-1)*3*N+1:(3*k-1)*N);  
  sigma(:,k) = sigDen((3*k-1)*N+1:3*k*N);
end


end % computeSigmaWRegul

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma,eta,etaInt,etaExt,RS,vesVel,iter,tt,sigDen] = ...
        computeSigAndEtaTToutInEx(vesicle,tt,walls,wallsInt,wallsExt)
% computes the tension, density function (if flow is confined), and
% velocity of the vesicle for the vesicle configuration given by
% vesicle.X and solid wall configuration given by walls.X, if there are two
% sets of walls with different disc. then we have wallsInt and wallsExt, 
% instead of walls.

% flag for whether there are two sets of walls with different
% discretizations
diffDiscWalls = tt.diffDiscWalls;

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles
if tt.confined
  if ~diffDiscWalls  
    Nbd = walls.N; % points per solid wall
    nvbd = walls.nv; % number of solid walls
    NbdInt = 0;
    NbdExt = 0;
    nvbdInt = 0;
    nvbdExt = 0;
  else
    NbdInt = wallsInt.N; % points per inner wall
    NbdExt = wallsExt.N; % points on exterior wall
    nvbdInt = wallsInt.nv; % number of inner walls
    nvbdExt = wallsExt.nv; % number of outer walls (assuming 1 now)
    Nbd = 0;
    nvbd = nvbdInt+nvbdExt;
  end
else
  Nbd = 0;
  nvbd = 0;
  NbdInt = 0;
  NbdExt = 0;
  nvbdInt = 0;
  nvbdExt = 0;
end

op = tt.op;

% compute self bending, tension, and divegence terms
[Ben,Ten,Div] = vesicle.computeDerivs;

% single-layer potential
tt.Galpert = op.stokesSLmatrix(vesicle);
%tt.Gbarnett = op.laplaceSLcomplexMatrix(vesicle);
tt.Gbarnett = [];

% double-layer potential
tt.D = op.stokesDLmatrix(vesicle);
% also build the DLP without correction if fmm is on
if tt.fmmDLP
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
    interpft(vesicle.X(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0);
  tt.DLPnoCorr = op.stokesDLmatrixNoCorr(vesicleUp);
  
else
  tt.DLPnoCorr = [];
end

if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(vesicle.X(1:end/2,:),Nup);...
    interpft(vesicle.X(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0);    
  tt.SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);  
else
  tt.SLPnoCorr = [];  
end

% near-singular integration structures
if tt.confined
 if ~diffDiscWalls
   [tt.NearV2V,tt.NearV2W] = vesicle.getZone(walls,3);
   [tt.NearW2W,tt.NearW2V] = walls.getZone(vesicle,3);
   
 else
   [tt.NearV2V,tt.NearV2Wint] = vesicle.getZone(wallsInt,3);
   [tt.NearWint2Wint,tt.NearWint2V] = wallsInt.getZone(vesicle,3);
   

   [~,tt.NearV2Wext] = vesicle.getZone(wallsExt,2);
   [~,tt.NearWext2V] = wallsExt.getZone(vesicle,2);

   [~,tt.NearWint2Wext] = wallsInt.getZone(wallsExt,2);
   
   [~,tt.NearWext2Wint] = wallsExt.getZone(wallsInt,2); 
 end % diffDiscWalls
else
  tt.NearV2V = vesicle.getZone([],1);
  tt.NearW2V = [];
  tt.NearV2W = [];
  tt.NearW2W = [];
  tt.NearV2Wint = [];
  tt.NearWint2Wint = [];
  tt.NearWint2V = [];
  tt.NearV2Wext = [];
  tt.NearWext2V = [];
  tt.NearWint2Wext = [];
  tt.NearWext2Wint = [];
end % confined

% initalize the right-hand side for [vesicle velocity;tension] stacked
% over each vesicle followed by the density function followed by the
% rotlets and stokeslets
if tt.confined
  if ~diffDiscWalls  
    rhs = zeros(3*N*nv+ 2*Nbd*nvbd + 3*(nvbd-1),1);
  else
    rhs = zeros(3*N*nv + 2*NbdExt*nvbdExt + 2*NbdInt*nvbdInt + 3*(nvbd-1),1);
  end
else
  rhs = zeros(3*N*nv,1);
end

% kernel for single-layer potential.  Only difference is if the FMM is
% used or not
if ~tt.fmm
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
  kernelDirect = @op.exactStokesSL;
end

% traction jump due only to the position
f = vesicle.tracJump(vesicle.X,zeros(N,nv));

% single-layer potential due to all vesicles except the diagonal one
SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);
Fslp = op.nearSingInt(vesicle,f,SLP,tt.SLPnoCorr,...
  tt.NearV2V,kernel,kernelDirect,vesicle,true,false);

% Evaluate single-layer potential due to all vesicles on
% the solid walls WITH near-singular integration
if tt.confined
  if ~diffDiscWalls
    FSLPwall = op.nearSingInt(vesicle,f,SLP,[],...
      tt.NearV2W,kernel,kernelDirect,walls,false,false);
  else
    FSLPwallInt = op.nearSingInt(vesicle,f,SLP,[],...
      tt.NearV2Wint,kernel,kernelDirect,wallsInt,false,false);
    FSLPwallExt = op.nearSingInt(vesicle,f,SLP,[],...
      tt.NearV2Wext,kernel,kernelDirect,wallsExt,false,false);    
  end
end


if ~tt.confined
  % If doing unbounded flow, add in the background velocity  
  Ffar = tt.farField(vesicle.X,[]);
else
  % no slip boundary condition for velocity on solid walls  
  Ffar = zeros(2*N,nv);
  if ~diffDiscWalls
    [U,~] = tt.farField(walls.X,[]);
  else
    [Uext,Uint] = tt.farField(wallsExt.X,wallsInt.X);  
  end  
end

velBen = Fslp + op.exactStokesSLdiag(vesicle,tt.Galpert,f) + Ffar;

for k = 1:nv
  istart = (k-1)*3*N + 1;
  iend = istart + 2*N - 1;
  rhs(istart:iend) = velBen(:,k);
end


if ~diffDiscWalls
  for k = 1:nvbd
    istart = 3*nv*N+1 + (k-1)*2*Nbd;
    iend = istart + 2*Nbd - 1;
    rhs(istart:iend) = U(:,k) - FSLPwall(:,k);
  end
else 
  % assuming one exterior wall  
  rhs(3*nv*N+1:3*nv*N+2*NbdExt) = Uext - FSLPwallExt; 
  for k = 1:nvbdInt
    istart = 3*nv*N+2*NbdExt + 1 + (k-1)*2*NbdInt;
    iend = istart + 2*NbdInt - 1;
    rhs(istart:iend) = Uint(:,k) - FSLPwallInt(:,k);
  end
end

% build vesicle part of the block-diagonal preconditioner
% GMRES
alpha = (1+vesicle.viscCont)/2;
for k = 1:nv
  if any(vesicle.viscCont ~= 1)
    [tt.bdiagVes.L(:,:,k),tt.bdiagVes.U(:,:,k)] = lu( ...
      [(alpha(k)*eye(2*N) - tt.D(:,:,k)) ...
          -tt.Galpert(:,:,k)*Ten(:,:,k); ...
      Div(:,:,k) zeros(N)]);
  else
    [tt.bdiagVes.L(:,:,k),tt.bdiagVes.U(:,:,k)] = lu(...
      [eye(2*N) -tt.Galpert(:,:,k)*Ten(:,:,k); ...
      Div(:,:,k) zeros(N)]);
  end
end

% solve for tension and density function with block-diagonal 
% preconditioning
maxGMRESit = min(tt.gmresMaxIter,N*nv+2*Nbd*nvbd + 2*NbdExt*nvbdExt + ...
    2*NbdInt*nvbdInt + 3*(nvbd-1));

[sigDen,F,R,I] = gmres(@(X) tt.sigDenMatVec(X,vesicle,walls,...
    wallsInt,wallsExt),rhs,[],tt.gmresTol,maxGMRESit,@tt.preconditionerBD);

iter = I(2);

vesVel = zeros(2*N,nv);
sigma = zeros(N,nv); % tension
eta = zeros(2*Nbd,nvbd); % solid wall density
etaInt = zeros(2*NbdInt,nvbdInt); % inner wall density
etaExt = zeros(2*NbdExt,nvbdExt); % outer wall density
RS = zeros(3,nvbd); % rotlets and stokeslets

% unstack the tension and density function from the GMRES solver
for k = 1:nv
  vesVel(:,k) = sigDen((k-1)*3*N+1:(3*k-1)*N);
  sigma(:,k) = sigDen((3*k-1)*N+1:3*k*N);
end

if ~diffDiscWalls
  for k = 1:nvbd
    eta(:,k) = sigDen(3*nv*N+(k-1)*2*Nbd+1:3*nv*N+k*2*Nbd);
  end
else
  % assuming one outer wall  
  etaExt = sigDen(3*nv*N+1:3*nv*N+2*NbdExt);
  for k = 1:nvbdInt
    istart = 3*nv*N+2*NbdExt+(k-1)*2*NbdInt+1;
    iend = istart+2*NbdInt-1;
    etaInt(:,k) = sigDen(istart:iend);
  end  
end

if ~diffDiscWalls
  sizeSoFar = 3*nv*N+2*nvbd*Nbd;
else
  sizeSoFar = 3*nv*N+2*NbdExt*nvbdExt+2*NbdInt*nvbdInt;  
end

for k = 2:nvbd  
  istart = sizeSoFar + 3*(k-2) + 1;
  iend = istart + 2;
  RS(:,k) = sigDen(istart:iend);
end


end % computeSigAndEtaTToutInEx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eta,etaInt,etaExt,RS,iter] = computeEta(walls,wallsExt,tt)
% [eta,etaInt,etaExt,RS,iter] = computeEta(walls,walls2,tt) 
% computes wall densities and RS when
% there is no vesicle but we want to compute the tracers. This is called by
% tstep.m/noVesTracers. It uses the the exact inverse of
% the double-layer potential for stokes flow in a bounded domain which is
% computed by tstep.m/wallsPrecond. if wallsExt is not empty, that means that
% there are two sets of walls with different discretizations. So this
% outputs etaInt and etaExt, instead of eta. In that case, walls is
% wallsInt.

% if walls2 is not empty, then there are two sets of walls with diff disc.
diffDiscWalls = 0;
if ~isempty(wallsExt)
  diffDiscWalls = 1;
end

if ~diffDiscWalls
  Nbd = walls.N; % points per solid wall
  nvbd = walls.nv; % number of solid walls
  NbdInt = 0;
  nvbdInt = 0;
  NbdExt = 0;
  nvbdExt = 0;
else
  NbdInt = walls.N; % points per inner wall
  nvbdInt = walls.nv; % number of inner walls
  NbdExt = wallsExt.N; % points on exterior wall
  nvbdExt = wallsExt.nv; % number of exterior walls (for now, it is 1)
  nvbd = nvbdInt+nvbdExt;
  Nbd = 0;
end

rhs = zeros(2*Nbd*nvbd+2*NbdInt*nvbdInt+2*NbdExt*nvbdExt+3*(nvbd-1),1);

% no slip boundary condition for velocity on solid walls
if ~diffDiscWalls
  U = tt.farField(walls.X,[]);
else
  [Uext, Uint] = tt.farField(wallsExt.X,walls.X);
end

% Form the right hand side
if ~diffDiscWalls
  for k = 1:nvbd
    istart = 1 + (k-1)*2*Nbd;
    iend = k*2*Nbd;
    rhs(istart:iend) = U(:,k);
  end
else
  rhs(1:2*NbdExt) = Uext;
  for k = 1 : nvbdInt
    istart = 2*NbdExt + 1 + (k-1)*2*NbdInt;
    iend = istart + 2*NbdInt - 1;
    rhs(istart:iend) = Uint(:,k);
  end  
end

% Direct solve for the system as tt.bdiagWall is an exact factorization
etaRS = tt.bdiagWall*rhs;
iter = 0;

% unstack the density and RS function from the GMRES solver
if ~diffDiscWalls
  eta = zeros(2*Nbd,nvbd); % solid wall density
else
  etaInt = zeros(2*NbdInt,nvbdInt);
  etaExt = zeros(2*NbdExt,nvbdExt);
end

RS = zeros(3,nvbd); % rotlets and stokeslets

if ~diffDiscWalls
  for k = 1:nvbd
    eta(:,k) = etaRS((k-1)*2*Nbd+1:k*2*Nbd);
  end
  sizeSoFar = 2*nvbd*Nbd;
  etaExt = [];
  etaInt = [];
else
  etaExt = etaRS(1:2*NbdExt); % assuming one outer wall
  for k = 1:nvbdInt
    istart = 2*NbdExt+(k-1)*2*NbdInt+1;
    iend = istart+2*NbdInt-1;
    etaInt(:,k) = etaRS(istart:iend);
  end
  sizeSoFar = 2*nvbdInt*NbdInt+2*nvbdExt*NbdExt;
  eta = [];
end

for k = 2:nvbd
  istart = sizeSoFar + 3*(k-2) + 1;
  iend = istart + 2;
  RS(:,k) = etaRS(istart:iend);
end


end % computeEta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [derivA,derivL] = computeAreaLengthDerivs(o,X,Xold,dt)
% [derivA,derivL] = computeAreaLengthDerivs(o,X,Xold,dt) computes dA/dt 
% and dL/dt based on the area, length computation through parameterization
oc = curve;

N = size(X,1)/2;
nv = size(X,2);
Nup = 96;

[x,y] = oc.getXY(X);
[xOld,yOld] = oc.getXY(Xold);

% Upsample
xUp = interpft(x,Nup); yUp = interpft(y,Nup);
xOldUp = interpft(xOld,Nup); yOldUp = interpft(yOld,Nup);


% Compute velocity on upsampled grid using finite differences
vel = ([xUp;yUp]-[xOldUp;yOldUp])/dt;

derivA = zeros(nv,1);
derivL = zeros(nv,1);

for k = 1 : nv  
    % Compute dx/da of vesicle k
    [DxOld,DyOld] = oc.getDXY([xOldUp(:,k);yOldUp(:,k)]);
    [Dx,Dy] = oc.getDXY([xUp(:,k);yUp(:,k)]);
    % Compute dv/da of vesicle k
    [DvelX,DvelY] = oc.getDXY(vel(:,k));
    % Get x and y components of velocity of vesicle k
    [velx,vely] = oc.getXY(vel(:,k));
    
    % For dA/dt
    % Compute the integrand and integrate it with equidistributed alpha
    integA = (velx.*DyOld + xOldUp(:,k).*DvelY - vely.*DxOld-yOldUp(:,k).*DvelX);
    
    derivA(k) = sum(integA)*pi/Nup;
    
    % For dL/dt
    % Compute the integrand and integrate it with equidistributed alpha
    integL = (DvelX.*DxOld + DvelY.*DyOld)./sqrt(DxOld.^2 + DyOld.^2);
    derivL(k) = sum(integL)*2*pi/Nup;
end
    
end %computeAreaLengthDerivs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NearSelf,NearOther] = getZone(vesicle1,vesicle2,relate)
% [NearSelf,NearOther] = getZone(vesicle1,vesicle2,relate) constructs
% each vesicle, index of the closest point, nearest point on a local
% interapolant, and argument of that nearest point.  vesicle1 contains
% the source points (which are also target points) and vesicle2 
% contains additional target points.  The 
% values of relate corresond to
% relate == 1  => only require NearSelf  (ie. vesicle to vesicle)
% relate == 2  => only require NearOther (ie. vesicle to wall)
% relate == 3  => require both NearSelf and NearOther
% THIS ROUTINE HAS A LOOP OVER THE TARGET POINTS WHICH SHOULD BE REMOVED
NearSelf = [];
NearOther = [];

N1 = vesicle1.N; % number of source/target points
nv1 = vesicle1.nv; % number of source/target vesicles
X1 = vesicle1.X; % source and target points
oc = curve;
[xsou,ysou] = oc.getXY(X1); 
% separate targets into x and y coordinates

h = vesicle1.length/N1; 
% smallest arclength over all vesicles
ptsperbox = 10; 
% Estimate for number of points per box.  This simply sets the 
% number of uniformly refined boxes we take.  Estimate is not very
% accurate.  What ptsperbox represents is the total number of points
% that could be put in each two-dimensional bin where no two are
% less than distance h from one another.  However, our points live
% on curves and thus will not fill up an entire bin

H = sqrt(ptsperbox)*h;
xmin = min(min(xsou));
xmax = max(max(xsou));
xmin = xmin - H;
xmax = xmax + H;
ymin = min(min(ysou));
ymax = max(max(ysou));
ymin = ymin - H;
ymax = ymax + H;
% Add a buffer around the points so that it is easier to
% work with vesicle2

Nx = ceil((xmax - xmin)/H);
Ny = ceil((ymax - ymin)/H);
% Find bounds for box that contains all points and add a buffer 
% so that all points are guaranteed to be in the box

Nbins = Nx * Ny; % Total number of bins

ii = ceil((xsou - xmin)/H);
jj = ceil((ysou - ymin)/H);
% Index in x and y direction of the box containing each point
bin = (jj-1)*Nx + ii;
% Find bin of each point using lexiographic ordering (x then y)


%figure(2);
%clf; hold on
%plot(xsou,ysou,'k.')
%axis equal
%axis([xmin xmin+Nx*H ymin ymin+Ny*H])
%set(gca,'xtick',linspace(xmin,xmin+Nx*H,Nx+1))
%set(gca,'ytick',linspace(ymin,ymin+Ny*H,Ny+1))
%grid on
%set(gca,'xticklabel',[])
%set(gca,'yticklabel',[])
%figure(1)
%pause
% DEBUG: This does a simple plot of the points with a grid that 
% aligns with the boundary of the boxes

%whos
%disp([Nbins nv1])
%disp([xmin xmax ymin ymax])
%disp([h H])
%pause
fpt = zeros(Nbins,nv1);
lpt = zeros(Nbins,nv1);
% allocate space for storing first and last points
[binsort,permute] = sort(bin);
% build permute.  Need binsort to find first and last points
% in each box

for k = 1:nv1 % Loop over vesicles
  for j = 1:N1 % Loop over bins
    ibox = binsort(j,k);
    if (fpt(ibox,k) == 0)
      fpt(ibox,k) = j;
      lpt(ibox,k) = 1;
    else
      lpt(ibox,k) = lpt(ibox,k) + 1;
    end
  end
  lpt(:,k) = fpt(:,k) + lpt(:,k) - 1;
end
% Construct first and last point in each box corresponding
% to each vesicle.  The order is based on permute.  For example,
% permute(fpt(ibox,k)),...,permute(lpt(ibox,k)) is the set of 
% all points from vesicle k contained in box ibox

neigh = zeros(Nbins,9);

%Do corners first
neigh(1,1:4) = [1 2 Nx+1 Nx+2]; 
% bottom left corner
neigh(Nx,1:4) = [Nx Nx-1 2*Nx 2*Nx-1]; 
% bottom right corner
neigh(Nbins-Nx+1,1:4) = [Nbins-Nx+1 Nbins-Nx+2 ...
    Nbins-2*Nx+1 Nbins-2*Nx+2];
% top left corner
neigh(Nbins,1:4) = [Nbins Nbins-1 Nbins-Nx Nbins-Nx-1]; 
% top right corner

for j = 2:Nx-1
  neigh(j,1:6) = j + [-1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of bottom row

for j = Nbins-Nx+2:Nbins-1
  neigh(j,1:6) = j + [-1 0 1 -Nx-1 -Nx -Nx+1];
end
% neighbors of top row

for j=Nx+1:Nx:Nbins-2*Nx+1
  neigh(j,1:6) = j + [-Nx -Nx+1 0 1 Nx Nx+1];
end
% neighbors of left column

for j=2*Nx:Nx:Nbins-Nx
  neigh(j,1:6) = j + [-Nx-1 -Nx -1 0 Nx-1 Nx];
end
% neighbors of right column

J = (Nx + 1:Nbins - Nx);
J = J(mod(J-1,Nx)~=0);
J = J(mod(J,Nx)~=0);
% J is the index of boxes that are not on the boundary
for j=J
  neigh(j,:) = j + [-Nx-1 -Nx -Nx+1 -1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of interior points
% TREE STRUCTURE IS COMPLETE


if (relate == 1 || relate == 3)
%  distSS = -ones(N1,nv1,nv1); 
%  % dist(n,k,j) is the distance of point n on vesicle k to
%  % vesicle j
%  zoneSS = -ones(N1,nv1,nv1); 
%  % near or far zone
%  nearestSS = -ones(2*N1,nv1,nv1); 
%  % nearest point using local interpolant
%  icpSS = -ones(N1,nv1,nv1); 
%  % index of closest discretization point
%  argnearSS = -ones(N1,nv1,nv1); 
%  % argument in [0,1] of local interpolant
  for k = 1:nv1
    distSS{k} = spalloc(N1,nv1,0);
    % dist(n,k,j) is the distance of point n on vesicle k to
    zoneSS{k} = spalloc(N1,nv1,0);
    % near or far zone
    nearestSS{k} = spalloc(2*N1,nv1,0);
    % nearest point using local interpolant
    icpSS{k} = spalloc(N1,nv1,0);
    % index of closest discretization point
    argnearSS{k} = spalloc(N1,nv1,0);
    % argument in [0,1] of local interpolant
  end
  % New way of representing near-singular integration structure so that
  % we can use sparse matricies.


  % begin classifying points where we are considering 
  % vesicle to vesicle relationships
  for k = 1:nv1
    boxes = unique(bin(:,k));
    % Find all boxes containing points of vesicle k
    boxes = neigh(boxes,:);
    % Look at all neighbors of boxes containing vesicle k
    boxes = unique(boxes(:));
    % Remove repetition
    boxes = boxes(boxes~=0);
    % Delete non-existent boxes that came up because of neigh

    K = [(1:k-1) (k+1:nv1)];
    for k2 = K
      istart = fpt(boxes,k2);
      iend = lpt(boxes,k2);
      istart = istart(istart ~= 0);
      iend = iend(iend ~= -1);
      % Find index of all points in neighboring boxes of vesicle k
      % that are in vesicle k2
      
      neighpts = zeros(sum(iend-istart+1),1);
      
      % Allocate space to assign possible near points
      is = 1;
      for j=1:numel(istart)
        ie = is + iend(j) - istart(j);
        neighpts(is:ie) = permute(istart(j):iend(j),k2);
        is = ie + 1;
      end
      % neighpts contains all points on vesicle k2 that are in 
      % neighboring boxes to vesicle k

      neighpts = sort(neighpts);
      % sorting should help speedup as we won't be jumping around
      % through different boxes

      n = 0;
      for i=1:numel(neighpts)
        ipt = neighpts(i);
        ibox = bin(ipt,k2);
        % box containing ipt on vesicle k2
        if (ibox ~= n)
          n = ibox;
          % Check if we need to move to a new box
          neighbors = neigh(ibox,:);
          % neighbors of this box
          neighbors = neighbors(neighbors~=0);
          % Remove non-existent neighbors
          istart = fpt(neighbors,k);
          iend = lpt(neighbors,k);
          istart = istart(istart ~= 0);
          iend = iend(iend ~= -1);
          % Find points on vesicle k in neighboring boxes
          neighpts2 = zeros(sum(iend-istart+1),1);
          is = 1;
          for j=1:numel(istart)
            ie = is + iend(j) - istart(j);
            neighpts2(is:ie) = permute(istart(j):iend(j),k);
            is = ie + 1;
          end
          % neighpts2 contains all points on vesicle k that 
          % are in neighboring box of ibox 
        end % decide if we need to switch boxes

        [d0,d0loc] = min((xsou(ipt,k2) - xsou(:,k)).^2 + ...
            (ysou(ipt,k2) - ysou(:,k)).^2);
        % Find minimum distance between ipt on vesicle k2 to
        % possible closest points on vesicle k
        d0 = sqrt(d0);
        % Save on not taking the square root on a vector but instead
        % on a single real number

        icpSS{k}(ipt,k2) = d0loc;
        if (d0 < 2*h);
          [distSS{k}(ipt,k2),nearestx,nearesty,argnearSS{k}(ipt,k2)] = ...
              vesicle1.closestPnt([xsou;ysou],xsou(ipt,k2),...
              ysou(ipt,k2),k,icpSS{k}(ipt,k2));
          nearestSS{k}(ipt,k2) = nearestx;
          nearestSS{k}(ipt+N1,k2) = nearesty;
          % Find closest point along a local interpolant using
          % Newton's method.

          if (distSS{k}(ipt,k2) < h)
            zoneSS{k}(ipt,k2) = 1;
          end
          % Point ipt of vesicle k2 is in the near zone of
          % vesicle k
        end


      end % ipt

    end % k2

  end % k

  NearSelf.dist = distSS;
  NearSelf.zone = zoneSS;
  NearSelf.nearest = nearestSS;
  NearSelf.icp = icpSS;
  NearSelf.argnear = argnearSS;
  % Store everything in the structure NearSelf.  This way it is 
  % much cleaner to pass everything around
end % relate == 1 || relate == 3

% Bin target points with respect to the source points
if (relate == 2 || relate == 3)
  N2 = vesicle2.N; % number of additional targets
  nv2 = vesicle2.nv; % number of additional vesicles
  X2 = vesicle2.X; % additional target points
  [xtar,ytar] = oc.getXY(X2);
  % separate additional target points into x and y coordinates
%  figure(2); clf
%  plot(xtar,ytar,'r.')
%  hold on
%  plot(xsou,ysou,'k.')
%  pause
% DEBUG: FOR SEEING TARGET AND SOURCE POINTS IN THE TREE STRUCTURE
% WHICH CAN BE PLOTTED ABOVE

  for k = 1:nv1
    distST{k} = spalloc(N1,nv2,0);
    % dist(n,k,j) is the distance of point n on vesicle k to
    zoneST{k} = spalloc(N1,nv2,0);
    % near or far zone
    nearestST{k} = spalloc(2*N1,nv2,0);
    % nearest point using local interpolant
    icpST{k} = spalloc(N1,nv2,0);
    % index of closest discretization point
    argnearST{k} = spalloc(N1,nv2,0);
    % argument in [0,1] of local interpolant
  end
  % Represent near-singular integration structure using sparse matricies

  itar = ceil((xtar - xmin)/H);
  jtar = ceil((ytar - ymin)/H);
  [indx,indy] = find((itar >= 1) & (itar <= Nx) & ...
      (jtar >= 1) & (jtar <= Ny));
  % Only have to consider xx(ind),yy(ind) since all other points
  % are not contained in the box [xmin xmax] x [ymin ymax]

  for k = 1:nv1 % loop over sources
    for nind = 1:numel(indx) 
      % loop over points that are not outside the box that surrounds
      % all target points with a sufficiently large buffer
      ii = indx(nind);
      jj = indy(nind);
      binTar = (jtar(ii,jj)-1)*Nx + itar(ii,jj);
      boxesTar = neigh(binTar,:);
      boxesTar = boxesTar(boxesTar~=0);
      istart = fpt(boxesTar,k);
      iend  = lpt(boxesTar,k);
      istart = istart(istart ~= 0);
      iend   = iend(iend ~= -1);
      
      neighpts = zeros(sum(iend-istart+1),1);
  
      % Allocate space to assign possible near points
      if numel(neighpts) > 0
        % it is possible of the neighboring boxes to contain
        % no points.
        is = 1;
        
        for j = 1:numel(istart)
          ie = is + iend(j) - istart(j);
          neighpts(is:ie) = permute(istart(j):iend(j),k);
          is = ie + 1;
        end
        % Set of potentially nearest points to 
        % (xtar(jj),ytar(jj))
        
        [d0,d0loc] = min((xtar(ii,jj) - xsou(neighpts,k)).^2 + ...
          (ytar(ii,jj) - ysou(neighpts,k)).^2);
        % find closest point and distance between (xtar(jj),ytar(jj))
        % and vesicle k.  Only need to look at points in neighboring
        % boxes
        
        d0 = d0.^0.5;
        icpST{k}(ii,jj) = neighpts(d0loc);

        if d0 < 2*h
          [distST{k}(ii,jj),nearestx,nearesty,argnearST{k}(ii,jj)] = ...
            vesicle1.closestPnt([xsou;ysou],xtar(ii,jj),...
            ytar(ii,jj),k,icpST{k}(ii,jj));
          nearestST{k}(ii,jj) = nearestx;
          nearestST{k}(ii+N2,jj) = nearesty;
          
          % DEBUG: CHECK THAT NEWTON'S METHOD HAS DONE A GOOD JOB
          % CONVERGING TO THE NEAREST POINT
          % compute distance and nearest point between 
          % (xtar(ii,jj),ytar(ii,jj)) and vesicle k
          if distST{k}(ii,jj) < h
            zoneST{k}(ii,jj) = 1;
            % (xtar(ii,jj),ytar(ii,jj)) is in the near zone of vesicle k
          end
        end % d0 < 2*h
      end % numel(neighpts) > 0

    end % ii and jj

  end % k

  NearOther.dist = distST;
  NearOther.zone = zoneST;
  NearOther.nearest = nearestST;
  NearOther.icp = icpST;
  NearOther.argnear = argnearST;
  % store near-singluar integration requirements in structure NearOther

end % relate == 2 || relate == 3
end % getZone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NearSelf,NearOther] = getZoneNew(vesicle1,vesicle2,relate)
% THIS IS GOKBERK'S FIRST IMPROVEMENT TO getZone.  I DON'T KNOW IF IT
% IS FASTER IN ALL CASES, SO THERE IS STILL WORK TO BE DONE
% [NearSelf,NearOther] = getZone(vesicle1,vesicle2,relate) constructs
% each vesicle, index of the closest point, nearest point on a local
% interapolant, and argument of that nearest point.  vesicle1 contains
% the source points (which are also target points) and vesicle2 
% contains additional target points.  The 
% values of relate corresond to
% relate == 1  => only require NearSelf  (ie. vesicle to vesicle)
% relate == 2  => only require NearOther (ie. vesicle to wall)
% relate == 3  => require both NearSelf and NearOther
NearSelf = [];
NearOther = [];

N1 = vesicle1.N; % number of source/target points
nv1 = vesicle1.nv; % number of source/target vesicles
X1 = vesicle1.X; % source and target points
oc = curve;
[xsou,ysou] = oc.getXY(X1); 
% separate targets into x and y coordinates

h = vesicle1.length/N1; 
% smallest arclength over all vesicles
ptsperbox = 10; 
% Estimate for number of points per box.  This simply sets the 
% number of uniformly refined boxes we take.  Estimate is not very
% accurate.  What ptsperbox represents is the total number of points
% that could be put in each two-dimensional bin where no two are
% less than distance h from one another.  However, our points live
% on curves and thus will not fill up an entire bin

H = sqrt(ptsperbox)*h;
xmin = min(min(xsou));
xmax = max(max(xsou));
xmin = xmin - H;
xmax = xmax + H;
ymin = min(min(ysou));
ymax = max(max(ysou));
ymin = ymin - H;
ymax = ymax + H;
% Add a buffer around the points so that it is easier to
% work with vesicle2

Nx = ceil((xmax - xmin)/H);
Ny = ceil((ymax - ymin)/H);
% Find bounds for box that contains all points and add a buffer 
% so that all points are guaranteed to be in the box

Nbins = Nx * Ny; % Total number of bins

ii = ceil((xsou - xmin)/H);
jj = ceil((ysou - ymin)/H);
% Index in x and y direction of the box containing each point
bin = (jj-1)*Nx + ii;
% Find bin of each point using lexiographic ordering (x then y)


%figure(2);
%clf; hold on
%plot(xsou,ysou,'k.')
%axis equal
%axis([xmin xmin+Nx*H ymin ymin+Ny*H])
%set(gca,'xtick',linspace(xmin,xmin+Nx*H,Nx+1))
%set(gca,'ytick',linspace(ymin,ymin+Ny*H,Ny+1))
%grid on
%figure(1)
%pause
% DEBUG: This does a simple plot of the points with a grid that 
% aligns with the boundary of the boxes

%whos
%disp([Nbins nv1])
%disp([xmin xmax ymin ymax])
%disp([h H])
fpt = zeros(Nbins,nv1);
lpt = zeros(Nbins,nv1);
% allocate space for storing first and last points
[binsort,permute] = sort(bin);
% build permute.  Need binsort to find first and last points
% in each box

for k = 1:nv1 % Loop over vesicles
  for j = 1:N1 % Loop over bins
    ibox = binsort(j,k);
    if (fpt(ibox,k) == 0)
      fpt(ibox,k) = j;
      lpt(ibox,k) = 1;
    else
      lpt(ibox,k) = lpt(ibox,k) + 1;
    end
  end
  lpt(:,k) = fpt(:,k) + lpt(:,k) - 1;
end
% Construct first and last point in each box corresponding
% to each vesicle.  The order is based on permute.  For example,
% permute(fpt(ibox,k)),...,permute(lpt(ibox,k)) is the set of 
% all points from vesicle k contained in box ibox

neigh = zeros(Nbins,9);

%Do corners first
neigh(1,1:4) = [1 2 Nx+1 Nx+2]; 
% bottom left corner
neigh(Nx,1:4) = [Nx Nx-1 2*Nx 2*Nx-1]; 
% bottom right corner
neigh(Nbins-Nx+1,1:4) = [Nbins-Nx+1 Nbins-Nx+2 ...
    Nbins-2*Nx+1 Nbins-2*Nx+2];
% top left corner
neigh(Nbins,1:4) = [Nbins Nbins-1 Nbins-Nx Nbins-Nx-1]; 
% top right corner

for j = 2:Nx-1
  neigh(j,1:6) = j + [-1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of bottom row

for j = Nbins-Nx+2:Nbins-1
  neigh(j,1:6) = j + [-1 0 1 -Nx-1 -Nx -Nx+1];
end
% neighbors of top row

for j=Nx+1:Nx:Nbins-2*Nx+1
  neigh(j,1:6) = j + [-Nx -Nx+1 0 1 Nx Nx+1];
end
% neighbors of left column

for j=2*Nx:Nx:Nbins-Nx
  neigh(j,1:6) = j + [-Nx-1 -Nx -1 0 Nx-1 Nx];
end
% neighbors of right column

J = (Nx + 1:Nbins - Nx);
J = J(mod(J-1,Nx)~=0);
J = J(mod(J,Nx)~=0);
% J is the index of boxes that are not on the boundary
for j=J
  neigh(j,:) = j + [-Nx-1 -Nx -Nx+1 -1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of interior points
% TREE STRUCTURE IS COMPLETE


if (relate == 1 || relate == 3)
%  distSS = -ones(N1,nv1,nv1); 
%  % dist(n,k,j) is the distance of point n on vesicle k to
%  % vesicle j
%  zoneSS = -ones(N1,nv1,nv1); 
%  % near or far zone
%  nearestSS = -ones(2*N1,nv1,nv1); 
%  % nearest point using local interpolant
%  icpSS = -ones(N1,nv1,nv1); 
%  % index of closest discretization point
%  argnearSS = -ones(N1,nv1,nv1); 
%  % argument in [0,1] of local interpolant
  for k = 1:nv1
    distSS{k} = spalloc(N1,nv1,0);
    % dist(n,k,j) is the distance of point n on vesicle k to
    zoneSS{k} = spalloc(N1,nv1,0);
    % near or far zone
    nearestSS{k} = spalloc(2*N1,nv1,0);
    % nearest point using local interpolant
    icpSS{k} = spalloc(N1,nv1,0);
    % index of closest discretization point
    argnearSS{k} = spalloc(N1,nv1,0);
    % argument in [0,1] of local interpolant
  end
  % New way of representing near-singular integration structure so that
  % we can use sparse matricies.


  % begin classifying points where we are considering 
  % vesicle to vesicle relationships
  for k = 1:nv1
    boxes = unique(bin(:,k));
    % Find all boxes containing points of vesicle k
    boxes = neigh(boxes,:);
    % Look at all neighbors of boxes containing vesicle k
    boxes = unique(boxes(:));
    % Remove repetition
    boxes = boxes(boxes~=0);
    % Delete non-existent boxes that came up because of neigh

    K = [(1:k-1) (k+1:nv1)];
    for k2 = K
      istart = fpt(boxes,k2);
      iend = lpt(boxes,k2);
      istart = istart(istart ~= 0);
      iend = iend(iend ~= -1);
      % Find index of all points in neighboring boxes of vesicle k
      % that are in vesicle k2
      
      neighpts = zeros(sum(iend-istart+1),1);
      
      % Allocate space to assign possible near points
      is = 1;
      for j=1:numel(istart)
        ie = is + iend(j) - istart(j);
        neighpts(is:ie) = permute(istart(j):iend(j),k2);
        is = ie + 1;
      end
      % neighpts contains all points on vesicle k2 that are in 
      % neighboring boxes to vesicle k

      neighpts = sort(neighpts);
      % sorting should help speedup as we won't be jumping around
      % through different boxes

      n = 0;
      for i=1:numel(neighpts)
        ipt = neighpts(i);
        ibox = bin(ipt,k2);
        % box containing ipt on vesicle k2
        if (ibox ~= n)
          n = ibox;
          % Check if we need to move to a new box
          neighbors = neigh(ibox,:);
          % neighbors of this box
          neighbors = neighbors(neighbors~=0);
          % Remove non-existent neighbors
          istart = fpt(neighbors,k);
          iend = lpt(neighbors,k);
          istart = istart(istart ~= 0);
          iend = iend(iend ~= -1);
          % Find points on vesicle k in neighboring boxes
          neighpts2 = zeros(sum(iend-istart+1),1);
          is = 1;
          for j=1:numel(istart)
            ie = is + iend(j) - istart(j);
            neighpts2(is:ie) = permute(istart(j):iend(j),k);
            is = ie + 1;
          end
          % neighpts2 contains all points on vesicle k that 
          % are in neighboring box of ibox 
        end % decide if we need to switch boxes

        [d0,d0loc] = min((xsou(ipt,k2) - xsou(:,k)).^2 + ...
            (ysou(ipt,k2) - ysou(:,k)).^2);
        % Find minimum distance between ipt on vesicle k2 to
        % possible closest points on vesicle k
        d0 = sqrt(d0);
        % Save on not taking the square root on a vector but instead
        % on a single real number

        icpSS{k}(ipt,k2) = d0loc;
        if (d0 < 2*h);
          [distSS{k}(ipt,k2),nearestx,nearesty,argnearSS{k}(ipt,k2)] = ...
              vesicle1.closestPnt([xsou;ysou],xsou(ipt,k2),...
              ysou(ipt,k2),k,icpSS{k}(ipt,k2));
          nearestSS{k}(ipt,k2) = nearestx;
          nearestSS{k}(ipt+N1,k2) = nearesty;
          % Find closest point along a local interpolant using
          % Newton's method.

          if (distSS{k}(ipt,k2) < h)
            zoneSS{k}(ipt,k2) = 1;
          end
          % Point ipt of vesicle k2 is in the near zone of
          % vesicle k
        end


      end % ipt

    end % k2

  end % k

  NearSelf.dist = distSS;
  NearSelf.zone = zoneSS;
  NearSelf.nearest = nearestSS;
  NearSelf.icp = icpSS;
  NearSelf.argnear = argnearSS;
  % Store everything in the structure NearSelf.  This way it is 
  % much cleaner to pass everything around

end % relate == 1 || relate == 3



% Bin target points with respect to the source points
if (relate == 2 || relate == 3)
  N2 = vesicle2.N; % number of additional targets
  nv2 = vesicle2.nv; % number of additional vesicles
  X2 = vesicle2.X; % additional target points
  [xtar,ytar] = oc.getXY(X2);
  % separate additional target points into x and y coordinates
%  figure(2); clf
%  plot(xtar,ytar,'r.')
%  hold on
%  plot(xsou,ysou,'k.')
% DEBUG: FOR SEEING TARGET AND SOURCE POINTS IN THE TREE STRUCTURE
% WHICH CAN BE PLOTTED ABOVE

%  distST = -ones(N2,nv2,nv1); 
%  % dist(n,k,j) is the distance of point n on vesicle k to
%  % vesicle j
%  zoneST = -ones(N2,nv2,nv1); 
%  % near or far zone
%  nearestST = -ones(2*N2,nv2,nv1); 
%  % nearest point using local interpolant
%  icpST = -ones(N2,nv2,nv1); 
%  % index of closest discretization point
%  argnearST = -ones(N2,nv2,nv1); 
%  % argument in [0,1] of local interpolant
  for k = 1:nv1
    distST{k} = spalloc(N1,nv2,0);
    % dist(n,k,j) is the distance of point n on vesicle k to
    zoneST{k} = spalloc(N1,nv2,0);
    % near or far zone
    nearestST{k} = spalloc(2*N1,nv2,0);
    % nearest point using local interpolant
    icpST{k} = spalloc(N1,nv2,0);
    % index of closest discretization point
    argnearST{k} = spalloc(N1,nv2,0);
    % argument in [0,1] of local interpolant
  end
  % New way of representing near-singular integration structure so that
  % we can use sparse matricies.

  itar = ceil((xtar - xmin)/H);
  jtar = ceil((ytar - ymin)/H);
  [indx,indy] = find((itar >= 1) & (itar <= Nx) & ...
      (jtar >= 1) & (jtar <= Ny));
  % Only have to consider xx(ind),yy(ind) since all other points
  % are not contained in the box [xmin xmax] x [ymin ymax]

  for k = 1:nv1 % loop over sources
      binTar = (jtar(indx,1)-1)*Nx + itar(indx,1);

      [uBinTar,~,idxBinTar] = unique(binTar);
      
      boxesTar = neigh(uBinTar,:)';

      [uBoxesTar,~,idxUBoxTar] = unique(boxesTar);
      
      istart = zeros(numel(uBoxesTar),1);
      iend   = zeros(numel(uBoxesTar),1);
      if uBoxesTar(1) == 0
        istart(2:end) = fpt(uBoxesTar(2:end),k);
        iend(2:end)   = lpt(uBoxesTar(2:end),k);
      else
        istart = fpt(uBoxesTar,k);
        iend   = lpt(uBoxesTar,k); 
      end
      
      NonZ_idx = find(istart);
      
      d0keep  = (1e9)*ones(max(indx),1);
      
      l = 1; 
      for j = 1 : numel(NonZ_idx)

        BoI = NonZ_idx(j);

        neighpts = permute(istart(BoI):iend(BoI),k);
        
        idx_neighTBox = ceil(find(idxUBoxTar == BoI)/9);
        %uBoxesTar(BoI)
        %uBinTar(idx_neighTBox)
        
        
        %idx = cIndx(BOXES==del_uBoxesTar(j));
        for box = 1 : numel(idx_neighTBox)
          %uBinTar(idx_neighTBox(box))
          idx = find(idxBinTar == idx_neighTBox(box));
          %binTar(idx)
          
          %(jtar(indx(idx),1)-1)*Nx + itar(indx(idx),1)
          %pause
          for ii = 1 : numel(idx)
            [d0,d0loc] = min((xtar(indx(idx(ii)),1)-xsou(neighpts,k)).^2 + ...
              (ytar(indx(idx(ii)),1)-ysou(neighpts,k)).^2);
            d0 = d0^0.5;
            
            if d0 < d0keep(indx(idx(ii)))
              d0keep(indx(idx(ii)))  = d0;
              gid(l) = indx(idx(ii));
              l = l + 1;
              icpST{k}(indx(idx(ii)),1) = neighpts(d0loc);  
            end
           end % for ii
         end % for box
       end % for j
       

       gid = unique(gid);
       for i = 1 : numel(gid)
         IDoI = gid(i);
         d0 = d0keep(IDoI);       
      
          if d0 < 2*h
            [distST{k}(IDoI,1),nearestx,nearesty,argnearST{k}(IDoI,1)] = ...
              vesicle1.closestPnt([xsou;ysou],xtar(IDoI,1),...
              ytar(IDoI,1),k,icpST{k}(IDoI,1));
           
            nearestST{k}(IDoI,1) = nearestx;
            nearestST{k}(IDoI+N2,1) = nearesty;

            if distST{k}(IDoI,1) < h
              zoneST{k}(IDoI,1) = 1;
            end

          end % d0 < 2*h
        end % for ll

          %end % for ii

         %end % for box

      %end % for j

  end % k

  NearOther.dist = distST;
  NearOther.zone = zoneST;
  NearOther.nearest = nearestST;
  NearOther.icp = icpST;
  NearOther.argnear = argnearST;
  % store near-singluar integration requirements in structure NearOther

end % relate == 2 || relate == 3
end % getZoneNew


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist,nearestx,nearesty,theta] = ...
    closestPnt(o,X,xtar,ytar,k,icp)
% [dist,nearestx,nearesty,theta] = closestPnt(X,xtar,ytar,k,icp)
% computes the closest point on vesicle k to (xtar,ytar)
% using a Lagrange interpolant.  icp is the index of the closest
% point on the discrete mesh which is used as an initial guess

N = size(X,1)/2; % Number of points per vesicle
A = poten.lagrangeInterp;
interpOrder = size(A,1);
% need interpolation matrix and its size

p = ceil((interpOrder+1)/2);
% Accommodate for either an even or odd number of interpolation points
pn = mod((icp-p+1:icp-p+interpOrder)' - 1,N) + 1;
% band of points around icp.  The -1,+1 combination sets index
% 0 to N as required by the code

px = A*X(pn,k); % polynomial interpolant of x-coordinate
py = A*X(pn+N,k); % polynomial interpolant of y-coordinate
Dpx = px(1:end-1).*(interpOrder-1:-1:1)';
Dpy = py(1:end-1).*(interpOrder-1:-1:1)';
D2px = Dpx(1:end-1).*(interpOrder-2:-1:1)';
D2py = Dpy(1:end-1).*(interpOrder-2:-1:1)';
% To do Newton's method, need two derivatives

theta = 1/2;
% midpoint is a good initial guess
for newton = 1:1
  zx = filter(1,[1 -theta],px);
  zx = zx(end);
  zy = filter(1,[1 -theta],py);
  zy = zy(end);
  Dzx = filter(1,[1 -theta],Dpx);
  Dzx = Dzx(end);
  Dzy = filter(1,[1 -theta],Dpy);
  Dzy = Dzy(end);
  D2zx = filter(1,[1 -theta],D2px);
  D2zx = D2zx(end);
  D2zy = filter(1,[1 -theta],D2py);
  D2zy = D2zy(end);
  % Using filter is the same as polyval, but it is much
  % faster when only requiring a single polyval such as here.

  newtonNum = (zx-xtar)*Dzx + (zy-ytar)*Dzy;
  % numerator of Newton's method
  newtonDen = (zx-xtar)*D2zx + (zy-ytar)*D2zy + ...
      Dzx^2 + Dzy^2;
  % denominator of Newton's method
  theta = theta - newtonNum/newtonDen;
  % one step of Newton's method
end
% Do a few (no more than 3) Newton iterations

nearestx = filter(1,[1,-theta],px);
nearestx = nearestx(end);
nearesty = filter(1,[1,-theta],py);
nearesty = nearesty(end);
dist = sqrt((nearestx - xtar)^2 + (nearesty - ytar)^2);
% Compute nearest point and its distance from the target point

end % closestPnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [icollisionVes,icollisionWall] = ...
    collision(vesicle,walls,NearV2V,NearV2W,fmm,op)
% [icollisionVes,icollisionWall] = ...
%    collision(vesicle,walls,NearV2V,NearV2W,fmm,op)
% uses double-layer potential for Laplace's equation to determing if
% vesicles have crossed (stored in icollisionVes) or if solid wall 
% point is inside a vesicle (icollisionWall).
% Can be done with or without the fmm.  vesicle stores all the
% vesicles, walls stores the solid walls (if they exist), NearV2V
% and NearV2W are near-singular integration structures, fmm and near
% are flags for using the FMM and near-singular integration

oc = curve;
[x,y] = oc.getXY(vesicle.X);

f = [ones(vesicle.N,vesicle.nv);zeros(vesicle.N,vesicle.nv)];
% Density function is constant.  Pad second half of it with zero


if ~fmm
  kernel = @op.exactLaplaceDL;
else
  kernel = @op.exactLaplaceDLfmm;
end
% kernel for laplace's double layer potential.  Only difference
% is if FMM is used or not


DLP = @(X) zeros(2*size(X,1),size(X,2));
% can cheat here because we know that the double-layer
% potential applied to our function f will always be 0
% This won't work if we are considering density functions
% that are not one everywhere

Fdlp = op.nearSingInt(vesicle,f,DLP,[],...
      NearV2V,kernel,kernel,vesicle,true,false);
% Apply the double-layer laplace potential with constant boundary
% condition.  Skips self-vesicle term.  This means that if the vesicles have
% crossed, the Fdlp will achieve a value of 1.  If they have not crossed, it
% will be 0
if ~isempty(walls)
  FDLPwall = op.nearSingInt(vesicle,f,DLP,[],...
        NearV2W,kernel,kernel,walls,false,false);
end
% Evaulate the laplace double-layer potential due to all vesicles with the
% solid walls acting as the targets.  Will be zero at all the targets if
% none of the vesicles have crossed the solid walls

Fdlp = Fdlp(1:vesicle.N,:);
% this layer-potential is not vector valued, so we throw away the second half.

bufferVes = 1e-4;
bufferWall = 1e-4;
% can't set buffer too big because near singular integration does not
% assign a value of 1 when near points cross.  This is because I did not
% swtich the direction of the normal for this case.  So, the lagrange
% interpolation points may look something like [0 1 1 1 1 ...] instead
% of [1 1 1 1 1 ...].  The local interpolant is affected by this and
% a value somewhere between 0 and 1 will be returned
icollisionVes = any(abs(Fdlp(:)) > bufferVes);

if ~isempty(walls)
  FDLPwall = FDLPwall(1:walls.N,:);
  icollisionWall = any(abs(FDLPwall(:)) > bufferWall);
else
  icollisionWall = false;
end

end % collision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [icollisionVes,collidingIdcs,safeIdcs] = ...
    collisionWVesOutput(vesicle,NearV2V,fmm,op)

oc = curve;
[x,y] = oc.getXY(vesicle.X);

f = [ones(vesicle.N,vesicle.nv);zeros(vesicle.N,vesicle.nv)];
% Density function is constant.  Pad second half of it with zero

if ~fmm
  kernel = @op.exactLaplaceDL;
else
  kernel = @op.exactLaplaceDLfmm;
end
% kernel for laplace's double layer potential.  Only difference
% is if FMM is used or not


DLP = @(X) zeros(2*size(X,1),size(X,2));
% can cheat here because we know that the double-layer
% potential applied to our function f will always be 0
% This won't work if we are considering density functions
% that are not one everywhere

Fdlp = op.nearSingInt(vesicle,f,DLP,[],...
      NearV2V,kernel,kernel,vesicle,true,false);
% Apply the double-layer laplace potential with constant boundary
% condition.  Skips self-vesicle term.  This means that if the vesicles have
% crossed, the Fdlp will achieve a value of 1.  If they have not crossed, it
% will be 0

Fdlp = Fdlp(1:vesicle.N,:);
% this layer-potential is not vector valued, so we throw away the second half.

bufferVes = 2e-3;
% can't set buffer too big because near singular integration does not
% assign a value of 1 when near points cross.  This is because I did not
% swtich the direction of the normal for this case.  So, the lagrange
% interpolation points may look something like [0 1 1 1 1 ...] instead
% of [1 1 1 1 1 ...].  The local interpolant is affected by this and
% a value somewhere between 0 and 1 will be returned
icollisionVes = any(abs(Fdlp(:)) > bufferVes);

% retrieve the indeces of colliding vesicles
safeIdcs = [];

ids = find(abs(Fdlp)>bufferVes);
collidingIdcs = ceil(ids/vesicle.N);
collidingIdcs = collidingIdcs(:);
collidingIdcs = unique(collidingIdcs);
for k = 1 : vesicle.nv
  if any(collidingIdcs) ~= k
    safeIdcs = [safeIdcs; k];
  end
end


% if there is a single vesicle colliding (!), that's a numerical error and
% there is no collision
if numel(collidingIdcs) == 1
  collidingIdcs = [];
  icollisionVes = 0;
else
  % if there are colliding vesicles far from each other, then there is no collision  
  % numerical error again
  % centers of vesicles
  centx = mean(x,1); centy = mean(y,1);  
  reallyCollIdcs = [];  
  for i = 1 : numel(collidingIdcs)
    dist2rest = sqrt((centx(collidingIdcs(i))-centx(collidingIdcs)).^2 + ...
      (centy(collidingIdcs(i))-centy(collidingIdcs)).^2);
    [~,sortId] = sort(dist2rest);     
    if dist2rest(sortId(2))<vesicle.length/2
       reallyCollIdcs = [reallyCollIdcs;collidingIdcs(i)];   
    end
  end
  collidingIdcs = reallyCollIdcs;
  if isempty(reallyCollIdcs)
    icollisionVes = 0;
  end  
end



end % collision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [icollisionVes,collidingIdcs,safeIdcs] = ...
    collisionWwallOutput(vesicle,walls,NearV2W,fmm,op)

oc = curve;
[x,y] = oc.getXY(vesicle.X);

f = [ones(vesicle.N,vesicle.nv);zeros(vesicle.N,vesicle.nv)];
% Density function is constant.  Pad second half of it with zero

if ~fmm
  kernel = @op.exactLaplaceDL;
else
  kernel = @op.exactLaplaceDLfmm;
end
% kernel for laplace's double layer potential.  Only difference
% is if FMM is used or not


DLP = @(X) zeros(2*size(X,1),size(X,2));
% can cheat here because we know that the double-layer
% potential applied to our function f will always be 0
% This won't work if we are considering density functions
% that are not one everywhere

Fdlp = op.nearSingInt(vesicle,f,DLP,[],...
      NearV2W,kernel,kernel,walls,false,false);
% Apply the double-layer laplace potential with constant boundary
% condition.  Skips self-vesicle term.  This means that if the vesicles have
% crossed, the Fdlp will achieve a value of 1.  If they have not crossed, it
% will be 0

Fdlp = Fdlp(1:vesicle.N,:);
% this layer-potential is not vector valued, so we throw away the second half.

bufferVes = 2e-3;
% can't set buffer too big because near singular integration does not
% assign a value of 1 when near points cross.  This is because I did not
% swtich the direction of the normal for this case.  So, the lagrange
% interpolation points may look something like [0 1 1 1 1 ...] instead
% of [1 1 1 1 1 ...].  The local interpolant is affected by this and
% a value somewhere between 0 and 1 will be returned
icollisionVes = any(abs(Fdlp(:)) > bufferVes);

% retrieve the indeces of colliding vesicles
safeIdcs = [];

ids = find(abs(Fdlp)>bufferVes);
collidingIdcs = ceil(ids/vesicle.N);
collidingIdcs = collidingIdcs(:);
collidingIdcs = unique(collidingIdcs);
for k = 1 : vesicle.nv
  if any(collidingIdcs) ~= k
    safeIdcs = [safeIdcs; k];
  end
end


% if there is a single vesicle colliding (!), that's a numerical error and
% there is no collision
if numel(collidingIdcs) == 1
  collidingIdcs = [];
  icollisionVes = 0;
else
  % if there are colliding vesicles far from each other, then there is no collision  
  % numerical error again
  % centers of vesicles
  centx = mean(x,1); centy = mean(y,1);  
  reallyCollIdcs = [];  
  for i = 1 : numel(collidingIdcs)
    dist2rest = sqrt((centx(collidingIdcs(i))-centx(collidingIdcs)).^2 + ...
      (centy(collidingIdcs(i))-centy(collidingIdcs)).^2);
    [~,sortId] = sort(dist2rest);     
    if dist2rest(sortId(2))<vesicle.length/2
       reallyCollIdcs = [reallyCollIdcs;collidingIdcs(i)];   
    end
  end
  collidingIdcs = reallyCollIdcs;
  if isempty(reallyCollIdcs)
    icollisionVes = 0;
  end  
end



end % collision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InOut = sortPts(vesicle,Xtar,fmm,NearV2T,op)
% InOut = sortPts(vesicle,Xtar,fmm,nearV2T) determines if the set of
% points in Xtar are inside or outside of a vesicle configuration

density = [ones(vesicle.N,vesicle.nv);...
           zeros(vesicle.N,vesicle.nv)];
% density function that is used to check if a point is inside or
% outside

tracers = capsules(Xtar,[],[],0,0,0);
if nargin == 3 
  [~,NearV2T] = vesicle.getZone(tracers,2);
end

if nargin ~= 5
  op = poten(vesicle.N);
end

if ~fmm
  kernel = @op.exactLaplaceDL;
else
  kernel = @op.exactLaplaceDLfmm;
end
% kernel for Laplace's double layer potential

DLP = @(X) zeros(2*size(X,1),size(X,2));
% can cheat here because we know that the double-layer potential
% applied to the constant density function will always return zero

kernelDirect = kernel;
InOut = op.nearSingInt(vesicle,density,DLP,[],...
    NearV2T,kernel,kernelDirect,tracers,false,false);

InOut = InOut(1:end/2);
% only care about the first half as we are working with laplace as
% opposed to the vector-valued stokes layer potentials

thresh = 1e-1;
InOut(abs(InOut) > thresh) = 1;
InOut(abs(InOut) <= thresh) = 0;
% for points inside a vesicle, but close to its boundary, we will be
% interpolating the function [0;1;1;1;1...] close to the initial
% point.  Therefore, the value returned by the DLP will be some value
% between 0 and 1, but it can be quite close to 0.  So, we threshold

end % sortPts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function press = pressure(vesicle,f,RS,pressTar,fmm,LP)
% press = pressure(vesicle,f,RS,pressTar,fmm,LP) computes the pressure
% due to vesicle at the locations pressTar with near-singular
% integration with traction jump f and Rotlets stored in RS.  LP is
% either SLP or DLP and tells this routine if it needs to evaluate the
% pressure due to the single-layer potential or the double-layer
% potential.  Note that vesicle may correspond to solid walls rather
% than vesicles

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles
op = poten(N);
Ntar = pressTar.N;
% number of target points where we want to compute the pressure
% and stress

[~,NearV2T] = vesicle.getZone(pressTar,2);
% build near-singular integration structure for vesicle
% to pressure target points
InOutFlag = find(vesicle.sortPts(pressTar.X,fmm,NearV2T) == 1);
%figure(1); clf; hold on
%plot(vesicle.X(1:end/2,:),vesicle.X(end/2+1:end,:),'r');
%oc = curve;
%[xx,yy] = oc.getXY(pressTar.X);
%plot(xx,yy,'b.');
%plot(xx(InOutFlag),yy(InOutFlag),'k.');
% DEBUG: MAKE SURE POINTS INSIDE AND OUTSIDE ARE LABELLED CORRECTLY

%fprintf('PRESSURE\n')
if strcmp(LP,'SLP')
  PdiagIn = op.pressSLmatrix(vesicle);
  PdiagOut = PdiagIn;
  % use odd-even integration to evaluate pressure due to
  % each vesicle independent of the others on its boundary
  for k = 1:nv
    PdiagOut(1:N,1:N,k) = PdiagOut(1:N,1:N,k) + ...
        1/2*diag(vesicle.xt(N+1:2*N,k));
    PdiagOut(1:N,N+1:2*N,k) = PdiagOut(1:N,N+1:2*N,k) + ...
        1/2*diag(-vesicle.xt(1:N,k));
    PdiagIn(1:N,1:N,k) = PdiagIn(1:N,1:N,k) - ...
        1/2*diag(vesicle.xt(N+1:2*N,k));
    PdiagIn(1:N,N+1:2*N,k) = PdiagIn(1:N,N+1:2*N,k) - ...
        1/2*diag(-vesicle.xt(1:N,k));
  end
  % Jump term that comes from pressure of single-layer potential is the
  % dot product of the traction jump with the normal vector multiplied
  % by 1/2
  Pdiag = @(X) op.exactPressureSLdiag(vesicle,[PdiagOut;PdiagIn],X);
  % nearSingInt assumes that input and output are vector-valued Take
  % adavntage of the fact that nearSingInt only works with
  % vector-valued densities by passing in two different jumps---one for
  % the limit from the inside and one from the outside.

  if ~fmm
    kernel = @op.exactPressSL;
  else
    kernel = @op.exactPressSLfmm;
  end
  kernelDirect = @op.exactPressSL;
elseif strcmp(LP,'DLP')
  oc = curve;
  PdiagIn = op.pressDLmatrix(vesicle);
  PdiagOut = PdiagIn;
  Deriv = fft1.D1(N);
  % spectral differentiation matrix
  for k = 1:nv
    [tanx,tany] = oc.getXY(vesicle.xt(:,k));
    % tangent vector
    sa = vesicle.sa(:,k);
    % Jacobian

    jump = -[diag(tanx./sa) diag(tany./sa)] * ...
        [Deriv zeros(N); zeros(N) Deriv];
    PdiagIn(:,:,k) = PdiagIn(:,:,k) - jump;
    PdiagOut(:,:,k) = PdiagOut(:,:,k) + jump;
    % add in the jump term.  Jump term is negative because target
    % points are interior to the solid walls
  end
  % Jump term that comes from pressure of double-layer potential
  % is the dot product of the tangential derivative of the traction 
  % jump with the tangent.  Assuming that all the pressure target
  % points are inside the physical domain so there is no use to
  % consider the limiting value from the other side
  Pdiag = @(X) op.exactPressureDLdiag(vesicle,[PdiagOut;PdiagIn],X);
  % nearSingInt assumes that input and output are vector-valued Take
  % adavntage of the fact that nearSingInt only works with
  % vector-valued densities by passing in two different jumps---one for
  % the limit from the inside and one from the outside.

  kernel = @op.exactPressDL;
  kernelDirect = @op.exactPressDL;
end

pressT = op.nearSingInt(vesicle,f,Pdiag,[],...
    NearV2T,kernel,kernelDirect,pressTar,false,false);
% compute the pressure of the single- or double-layer 
% potential using near-singular integration.  First row 
% corresponds to using the limiting value for the pressure 
% exterior of the vesicle and the second row corresponds to 
% using the limiting value for the pressure interior of the 
% vesicle
press = pressT(1:Ntar);
press1 = press;
press(InOutFlag) = pressT(Ntar+InOutFlag);
% At the interior points, take the second row

if ~isempty(RS)
  for k = 2:vesicle.nv
    cx = vesicle.center(1,k);
    cy = vesicle.center(2,k);
    % center of interior solid wall k
    xi1 = RS(3*(k-2) + 1,k);
    xi2 = RS(3*(k-2) + 2,k);
    % first and second components of the stokeslet on 
    % interior solid wall k 
    [x,y] = oc.getXY(pressTar.X);
    rx = x - cx;
    ry = y - cy;
    press = press +  2./(rx.^2 + ry.^2).*(rx*xi1 + ry*xi2);
    % add in pressure of stokeslet term
    % rotlet has a vanishing pressure gradient
  end
end

end % pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stress1,stress2] = stressSelfTensor(vesicle,f,LP,op)
    
N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles

% decompose tangent and normal vectors into their x and y
% coordinates
tangent = vesicle.xt;
oc = curve;
[tanx,tany] = oc.getXY(tangent);
nx = tany;
ny = -tanx;

if strcmp(LP,'SLP')
  
  % use odd-even integration to evaluate stress due to
  % each vesicle independent of the others on its boundary
  % S1diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [1;0] evaluated on the boundary of
  % the vesicle
  % S2diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [0;1] evaluated on the boundary of
  % the vesicle  
  [S1diagIn,S2diagIn] = op.stressSLmatrix(vesicle);
  
  for k = 1:nv
    % add in (1,1) jump component to stress applied to [1;0]   
    jump = +0.5*nx(:,k) + tanx(:,k).^2.*tany(:,k);
    S1diagIn(1:N,1:N,k) = S1diagIn(1:N,1:N,k) + diag(jump);
    
    % add in (2,1) jump component to stress applied to [1;0] 
    jump = +0.5*ny(:,k) + tanx(:,k).*tany(:,k).^2;
    S1diagIn(N+1:2*N,1:N,k) = S1diagIn(N+1:2*N,1:N,k) + diag(jump);
    
    % add in (1,2) jump component to stress applied to [1;0]
    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tanx(:,k);
    S1diagIn(1:N,N+1:2*N,k) = S1diagIn(1:N,N+1:2*N,k) + diag(jump);
     

    % add in (2,2) jump component to stress applied to [1;0]
    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tany(:,k);
    S1diagIn(N+1:2*N,N+1:2*N,k) = ...
        S1diagIn(N+1:2*N,N+1:2*N,k) + diag(jump);
     

    % add in (1,1) jump component to stress applied to [0;1]
    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tanx(:,k);
    S2diagIn(1:N,1:N,k) = S2diagIn(1:N,1:N,k) + diag(jump);
     
    % add in (2,1) jump component to stress applied to [0;1]
    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tany(:,k);
    S2diagIn(N+1:2*N,1:N,k) = S2diagIn(N+1:2*N,1:N,k) + diag(jump);
     

    % add in (1,2) jump component to stress applied to [0;1]
    jump = +0.5*nx(:,k) - tanx(:,k).^2.*tany(:,k);
    S2diagIn(1:N,N+1:2*N,k) = S2diagIn(1:N,N+1:2*N,k) + diag(jump);
     

    % add in (2,2) jump component to stress applied to [0;1]
    jump = +0.5*ny(:,k) - tanx(:,k).*tany(:,k).^2;
    S2diagIn(N+1:2*N,N+1:2*N,k) = ...
        S2diagIn(N+1:2*N,N+1:2*N,k) + diag(jump);
  end
  % Jump term that comes from stress of single-layer potential
  % is the 
  S1diagInFn = @(X) op.exactStressSLdiag(vesicle,S1diagIn,X);
  S2diagInFn = @(X) op.exactStressSLdiag(vesicle,S2diagIn,X);


elseif strcmp(LP,'DLP')
  % use odd-even integration to evaluate stress due to
  % each vesicle independent of the others on its boundary
  % S1diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [1;0] evaluated on the boundary of
  % the vesicle
  % S2diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [0;1] evaluated on the boundary of
  % the vesicle  
  [S1diagIn,S2diagIn] = op.stressDLmatrix(vesicle);
  
  % Fourier differentiation matrix
  Deriv = fft1.D1(N);
  

  for k = 1:nv
    % tangent vector  
    [tanx,tany] = oc.getXY(vesicle.xt(:,k));
    % Jacobian
    sa = vesicle.sa(:,k);
    
    % add in (1,1) jump component to stress applied to [1;0] 
    jump = diag(1./sa.*(1+tanx.^2-tany.^2).*tanx)*Deriv;
    S1diagIn(1:N,1:N,k) = S1diagIn(1:N,1:N,k) - jump;
    
    % add in (1,2) jump component to stress applied to [1;0] 
    jump = diag(1./sa.*(1+tanx.^2-tany.^2).*tany)*Deriv;
    S1diagIn(1:N,N+1:2*N,k) = S1diagIn(1:N,N+1:2*N,k) - jump;
    
    % add in (2,1) jump component to stress applied to [1;0] 
    jump = diag(1./sa.*(2*tanx.*tany).*tanx)*Deriv;
    S1diagIn(N+1:2*N,1:N,k) = S1diagIn(N+1:2*N,1:N,k) - jump;
    
    % add in (2,2) jump component to stress applied to [1;0]
    jump = diag(1./sa.*(2*tanx.*tany).*tany)*Deriv;
    S1diagIn(N+1:2*N,N+1:2*N,k) = ...
        S1diagIn(N+1:2*N,N+1:2*N,k) - jump;
     

    % add in (1,1) jump component to stress applied to [0;1]
    jump = diag(1./sa.*(2*tanx.*tany).*tanx)*Deriv;
    S2diagIn(1:N,1:N,k) = S2diagIn(1:N,1:N,k) - jump;
     
    % add in (1,2) jump component to stress applied to [0;1]
    jump = diag(1./sa.*(2*tanx.*tany).*tany)*Deriv;
    S2diagIn(1:N,N+1:2*N,k) = S2diagIn(1:N,N+1:2*N,k) - jump;
     
    % add in (2,1) jump component to stress applied to [0;1] 
    jump = diag(1./sa.*(1-tanx.^2+tany.^2).*tanx)*Deriv;
    S2diagIn(N+1:2*N,1:N,k) = S2diagIn(N+1:2*N,1:N,k) - jump;
    
    % add in (2,2) jump component to stress applied to [0;1]
    jump = diag(1./sa.*(1-tanx.^2+tany.^2).*tany)*Deriv;
    S2diagIn(N+1:2*N,N+1:2*N,k) = ...
        S2diagIn(N+1:2*N,N+1:2*N,k) - jump;
     
  end
  % Jump term that comes from stress of double-layer potential
  S1diagInFn = @(X) op.exactStressDLdiag(vesicle,S1diagIn,X);
  S2diagInFn = @(X) op.exactStressDLdiag(vesicle,S2diagIn,X);
  
end

% Points on the membrane are considered outside
stress1 = S1diagInFn(f);
stress2 = S2diagInFn(f);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stress1,stress2] = stressTensor(vesicle,f,...
    RS,stressTar,fmm,LP,op,NearV2T)
% For the details of computing stressTensor, refer to Quaife & Biros, JCP,
% 2014.
% [stress1 stress2] = stressTensor(vesicle,f,RS,stressTar,fmm,LP)
% computes the stress tensor due to vesicle at the locations stressTar
% with or without near-singular integration with traction jump f.
% Returns in two components.  stress1 is the stress tensor applied to
% [1;0] and stress2 is the stress tensor applied to [0;1].  LP is
% either SLP or DLP and tells this routine if it needs to evaluate the
% stress due to the single-layer potential or the double-layer
% potential.  Note that vesicle may correspond to solid walls rather
% than vesicles.  If doing solid walls, there is a stress due to the
% rotlet and stokeslet terms that needs to be added in

N = vesicle.N; % points per vesicle
nv = vesicle.nv; % number of vesicles

% number of points where we want to evaluate the stress
Ntar = stressTar.N*stressTar.nv;

% decompose tangent and normal vectors into their x and y
% coordinates
tangent = vesicle.xt;
oc = curve;
[tanx,tany] = oc.getXY(tangent);
nx = tany;
ny = -tanx;

% build near-singular integration structure for vesicle
% to pressure target points
if isempty(NearV2T)
  [~,NearV2T] = vesicle.getZone(stressTar,2);
end

% InOutFlag = find(vesicle.sortPts(stressTar.X,fmm,NearV2T) == 1);

InOutFlag = vesicle.sortPts(stressTar.X,fmm,NearV2T)==1;

if strcmp(LP,'SLP')
  
  % use odd-even integration to evaluate stress due to
  % each vesicle independent of the others on its boundary
  % S1diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [1;0] evaluated on the boundary of
  % the vesicle
  % S2diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [0;1] evaluated on the boundary of
  % the vesicle  
  [S1diagIn,S2diagIn] = op.stressSLmatrix(vesicle);
  
  % Use a second copy for the jump from the outside
  S1diagOut = S1diagIn;
  S2diagOut = S2diagIn;
  
  for k = 1:nv
    % add in (1,1) jump component to stress applied to [1;0]   
    jump = +0.5*nx(:,k) + tanx(:,k).^2.*tany(:,k);
    S1diagIn(1:N,1:N,k) = S1diagIn(1:N,1:N,k) + diag(jump);
    S1diagOut(1:N,1:N,k) = S1diagOut(1:N,1:N,k) - diag(jump);
    
    % add in (2,1) jump component to stress applied to [1;0] 
    jump = +0.5*ny(:,k) + tanx(:,k).*tany(:,k).^2;
    S1diagIn(N+1:2*N,1:N,k) = S1diagIn(N+1:2*N,1:N,k) + diag(jump);
    S1diagOut(N+1:2*N,1:N,k) = S1diagOut(N+1:2*N,1:N,k) - diag(jump);
    
    % add in (1,2) jump component to stress applied to [1;0]
    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tanx(:,k);
    S1diagIn(1:N,N+1:2*N,k) = S1diagIn(1:N,N+1:2*N,k) + diag(jump);
    S1diagOut(1:N,N+1:2*N,k) = S1diagOut(1:N,N+1:2*N,k) - diag(jump);
     

    % add in (2,2) jump component to stress applied to [1;0]
    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tany(:,k);
    S1diagIn(N+1:2*N,N+1:2*N,k) = ...
        S1diagIn(N+1:2*N,N+1:2*N,k) + diag(jump);
    S1diagOut(N+1:2*N,N+1:2*N,k) = ...
        S1diagOut(N+1:2*N,N+1:2*N,k) - diag(jump);
     

    % add in (1,1) jump component to stress applied to [0;1]
    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tanx(:,k);
    S2diagIn(1:N,1:N,k) = S2diagIn(1:N,1:N,k) + diag(jump);
    S2diagOut(1:N,1:N,k) = S2diagOut(1:N,1:N,k) - diag(jump);
     
    % add in (2,1) jump component to stress applied to [0;1]
    jump = +0.5*(tany(:,k).^2 - tanx(:,k).^2).*tany(:,k);
    S2diagIn(N+1:2*N,1:N,k) = S2diagIn(N+1:2*N,1:N,k) + diag(jump);
    S2diagOut(N+1:2*N,1:N,k) = S2diagOut(N+1:2*N,1:N,k) - diag(jump);
     

    % add in (1,2) jump component to stress applied to [0;1]
    jump = +0.5*nx(:,k) - tanx(:,k).^2.*tany(:,k);
    S2diagIn(1:N,N+1:2*N,k) = S2diagIn(1:N,N+1:2*N,k) + diag(jump);
    S2diagOut(1:N,N+1:2*N,k) = S2diagOut(1:N,N+1:2*N,k) - diag(jump);
     

    % add in (2,2) jump component to stress applied to [0;1]
    jump = +0.5*ny(:,k) - tanx(:,k).*tany(:,k).^2;
    S2diagIn(N+1:2*N,N+1:2*N,k) = ...
        S2diagIn(N+1:2*N,N+1:2*N,k) + diag(jump);
    S2diagOut(N+1:2*N,N+1:2*N,k) = ...
        S2diagOut(N+1:2*N,N+1:2*N,k) - diag(jump);
  end
  % Jump term that comes from stress of single-layer potential
  % is the 
  S1diagInFn = @(X) op.exactStressSLdiag(vesicle,S1diagIn,X);
  S2diagInFn = @(X) op.exactStressSLdiag(vesicle,S2diagIn,X);
  S1diagOutFn = @(X) op.exactStressSLdiag(vesicle,S1diagOut,X);
  S2diagOutFn = @(X) op.exactStressSLdiag(vesicle,S2diagOut,X);
  % nearSingInt assumes that input and output are vector-valued
  kernel1 = @op.exactStressSL1;
  kernel2 = @op.exactStressSL2;

elseif strcmp(LP,'DLP')
  % use odd-even integration to evaluate stress due to
  % each vesicle independent of the others on its boundary
  % S1diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [1;0] evaluated on the boundary of
  % the vesicle
  % S2diagIn is the matrix that takes the traction jump and returns
  % the stress tensor applied to [0;1] evaluated on the boundary of
  % the vesicle  
  [S1diagIn,S2diagIn] = op.stressDLmatrix(vesicle);
  
  % Use a second copy for the jump from the outside
  S1diagOut = S1diagIn;
  S2diagOut = S2diagIn;
  
  % Fourier differentiation matrix
  Deriv = fft1.D1(N);
  

  for k = 1:nv
    % tangent vector  
    [tanx,tany] = oc.getXY(vesicle.xt(:,k));
    % Jacobian
    sa = vesicle.sa(:,k);
    
    % add in (1,1) jump component to stress applied to [1;0] 
    jump = diag(1./sa.*(1+tanx.^2-tany.^2).*tanx)*Deriv;
    S1diagIn(1:N,1:N,k) = S1diagIn(1:N,1:N,k) - jump;
    S1diagOut(1:N,1:N,k) = S1diagOut(1:N,1:N,k) + jump;
    
    % add in (1,2) jump component to stress applied to [1;0] 
    jump = diag(1./sa.*(1+tanx.^2-tany.^2).*tany)*Deriv;
    S1diagIn(1:N,N+1:2*N,k) = S1diagIn(1:N,N+1:2*N,k) - jump;
    S1diagOut(1:N,N+1:2*N,k) = S1diagOut(1:N,N+1:2*N,k) + jump;
    
    % add in (2,1) jump component to stress applied to [1;0] 
    jump = diag(1./sa.*(2*tanx.*tany).*tanx)*Deriv;
    S1diagIn(N+1:2*N,1:N,k) = S1diagIn(N+1:2*N,1:N,k) - jump;
    S1diagOut(N+1:2*N,1:N,k) = S1diagOut(N+1:2*N,1:N,k) + jump;
    
    % add in (2,2) jump component to stress applied to [1;0]
    jump = diag(1./sa.*(2*tanx.*tany).*tany)*Deriv;
    S1diagIn(N+1:2*N,N+1:2*N,k) = ...
        S1diagIn(N+1:2*N,N+1:2*N,k) - jump;
    S1diagOut(N+1:2*N,N+1:2*N,k) = ...
        S1diagOut(N+1:2*N,N+1:2*N,k) + jump;
     

    % add in (1,1) jump component to stress applied to [0;1]
    jump = diag(1./sa.*(2*tanx.*tany).*tanx)*Deriv;
    S2diagIn(1:N,1:N,k) = S2diagIn(1:N,1:N,k) - jump;
    S2diagOut(1:N,1:N,k) = S2diagOut(1:N,1:N,k) + jump;
     
    % add in (1,2) jump component to stress applied to [0;1]
    jump = diag(1./sa.*(2*tanx.*tany).*tany)*Deriv;
    S2diagIn(1:N,N+1:2*N,k) = S2diagIn(1:N,N+1:2*N,k) - jump;
    S2diagOut(1:N,N+1:2*N,k) = S2diagOut(1:N,N+1:2*N,k) + jump;
     
    % add in (2,1) jump component to stress applied to [0;1] 
    jump = diag(1./sa.*(1-tanx.^2+tany.^2).*tanx)*Deriv;
    S2diagIn(N+1:2*N,1:N,k) = S2diagIn(N+1:2*N,1:N,k) - jump;
    S2diagOut(N+1:2*N,1:N,k) = S2diagOut(N+1:2*N,1:N,k) + jump;
    
    % add in (2,2) jump component to stress applied to [0;1]
    jump = diag(1./sa.*(1-tanx.^2+tany.^2).*tany)*Deriv;
    S2diagIn(N+1:2*N,N+1:2*N,k) = ...
        S2diagIn(N+1:2*N,N+1:2*N,k) - jump;
    S2diagOut(N+1:2*N,N+1:2*N,k) = ...
        S2diagOut(N+1:2*N,N+1:2*N,k) + jump;
     
  end
  % Jump term that comes from stress of double-layer potential
  S1diagInFn = @(X) op.exactStressDLdiag(vesicle,S1diagIn,X);
  S2diagInFn = @(X) op.exactStressDLdiag(vesicle,S2diagIn,X);
  S1diagOutFn = @(X) op.exactStressDLdiag(vesicle,S1diagOut,X);
  S2diagOutFn = @(X) op.exactStressDLdiag(vesicle,S2diagOut,X);
  % nearSingInt assumes that input and output are vector-valued

  kernel1 = @op.exactStressDL1;
  kernel2 = @op.exactStressDL2;
end
% Have built all single- or double-layer potentials on boundary so that
% we can compute diagonal value to use near-singular integration

% compute stress1 for all points
stress1 = op.nearSingInt(vesicle,f,S1diagInFn,[],...
      NearV2T,kernel1,kernel1,stressTar,false,false);
% if there are points outside then correct stress at those points
if any(InOutFlag==0)
  stressT = op.nearSingInt(vesicle,f,S1diagOutFn,[],...
      NearV2T,kernel1,kernel1,stressTar,false,false);      
  stress1(InOutFlag) = stressT(InOutFlag);
  stress1(Ntar+InOutFlag) = stressT(Ntar+InOutFlag);
end

% compute stress2 for all points
stress2 = op.nearSingInt(vesicle,f,S2diagInFn,[],...
      NearV2T,kernel2,kernel2,stressTar,false,false);
% if there are points outside then correct stress at those points
if any(InOutFlag==0)
  stressT = op.nearSingInt(vesicle,f,S2diagOutFn,[],...
      NearV2T,kernel2,kernel2,stressTar,false,false);
  stress2(InOutFlag) = stressT(InOutFlag);
  stress2(Ntar+InOutFlag) = stressT(Ntar+InOutFlag);
end

% if numel(InOutFlag) > 0
% %    fprintf('Stress1 Inside\n')
%  stress1 = op.nearSingInt(vesicle,f,S1diagInFn,...
%       NearV2T,kernel1,kernel1,stressTar,false,false);
% else
%   stress1 = zeros(2*Ntar,1);
% end
% if Ntar - numel(InOutFlag) > 0
% %    fprintf('Stress1 Outside\n')
%   
%   stressT = op.nearSingInt(vesicle,f,S1diagOutFn,...
%       NearV2T,kernel1,kernel1,stressTar,false,false);
%   
%   stress1(InOutFlag) = stressT(InOutFlag);
%   stress1(Ntar+InOutFlag) = stressT(Ntar+InOutFlag);
%   % correct stress at exterior points
% end
% Use near-singular integration to compute first component of
% the stress due to the single-layer potential

% if numel(InOutFlag) > 0
% %    fprintf('Stress2 Inside\n')
%   stress2 = op.nearSingInt(vesicle,f,S2diagInFn,...
%       NearV2T,kernel2,kernel2,stressTar,false,false);
% else
%   stress2 = zeros(2*Ntar,1);
% end
% if Ntar - numel(InOutFlag) > 0
% %    fprintf('Stress2 Outside\n')
%   stressT = op.nearSingInt(vesicle,f,S2diagOutFn,...
%       NearV2T,kernel2,kernel2,stressTar,false,false);
%   stress2(InOutFlag) = stressT(InOutFlag);
%   stress2(Ntar+InOutFlag) = stressT(Ntar+InOutFlag);
%   % correct stress at exterior points
% end
% Use near-singular integration to compute second component of
% the stress due to the single-layer potential

if ~isempty(RS)
  for k = 1:vesicle.nv
    cx = vesicle.center(1,k);
    cy = vesicle.center(2,k);
    % center of interior solid wall k
    xi = RS(3,k);
    % rotlet on interior solid wall k
    xi1 = RS(1,k);
    xi2 = RS(2,k);
    % first and second components of the stokeslet on 
    % interior solid wall k 
    [x,y] = oc.getXY(stressTar.X);
    rx = x - cx;
    ry = y - cy;
    rho4 = (rx.^2 + ry.^2).^2;
    stress1(1:Ntar) = stress1(1:Ntar) - 4*xi*rx.*ry./rho4;
    % add in (1,1) component of stress due to Rotlet
    stress1(1+Ntar:2*Ntar) = stress1(1+Ntar:2*Ntar) + ...
        2*xi*(rx.^2 - ry.^2)./rho4;
    % add in (1,2) component of stress due to Rotlet
    stress2(1:Ntar) = stress2(1:Ntar) + ...
        2*xi*(rx.^2 - ry.^2)./rho4;
    % add in (2,1) component of stress due to Rotlet
    stress2(1+Ntar:2*Ntar) = stress2(1+Ntar:2*Ntar) + ...
        4*xi*rx.*ry./rho4;
    % add in (2,2) component of stress due to Rotlet


    stress1(1:Ntar) = stress1(1:Ntar) + ...
        2*(rx*xi1 + ry*xi2).*(-rx.^2 + ry.^2)./rho4;
    stress1(1:Ntar) = stress1(1:Ntar) - ...
        2*(rx*xi1 + ry*xi2)./(rx.^2+ry.^2);
    % add in (1,1) component of stress due to Stokeslet
    % need to add in -1*pressure term as well
    stress1(1+Ntar:2*Ntar) = stress1(1+Ntar:2*Ntar) - ...
        4*rx.*ry.*(rx*xi1 + ry*xi2)./rho4;
    % add in (1,2) component of stress due to Stokeslet
    stress2(1:Ntar) = stress2(1:Ntar) - ...
        4*rx.*ry.*(rx*xi1 + ry*xi2)./rho4;
    % add in (1,2) component of stress due to Stokeslet
    stress2(1+Ntar:2*Ntar) = stress1(1+Ntar:2*Ntar) + ...
        2*(rx*xi1 + ry*xi2).*(rx.^2 - ry.^2)./rho4;
    stress2(1+Ntar:2*Ntar) = stress2(1+Ntar:2*Ntar) - ...
        2*(rx*xi1 + ry*xi2)./(rx.^2+ry.^2);
    % add in (2,2) component of stress due to Stokeslet
    % need to add in -1*pressure term as well
  end
end
% Add in contributions due to Rotlets and Stokeslets

end % stressTensor

end % methods

end %capsules



