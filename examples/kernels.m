classdef kernels 
% this class defines single and double layers for various kernels 
% (stokes, laplace) on 2D periodic curves.  Also defines the 
% integrals required for pressure and stress.
% Defines the matricies that map a density function defined on the
% boundary of a curve to the layer potential evaluated on the 
% curve, and defines the operators that takes a density function 
% defined on the boundary of a curve and returns the layer 
% potential at arbitrary target points.
% This class also has the main routine that evaluates layer
% potentials using near-singular integration.
    
properties
  N; % points per curve
  qw; % quadrature weights for logarithmic singularity
  qp; % quadrature points for logarithmic singularity (Alpert's rule)
  interpMat;  
  % interpolation matrix used for near-singular integration
  % This matrix replaces the need to use polyfit
  fmmPrecision;
  % precision of the fmm
  % fmmPrecision   Accuracy
  %  -2             5e-1
  %  -1             5e-2
  %   0             5e-3
  %   1             5e-4
  %   2             5e-7
  %   3             5e-10
  %   4             5e-13
  %   5             5e-16
  Nup;
  qwUp;
  qpUp;
  % upsampled quadrature rules for Alpert's quadrature rule.
  Rfor;
  Rbac;
  RforUp;
  RbacUp;
  % indicies that are used to do a rotation of columns of matricies for
  % implementing Alpert's quadrate rule
  gamma
  % indicies for Kapur-Rokhlin quadrature
  LPuprate
  % upsampling rate for layer potentials
  Prolong_LP
  Restrict_LP
  % prolongation and restriction matrices for layer potentitals

end % properties

methods 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = kernels(N)
% o = poten(N,Nup): constructor; N is the number of points per
% curve; Nup is the number of points on an upsampled grid that is used to remove antialiasing.  initialize class.

o.N = N;
o.interpMat = o.lagrangeInterp;
% load in the interpolation matrix which is precomputed
% with 7 interpolation points
o.fmmPrecision = 5;

accuracyOrder = 8;
o.qw = o.quadratureS(accuracyOrder);
o.qp = o.qw(:,2:end);
o.qw = o.qw(:,1);
[o.Rfor,o.Rbac] = o.rotationIndicies;

if o.N < 32 
  o.gamma = [ +1.825748064736159e0;...
            -1.325748064736161e0];
  % second-order accurate
else
  o.gamma = [+4.967362978287758e+0; ...
           -1.620501504859126e+1; ...
           +2.585153761832639e+1; ...
           -2.222599466791883e+1; ...
           +9.930104998037539e+0; ...
           -1.817995878141594e+0];
  % sixth-order accurate
end
% Load weights that are required near the log singularity

% For now do not antialias
o.LPuprate = 1;
o.Restrict_LP = [];
o.Prolong_LP = [];


if o.LPuprate~=1
  % need to build quadrature weights for upsampled grid  
  o.N = o.N*o.LPuprate;
  o.qwUp = o.quadratureS(accuracyOrder);
  o.qpUp = o.qwUp(:,2:end);
  % matrix that takes values on equispaced grid and assignvs values at
  % Alpert's quadrature nodes
  o.qwUp = o.qwUp(:,1);
  % Alpret's quadrature weights
  [o.RforUp,o.RbacUp] = o.rotationIndicies;
  % build rotation indices for upsampled grid
  o.Nup = o.N;
  o.N = N;
  % move back to the number of points that we have on each vesicle
end


end % kernels: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = stokesSLmatrixAlpert(o,vesicle,viscosity)
% G = stokesSLmatrix(vesicle) generates the single-layer potential for
% Stokes vesicle is a data structure defined as in the curve class G is
% (2N,2N,nv) array where N is the number of points per curve and nv is
% the number of curves in X 

% It can be upsampled or not, if vesicle.N > o.N then upsampled vesicle is
% sent here.

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions
constTerm = 1/(viscosity);

if vesicle.N > o.N
  Nquad = numel(o.qwUp);
  qw = o.qwUp(:,ones(vesicle.N,1));
  qp = o.qpUp;
  Rbac = o.RbacUp;
  Rfor = o.RforUp;
else
  Nquad = numel(o.qw);
  qw = o.qw(:,ones(vesicle.N,1));
  qp = o.qp;
  Rbac = o.Rbac;
  Rfor = o.Rfor;
end
% number of quadrature points including the extra terms that Alpert's
% rule brings into the quadrature


G = zeros(2*o.N,2*o.N,vesicle.nv);
for k=1:vesicle.nv  % Loop over curves
  xx = x(:,k);
  yy = y(:,k);
  % locations
  sa = vesicle.sa(:,k)';
  sa = sa(ones(vesicle.N,1),:);
  % Jacobian

  xtar = xx(:,ones(Nquad,1))'; 
  ytar = yy(:,ones(Nquad,1))'; 
  % target points

  xsou = xx(:,ones(vesicle.N,1)); 
  ysou = yy(:,ones(vesicle.N,1));
  % source points
  
  xsou = xsou(Rfor);
  ysou = ysou(Rfor);
  % have to rotate each column so that it is compatiable with o.qp
  % which is the matrix that takes function values and maps them to the
  % intermediate values required for Alpert quadrature
  
  diffx = xtar - qp*xsou;
  diffy = ytar - qp*ysou;
  rho2 = (diffx.^2 + diffy.^2).^(-1);
  % one over distance squared

  logpart = 0.5*qp'*(qw .* log(rho2));
  % sign changed to positive because rho2 is one over distance squared

  Gves = logpart + qp'*(qw.*diffx.^2.*rho2);
  Gves = Gves(Rbac);
  G11 = Gves'.*sa;
  % (1,1)-term
 
  Gves = logpart + qp'*(qw.*diffy.^2.*rho2);
  Gves = Gves(Rbac);
  G22 = Gves'.*sa;
  % (2,2)-term

  Gves = qp'*(qw.*diffx.*diffy.*rho2);
  Gves = Gves(Rbac);
  G12 = Gves'.*sa;
  % (1,2)-term
  
  if vesicle.N > o.N
    G(1:o.N,1:o.N,k) = o.Restrict_LP * G11 * o.Prolong_LP;
    G(1:o.N,o.N+1:end,k) = o.Restrict_LP * G12 * o.Prolong_LP;
    G(o.N+1:end,1:o.N,k) = G(1:o.N,o.N+1:end,k);
    G(o.N+1:end,o.N+1:end,k) = o.Restrict_LP * G22 * o.Prolong_LP;
  else
    G(1:o.N,1:o.N,k) = G11;
    G(1:o.N,o.N+1:end,k) = G12;
    G(o.N+1:end,1:o.N,k) = G(1:o.N,o.N+1:end,k);
    G(o.N+1:end,o.N+1:end,k) = G22;
  end
end
G = constTerm * G;
end % stokesSLmatrixAlpert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLmatrix(o,vesicle)
% D = stokesDLmatrix(vesicle), generate double-layer potential for 
% Stokes vesicle is a data structure defined as in the capsules class
% D is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions
N = vesicle.N;
% number of points per vesicles
D = zeros(2*N,2*N,vesicle.nv);
% initialize space for double-layer potential matrix

if vesicle.N > o.N
  N = o.LPuprate*N;
end

for k=1:vesicle.nv  % Loop over curves
  if (vesicle.viscCont(k) ~= 1)
    const_coeff = -(1-vesicle.viscCont(k));
    % constant that has the d\theta and scaling with the viscosity
    % contrast
    xx = x(:,k);
    yy = y(:,k);
    % locations

    if vesicle.N > o.N
      xx = interpft(xx,N); yy = interpft(yy,N);
      vesicleUp = capsules([xx;yy],[],[],...
          vesicle.kappa,vesicle.viscCont(k),true);
      % upsampled single versicle
      [tx,ty] = oc.getXY(vesicleUp.xt);
      % Vesicle tangent
      sa = vesicleUp.sa';
      % Jacobian
      cur = vesicleUp.cur';
    else
      [tx,ty] = oc.getXY(vesicle.xt);
      tx = tx(:,k); ty = ty(:,k);
      % Vesicle tangent
      sa = vesicle.sa(:,k)';
      % Jacobian
      cur = vesicle.cur(:,k)';
      % curvature
    end

    xtar = xx(:,ones(N,1))';
    ytar = yy(:,ones(N,1))';
    % target points

    xsou = xx(:,ones(N,1));
    ysou = yy(:,ones(N,1));
    % source points

    txsou = tx';
    tysou = ty';
    % tangent at srouces
    sa = sa(ones(N,1),:);
    % Jacobian

    diffx = xtar - xsou;
    diffy = ytar - ysou;
    rho4 = (diffx.^2 + diffy.^2).^(-2);
    rho4(1:N+1:N.^2) = 0;
    % set diagonal terms to 0

    kernel = diffx.*(tysou(ones(N,1),:)) - ...
            diffy.*(txsou(ones(N,1),:));
    kernel = kernel.*rho4.*sa;
    kernel = const_coeff*kernel;

    D11 = kernel.*diffx.^2;
    % (1,1) component
    D11(1:N+1:N.^2) = 0.5*const_coeff*cur.*sa(1,:).*txsou.^2;
    % diagonal limiting term

    D12 = kernel.*diffx.*diffy;
    % (1,2) component
    D12(1:N+1:N.^2) = 0.5*const_coeff*cur.*sa(1,:).*txsou.*tysou;
    % diagonal limiting term

    D22 = kernel.*diffy.^2;
    % (2,2) component
    D22(1:N+1:N.^2) = 0.5*const_coeff*cur.*sa(1,:).*tysou.^2;
    % diagonal limiting term

    if vesicle.N > o.N
      D11 = o.Restrict_LP * D11 * o.Prolong_LP; 
      D12 = o.Restrict_LP * D12 * o.Prolong_LP; 
      D22 = o.Restrict_LP * D22 * o.Prolong_LP; 
    end
    % move to grid with vesicle.N points by applying prolongation and
    % restriction operators

    D(:,:,k) = [D11 D12; D12 D22];
    % build matrix with four blocks
    D(:,:,k) = D(:,:,k)*2/N;
    % scale with the arclength spacing and divide by pi

  end
end % k

end % stokesDLmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLTmatrix(o,vesicle)
% D = stokesDLTmatrix(vesicle), generate adjoint of double-layer potential for 
% Stokes vesicle is a data structure defined as in the capsules class
% D is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions
N = vesicle.N;
% number of points per vesicles
D = zeros(2*N,2*N,vesicle.nv);
% initialize space for double-layer potential matrix

for k=1:vesicle.nv  % Loop over curves
  if (vesicle.viscCont(k) ~= 1)
    const_coeff = -1;
    % constant that has the d\theta and scaling with the viscosity
    % contrast
    xx = x(:,k);
    yy = y(:,k);
    % locations

    if vesicle.N > o.N
      xx = interpft(xx,N); yy = interpft(yy,N);
      vesicleUp = capsules([xx;yy],[],[],...
          vesicle.kappa,vesicle.viscCont(k),true);
      % upsampled single versicle
      [tx,ty] = oc.getXY(vesicleUp.xt);
      % Vesicle tangent
      sa = vesicleUp.sa';
      % Jacobian
      cur = vesicleUp.cur;
      normal = vesicleUp.normal;
    else
      [tx,ty] = oc.getXY(vesicle.xt);
      tx = tx(:,k); ty = ty(:,k);
      % Vesicle tangent
      sa = vesicle.sa(:,k)';
      % Jacobian
      cur = vesicle.cur(:,k);
      % curvature
      normal = vesicle.normal(:,k);
    end

    xtar = xx(:,ones(N,1))';
    ytar = yy(:,ones(N,1))';
    % target points
    
    normx = normal(1:N,ones(N,1));
    normy = normal(N+1:2*N,ones(N,1));

    % normal points
    xsou = xx(:,ones(N,1));
    ysou = yy(:,ones(N,1));
    % source points

    txtar = tx;
    tytar = ty;
    % tangent at srouces
    sa = sa(ones(N,1),:);
    % Jacobian

    diffx = xtar - xsou;
    diffy = ytar - ysou;
    rho4 = (diffx.^2 + diffy.^2).^(-2);
    rho4(1:N+1:N.^2) = 0;
    % set diagonal terms to 0
    
    kernel = diffx.*normx + diffy.*normy;
    kernel = kernel.*rho4.*sa;
    kernel = const_coeff*kernel;

    D11 = kernel.*diffx.^2;
    % (1,1) component

    D11(1:N+1:N.^2) = 0.5*const_coeff*cur'.*sa(1,:).*txtar'.^2;
    % diagonal limiting term

    D12 = kernel.*diffx.*diffy;
    % (1,2) component
    D12(1:N+1:N.^2) = 0.5*const_coeff*cur'.*sa(1,:).*txtar'.*tytar';
    % diagonal limiting term

    D22 = kernel.*diffy.^2;
    % (2,2) component
    D22(1:N+1:N.^2) = 0.5*const_coeff*cur'.*sa(1,:).*tytar'.^2;
    % diagonal limiting term

    if vesicle.N > o.N
      D11 = o.Restrict_LP * D11 * o.Prolong_LP; 
      D12 = o.Restrict_LP * D12 * o.Prolong_LP; 
      D22 = o.Restrict_LP * D22 * o.Prolong_LP; 
    end
    % move to grid with vesicle.N points by applying prolongation and
    % restriction operators

    D(:,:,k) = [D11 D12; D12 D22];
    % build matrix with four blocks
    D(:,:,k) = -D(:,:,k)*2/N;
    % scale with the arclength spacing and divide by pi

  end
end % k

end % stokesDLTmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLmatrix2(o,vesicle)
% D = stokesDLmatrix(vesicle), generate double-layer potential for 
% Stokes vesicle is a data structure defined as in the capsules class
% D is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions
N = vesicle.N;
% number of points per vesicles
D = zeros(2*N,2*N,vesicle.nv);
% initialize space for double-layer potential matrix


for k=1:vesicle.nv  % Loop over curves
  if (vesicle.viscCont(k) ~= 1)
    
    % constant that has the d\theta and scaling with the viscosity
    % contrast
    xx = x(:,k);
    yy = y(:,k);
    % locations

    
    [tx,ty] = oc.getXY(vesicle.xt);
    tx = tx(:,k); ty = ty(:,k);
    % Vesicle tangent
    sa = vesicle.sa(:,k)';
    % Jacobian
    cur = vesicle.cur(:,k)';
    % curvature
   
    normx = vesicle.normal(1:N,k)';
    normy = vesicle.normal(N+1:2*N,k)';

    xtar = xx(:,ones(N,1))';
    ytar = yy(:,ones(N,1))';
    % target points

    xsou = xx(:,ones(N,1));
    ysou = yy(:,ones(N,1));
    % source points

    txsou = tx';
    tysou = ty';
    % tangent at srouces
    sa = sa(ones(N,1),:);
    % Jacobian

    diffx = xtar - xsou;
    diffy = ytar - ysou;
    rho4 = (diffx.^2 + diffy.^2).^(-2);
    rho4(1:N+1:N.^2) = 0;
    % set diagonal terms to 0

    kernel = diffx.*normx(ones(N,1),:) + diffy.*(normy(ones(N,1),:));
    kernel = kernel.*rho4.*sa;

    D11 = 1/pi*kernel.*diffx.^2;
    % (1,1) component
    D11(1:N+1:N.^2) = -0.5/pi*cur.*sa(1,:).*txsou.^2;
    % diagonal limiting term

    D12 = 1/pi*kernel.*diffx.*diffy;
    % (1,2) component
    D12(1:N+1:N.^2) = -0.5/pi*cur.*sa(1,:).*txsou.*tysou;
    % diagonal limiting term

    D22 = 1/pi*kernel.*diffy.^2;
    % (2,2) component
    D22(1:N+1:N.^2) = -0.5/pi*cur.*sa(1,:).*tysou.^2;
    % diagonal limiting term

    D(:,:,k) = [D11 D12; D12 D22];
    % build matrix with four blocks
    D(:,:,k) = D(:,:,k)*2*pi/N;
    % scale with the arclength spacing and divide by pi

  end
end % k

end % stokesDLmatrix2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLTmatrix2(o,vesicle)
% D = stokesDLTmatrix(vesicle), generate adjoint of double-layer potential for 
% Stokes vesicle is a data structure defined as in the capsules class
% D is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions
N = vesicle.N;
% number of points per vesicles
D = zeros(2*N,2*N,vesicle.nv);
% initialize space for double-layer potential matrix

for k=1:vesicle.nv  % Loop over curves
  if (vesicle.viscCont(k) ~= 1)
 
    % constant that has the d\theta and scaling with the viscosity
    % contrast
    xx = x(:,k);
    yy = y(:,k);
    % locations

    
    [tx,ty] = oc.getXY(vesicle.xt);
    tx = tx(:,k); ty = ty(:,k);
    % Vesicle tangent
    sa = vesicle.sa(:,k)';
    % Jacobian
    cur = vesicle.cur(:,k);
    % curvature
    normal = vesicle.normal(:,k);
    

    xtar = xx(:,ones(N,1))';
    ytar = yy(:,ones(N,1))';
    % target points
    
    normx = normal(1:N,ones(N,1));
    normy = normal(N+1:2*N,ones(N,1));

    % normal points
    xsou = xx(:,ones(N,1));
    ysou = yy(:,ones(N,1));
    % source points

    txtar = tx;
    tytar = ty;
    % tangent at srouces
    sa = sa(ones(N,1),:);
    % Jacobian

    diffx = xtar - xsou;
    diffy = ytar - ysou;
    rho4 = (diffx.^2 + diffy.^2).^(-2);
    rho4(1:N+1:N.^2) = 0;
    % set diagonal terms to 0
    
    kernel = diffx.*normx + diffy.*normy;
    kernel = kernel.*rho4.*sa;

    D11 = 1/pi*kernel.*diffx.^2;
    % (1,1) component

    D11(1:N+1:N.^2) = 0.5/pi*cur'.*sa(1,:).*txtar'.^2;
    % diagonal limiting term

    D12 = 1/pi*kernel.*diffx.*diffy;
    % (1,2) component
    D12(1:N+1:N.^2) = 0.5/pi*cur'.*sa(1,:).*txtar'.*tytar';
    % diagonal limiting term

    D22 = 1/pi*kernel.*diffy.^2;
    % (2,2) component
    D22(1:N+1:N.^2) = 0.5/pi*cur'.*sa(1,:).*tytar'.^2;
    % diagonal limiting term

    D(:,:,k) = [D11 D12; D12 D22];
    % build matrix with four blocks
    D(:,:,k) = -1 * D(:,:,k)*2*pi/N;
    % scale with the arclength spacing and divide by pi

  end
end % k

end % stokesDLTmatrix2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qw = quadratureS(o,q)
% qw = quadratureS(q) generates the quadrature rules for a function
% with o.N points and a logarithmic singularity at the origin.  q
% controls the accuracy.  All rules are from Alpert 1999.  This is
% Bryan's reformulation which uses Alpert's quadrature rules as
% described in Section 7, but with Fourier interpolation

[v,u,a] = o.getWeights(q);
% get the weights coming from Table 8 of Alpert's 1999 paper

h = 2*pi/o.N;
n = o.N - 2*a + 1;

of = fft1;
A1 = of.sinterpS(o.N,v*h);
A2 = of.sinterpS(o.N,2*pi-flipud(v*h));
yt = h*(a:n-1+a)';
% regular points away from the singularity
wt = [h*u; h*ones(length(yt),1); h*flipud(u)]/4/pi;
% quadrature points away from singularity

B = sparse(length(yt),o.N);
pos = 1 + (a:n-1+a)';

for k = 1:length(yt)
  B(k, pos(k)) = 1;
end
A = [sparse(A1); B; sparse(A2)];
qw = [wt, A];

end % quadratureS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,u,a] = getWeights(o,q)
% [v,u,a] = getWeights(q) loads quadrature rules for different types of
% singularities.  All rules come from Bradley Alpert's papers.  We are
% interested in nodesLogx.dat

xp = o.nodesLog;
lenth = [3;7;15];
par = [2;5;10];


switch q
case 4
  v = xp(1:lenth(1), 1);
  u = xp(1:lenth(1), 2);
  a = par(1);
case 8
  v = xp(1+lenth(1):lenth(1)+lenth(2),1);
  u = xp(1+lenth(1):lenth(1)+lenth(2),2);
  a = par(2);
case 16
  v = xp(1+lenth(2)+lenth(1):sum(lenth),1);
  u = xp(1+lenth(2)+lenth(1):sum(lenth),2);
  a = par(3);
end

end % getWeights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rfor,Rbac] = rotationIndicies(o)
% [Rfor,Rbac] = rotationIndicies generates two matricies of indicies
% that are used in stokesSLmatrix to implement Alpert's quadrature
% rules

ind = (1:o.N)';
Rfor = zeros(o.N);
Rbac = zeros(o.N);
% vector of indicies so that we can apply circshift to each column
% efficiently.  Need one for going 'forward' and one for going
% 'backwards'
Rfor(:,1) = ind;
Rbac(:,1) = ind;
for k = 2:o.N
  Rfor(:,k) = (k-1)*o.N + [ind(k:o.N);ind(1:k-1)];
  Rbac(:,k) = (k-1)*o.N + [...
      ind(o.N-k+2:o.N);ind(1:o.N-k+1)];
end

end % rotationIndicies

end % methods


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)

function LP = lagrangeInterp
% interpMap = lagrangeInterp builds the Lagrange interpolation
% matrix that takes seven function values equally distributed
% in [0,1] and returns the seven polynomial coefficients

interpMat = zeros(7);
LP(1,1) = 6.48e1;
LP(1,2) = -3.888e2;
LP(1,3) = 9.72e2;
LP(1,4) = -1.296e3;
LP(1,5) = 9.72e2;
LP(1,6) = -3.888e2;
LP(1,7) = 6.48e1;

LP(2,1) = -2.268e2;
LP(2,2) = 1.296e3;
LP(2,3) = -3.078e3;
LP(2,4) = 3.888e3;
LP(2,5) = -2.754e3;
LP(2,6) = 1.0368e3;
LP(2,7) = -1.62e2;

LP(3,1) = 3.15e2;
LP(3,2) = -1.674e3;
LP(3,3) = 3.699e3;
LP(3,4) = -4.356e3;
LP(3,5) = 2.889e3;
LP(3,6) = -1.026e3;
LP(3,7) = 1.53e2;

LP(4,1) = -2.205e2;
LP(4,2) = 1.044e3;
LP(4,3) = -2.0745e3;
LP(4,4) = 2.232e3;
LP(4,5) = -1.3815e3;
LP(4,6) = 4.68e2;
LP(4,7) = -6.75e1;

LP(5,1) = 8.12e1;
LP(5,2) = -3.132e2;
LP(5,3) = 5.265e2;
LP(5,4) = -5.08e2;
LP(5,5) = 2.97e2;
LP(5,6) = -9.72e1;
LP(5,7) = 1.37e1;

LP(6,1) = -1.47e1;
LP(6,2) = 3.6e1;
LP(6,3) = -4.5e1;
LP(6,4) = 4.0e1;
LP(6,5) = -2.25e1;
LP(6,6) = 7.2e0;
LP(6,7) = -1e0;

LP(7,1) = 1e0;
% rest of the coefficients are zero

end % lagrangeInterp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xp = nodesLog

xp = zeros(25,2);
xp(1,1) = 2.379647284118974e-2;
xp(1,2) = 8.795942675593887e-2;
xp(2,1) = 2.935370741501914e-1;
xp(2,2) = 4.989017152913699e-1;
xp(3,1) = 1.023715124251890e+0;
xp(3,2) = 9.131388579526912e-1;
xp(4,1) = 6.531815708567918e-3;
xp(4,2) = 2.462194198995203e-2;
xp(5,1) = 9.086744584657729e-2;
xp(5,2) = 1.701315866854178e-1;
xp(6,1) = 3.967966533375878e-1;
xp(6,2) = 4.609256358650077e-1;
xp(7,1) = 1.027856640525646e+0;
xp(7,2) = 7.947291148621895e-1;
xp(8,1) = 1.945288592909266e+0;
xp(8,2) = 1.008710414337933e+0;
xp(9,1) = 2.980147933889640e+0;
xp(9,2) = 1.036093649726216e+0;
xp(10,1) = 3.998861349951123e+0;
xp(10,2) = 1.004787656533285e+0;
xp(11,1) = 8.371529832014113e-4;
xp(11,2) = 3.190919086626234e-3;
xp(12,1) = 1.239382725542637e-2;
xp(12,2) = 2.423621380426338e-2;
xp(13,1) = 6.009290785739468e-2;
xp(13,2) = 7.740135521653088e-2;
xp(14,1) = 1.805991249601928e-1;
xp(14,2) = 1.704889420286369e-1;
xp(15,1) = 4.142832599028031e-1;
xp(15,2) = 3.029123478511309e-1;
xp(16,1) = 7.964747731112430e-1;
xp(16,2) = 4.652220834914617e-1;
xp(17,1) = 1.348993882467059e+0;
xp(17,2) = 6.401489637096768e-1;
xp(18,1) = 2.073471660264395e+0;
xp(18,2) = 8.051212946181061e-1;
xp(19,1) = 2.947904939031494e+0;
xp(19,2) = 9.362411945698647e-1;
xp(20,1) = 3.928129252248612e+0;
xp(20,2) = 1.014359775369075e+0;
xp(21,1) = 4.957203086563112e+0;
xp(21,2) = 1.035167721053657e+0;
xp(22,1) = 5.986360113977494e+0;
xp(22,2) = 1.020308624984610e+0;
xp(23,1) = 6.997957704791519e+0;
xp(23,2) = 1.004798397441514e+0;
xp(24,1) = 7.999888757524622e+0;
xp(24,2) = 1.000395017352309e+0;
xp(25,1) = 8.999998754306120e+0;
xp(25,2) = 1.000007149422537e+0;

end % nodesLog

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xp = nodesRegular

xp = zeros(35,2);
xp(1,1) = 2.000000000000000e-1;
xp(1,2) = 5.208333333333333e-1;
xp(2,1) = 1.000000000000000e+0;
xp(2,2) = 9.791666666666667e-1;
xp(3,1) = 2.250991042610971e-1;
xp(3,2) = 5.549724327164181e-1;
xp(4,1) = 1.014269060987992e+0;
xp(4,2) = 9.451317411845474e-1;
xp(5,1) = 2.000000000000000e+0;
xp(5,2) = 9.998958260990347e-1;
xp(6,1) = 2.087647422032129e-1;
xp(6,2) = 5.207988277246498e-1;
xp(7,1) = 9.786087373714483e-1;
xp(7,2) = 9.535038018555888e-1;
xp(8,1) = 1.989541386579751e+0;
xp(8,2) = 1.024871626402471e+0;
xp(9,1) = 3.000000000000000e+0;
xp(9,2) = 1.000825744017291e+0;
xp(10,1) = 7.023955461621939e-2;
xp(10,2) = 1.922315977843698e-1;
xp(11,1) = 4.312297857227970e-1;
xp(11,2) = 5.348399530514687e-1;
xp(12,1) = 1.117752734518115e+0;
xp(12,2) = 8.170209442488760e-1;
xp(13,1) = 2.017343724572518e+0;
xp(13,2) = 9.592111521445966e-1;
xp(14,1) = 3.000837842847590e+0;
xp(14,2) = 9.967143408044999e-1;
xp(15,1) = 4.000000000000000e+0;
xp(15,2) = 9.999820119661890e-1;
xp(16,1) = 9.919337841451029e-2;
xp(16,2) = 2.528198928766921e-1;
xp(17,1) = 5.076592669645529e-1;
xp(17,2) = 5.550158230159487e-1;
xp(18,1) = 1.184972925827278e+0;
xp(18,2) = 7.852321453615224e-1;
xp(19,1) = 2.047493467134072e+0;
xp(19,2) = 9.245915673876715e-1;
xp(20,1) = 3.007168911869310e+0;
xp(20,2) = 9.839350200445296e-1;
xp(21,1) = 4.000474996776184e+0;
xp(21,2) = 9.984463448413151e-1;
xp(22,1) = 5.000007879022339e+0;
xp(22,2) = 9.999592378464547e-1;
xp(23,1) = 6.000000000000000e+0;
xp(23,2) = 9.999999686258662e-1;
xp(24,1) = 6.001064731474805e-2;
xp(24,2) = 1.538932104518340e-1;
xp(25,1) = 3.149685016229433e-1;
xp(25,2) = 3.551058128559424e-1;
xp(26,1) = 7.664508240518316e-1;
xp(26,2) = 5.449200036280008e-1;
xp(27,1) = 1.396685781342510e+0;
xp(27,2) = 7.104078497715549e-1;
xp(28,1) = 2.175195903206602e+0;
xp(28,2) = 8.398780940253654e-1;
xp(29,1) = 3.062320575880355e+0;
xp(29,2) = 9.272767950890611e-1;
xp(30,1) = 4.016440988792476e+0;
xp(30,2) = 9.750605697371132e-1;
xp(31,1) = 5.002872064275734e+0;
xp(31,2) = 9.942629650823470e-1;
xp(32,1) = 6.000285453310164e+0;
xp(32,2) = 9.992421778421898e-1;
xp(33,1) = 7.000012964962529e+0;
xp(33,2) = 9.999534370786161e-1;
xp(34,1) = 8.000000175554469e+0;
xp(34,2) = 9.999990854912925e-1;
xp(35,1) = 9.000000000000000e+0;
xp(35,2) = 9.999999989466828e-1;

end % nodesRegular

end % methods (Static)

end % classdef kernels