classdef curve
% This class implements that basic calculus on the curve.
% The basic data structure is a matrix X in which the columns 
% represent periodic C^{\infty} closed curves with N points, 
% X(1:n,j) is the x-coordinate of the j_th curve and X(n+1:N,j) 
% is the y-coordinate of the j_th curve; here n=N/2
% X coordinates do not have to be periodic, but the curvature,
% normals, etc that they compute will be garbage.  This allows
% us to store target points for tracers or the pressure using
% this class and then using near-singular integration is easy
% to implement

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=getXY(o,X)
% [x,y] = getXY(X) get the [x,y] component of curves X
N = size(X,1)/2;
x = X(1:N,:);
y = X(N+1:end,:);

end % getXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = setXY(o,x,y)
% V = setXY(x,y) set the [x,y] component of vector V on the curve
N = size(x,1);
V=zeros(2*N,size(x,2));
V(1:N,:) = x;
V(N+1:end,:) = y;

end % setXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function center = getCenter(o,X)
% center = getCenter(o,X) finds the center of each capsule
N = size(X,1)/2;
nv = size(X,2);

% upsample
Nup = max(N,256);
x = interpft(X(1:end/2,:),Nup);
y = interpft(X(end/2+1:end,:),Nup);

center = zeros(nv,1);

for k = 1 : nv
  center(k) = sqrt(mean(x)^2 + mean(y)^2);
end

end % getCenter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IA = getIncAngle(o,X)
% IA = getIncAngle(o,X) finds the inclination angle of each capsule
% The inclination angle (IA) is the angle between the x-axis and the 
% principal axis corresponding to the smallest principal moment of inertia
nv = size(X,2);
IA = zeros(nv,1);

% compute inclination angle on an upsampled grid
N = max(256,size(X,1)/2);
modes = [(0:N/2-1)';0;(-N/2+1:-1)'];


for k = 1 : nv
  X(:,k) = [X(1:end/2,k)-mean(interpft(X(1:end/2,k),N));...
      X(end/2+1:end,k)-mean(interpft(X(end/2+1:end,k),N))];
end

for k = 1 : nv
    x = interpft(X(1:end/2,k),N); 
    y = interpft(X(end/2+1:end,k),N);
    
    Dx = real(ifft(1i*modes.*fft(x)));
    Dy = real(ifft(1i*modes.*fft(y)));
    jac = sqrt(Dx.^2 + Dy.^2);
    tx = Dx./jac; ty = Dy./jac;
    nx = ty; ny = -tx;
    rdotn = x.*nx + y.*ny;
    rho2 = x.^2 + y.^2;

    J11 = 0.25*sum(rdotn.*(rho2 - x.*x).*jac)*2*pi/N;
    J12 = 0.25*sum(rdotn.*(-x.*y).*jac)*2*pi/N;
    J21 = 0.25*sum(rdotn.*(-y.*x).*jac)*2*pi/N;
    J22 = 0.25*sum(rdotn.*(rho2 - y.*y).*jac)*2*pi/N;

    J = [J11 J12; J21 J22];
    [V,D] = eig(J);
    
    [~,ind] = min(abs(diag(D)));
    % make sure that the first components of e-vectors have the same sign
    if V(2,ind)<0
      V(:,ind) = -1*V(:,ind);
    end
    % IA is in [-pi, pi]
    IA(k) = atan2(V(2,ind),V(1,ind));
    % move IA to [0,2*pi];
    if IA(k) < 0
      IA(k) = IA(k)+2*pi;
    end
end

end % getIncAngle



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = getPrincAxes(o,X)
% IA = getIncAngle(o,X) finds the inclination angle of each capsule
% The inclination angle (IA) is the angle between the x-axis and the 
% principal axis corresponding to the smallest principal moment of inertia

% compute inclination angle on an upsampled grid
N = max(256,size(X,1)/2);
modes = [(0:N/2-1)';0;(-N/2+1:-1)'];



X = [X(1:end/2)-mean(interpft(X(1:end/2),N));...
    X(end/2+1:end)-mean(interpft(X(end/2+1:end),N))];



x = interpft(X(1:end/2),N); 
y = interpft(X(end/2+1:end),N);

Dx = real(ifft(1i*modes.*fft(x)));
Dy = real(ifft(1i*modes.*fft(y)));
jac = sqrt(Dx.^2 + Dy.^2);
tx = Dx./jac; ty = Dy./jac;
nx = ty; ny = -tx;
rdotn = x.*nx + y.*ny;
rho2 = x.^2 + y.^2;

J11 = 0.25*sum(rdotn.*(rho2 - x.*x).*jac)*2*pi/N;
J12 = 0.25*sum(rdotn.*(-x.*y).*jac)*2*pi/N;
J21 = 0.25*sum(rdotn.*(-y.*x).*jac)*2*pi/N;
J22 = 0.25*sum(rdotn.*(rho2 - y.*y).*jac)*2*pi/N;

J = [J11 J12; J21 J22];
[V,D] = eig(J);

[~,ind] = min(abs(diag(D)));
V = V(:,ind);

end % getPrincAxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IA = getIncAngle2(o,X)
% IA = getIncAngle(o,X) finds the inclination angle of each capsule
% The inclination angle (IA) is the angle between the x-axis and the 
% principal axis corresponding to the smallest principal moment of inertia
nv = size(X,2);
IA = zeros(nv,1);

% compute inclination angle on an upsampled grid
N = max(256,size(X,1)/2);
modes = [(0:N/2-1)';0;(-N/2+1:-1)'];

areaDependentAngle = true;

for k = 1 : nv
  X(:,k) = [X(1:end/2,k)-mean(interpft(X(1:end/2,k),N));...
      X(end/2+1:end,k)-mean(interpft(X(end/2+1:end,k),N))];
end

for k = 1 : nv
    x = interpft(X(1:end/2,k),N); 
    y = interpft(X(end/2+1:end,k),N);
    
    Dx = real(ifft(1i*modes.*fft(x)));
    Dy = real(ifft(1i*modes.*fft(y)));
    jac = sqrt(Dx.^2 + Dy.^2);
    tx = Dx./jac; ty = Dy./jac;
    nx = ty; ny = -tx;
    rdotn = x.*nx + y.*ny;
    rho2 = x.^2 + y.^2;

    J11 = 0.25*sum(rdotn.*(rho2 - x.*x).*jac)*2*pi/N;
    J12 = 0.25*sum(rdotn.*(-x.*y).*jac)*2*pi/N;
    J21 = 0.25*sum(rdotn.*(-y.*x).*jac)*2*pi/N;
    J22 = 0.25*sum(rdotn.*(rho2 - y.*y).*jac)*2*pi/N;

    J = [J11 J12; J21 J22];
    [V,D] = eig(J);
    
    [~,ind] = min(abs(diag(D)));
    % make sure that the first components of e-vectors have the same sign
    if V(2,ind)<0
      V(:,ind) = -1*V(:,ind);
    end
    % since V(2,ind) > 0, this will give angle between [0, pi]
    IA(k) = atan2(V(2,ind),V(1,ind));
    
    if areaDependentAngle 
    % FIND DIRECTION OF PRIN. AXIS DEPENDING ON HEAD-TAIL
    % 1) already translated to (0,0), so rotate to pi/2
    x0rot = x*cos(-IA(k)+pi/2) - y*sin(-IA(k)+pi/2);
    y0rot = x*sin(-IA(k)+pi/2) + y*cos(-IA(k)+pi/2);
    
    % 2) find areas (top, bottom)
    % need derivatives, so rotate the computed ones
    Dx = Dx*cos(-IA(k)+pi/2) - Dy*sin(-IA(k)+pi/2);
    Dy = Dx*sin(-IA(k)+pi/2) + Dy*cos(-IA(k)+pi/2);
    % Compute again, that is also fast
    %Dx = real(ifft(1i*modes.*fft(x0rot)));
    %Dy = real(ifft(1i*modes.*fft(y0rot)));
    
    idcsTop = find(y0rot>=0); idcsBot = find(y0rot<0);
    areaTop = sum(x0rot(idcsTop).*Dy(idcsTop)-y0rot(idcsTop).*...
        Dx(idcsTop))/N*pi;
    
    areaBot = sum(x0rot(idcsBot).*Dy(idcsBot)-y0rot(idcsBot).*...
        Dx(idcsBot))/N*pi;
    % DEBUG 
    idebug = false;
    if idebug
    figure(2); clf;
    plot(x,y,'b')
    axis equal
    figure(3); clf; hold on;
    plot(x0rot,y0rot,'b-o')
    axis equal
    plot(x0rot(idcsTop),y0rot(idcsTop),'ro')
    plot(x0rot(idcsBot),y0rot(idcsBot),'go')
    legend('vesicle at IA = pi/2','top idcs','bot idcs'); legend boxoff;
    disp(['area of top section: ' num2str(areaTop)])
    disp(['area of bottom section: ' num2str(areaBot)])
    end % if idebug
    
    if areaBot >= 1.1*areaTop
      IA(k) = IA(k) + pi;
    elseif areaTop < 1.1*areaBot  
      % if areaTop ~ areaBot, then check areaRight, areaLeft  
      idcsLeft = find(x0rot<0); idcsRight = find(x0rot>=0);
      areaRight = sum(x0rot(idcsRight).*Dy(idcsRight)-y0rot(idcsRight).*...
        Dx(idcsRight))/N*pi;
      areaLeft = sum(x0rot(idcsLeft).*Dy(idcsLeft)-y0rot(idcsLeft).*...
        Dx(idcsLeft))/N*pi;
      if areaLeft >= 1.1*areaRight
        IA(k) = IA(k) + pi;
      end
      % debug
    end
    end % if areaDependentAngle

end
end % getIncAngle2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dx,Dy]=getDXY(o,X)
% [Dx,Dy]=getDXY(X), compute the derivatives of each component of X 
% these are the derivatives with respect the parameterization 
% not arclength
x = X(1:end/2,:);
y = X(end/2+1:end,:);
N = size(x,1);
nv = size(x,2);
IK = fft1.modes(N,nv);
Dx = fft1.diffFT(x,IK);
Dy = fft1.diffFT(y,IK);

end % getDXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jacobian,tangent,curvature] = diffProp(o,X)
% [jacobian,tangent,curvature] = diffProp(X) returns differential
% properties of the curve each column of the matrix X. Each column of 
% X should be a closed curve defined in plane. The tangent is the 
% normalized tangent.
%
% EXAMPLE:
%    n = 128; nv = 3;
%    X = boundary(n,'nv',nv,'curly');
%    c = curve;
%    [k t s] = c.diffProp(X);

N = size(X,1)/2;
nv = size(X,2);

% get the x y components
[Dx,Dy] = o.getDXY(X);

jacobian = sqrt(Dx.^2 + Dy.^2); 

if nargout>1  % if user requires tangent
  tangent = o.setXY( Dx./jacobian, Dy./jacobian);
end

if nargout>2  % if user requires curvature
  IK = fft1.modes(N,nv);
  DDx = curve.arcDeriv(Dx,1,ones(N,nv),IK);
  DDy = curve.arcDeriv(Dy,1,ones(N,nv),IK);
  curvature = (Dx.*DDy - Dy.*DDx)./(jacobian.^3);
end
% curvature of the curve

end % diffProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [reducedArea,area,length] = geomProp(o,X)
% [reducedArea area length] = geomProp(X) calculate the length, area 
% and the reduced volume of domains inclose by columns of X. 
% Reduced volume is defined as 4*pi*A/L. 
% EXAMPLE:
%   X = boundary(64,'nv',3,'curly');
%   c = curve(X);
%   [rv A L] = c.geomProp(X);

if size(X,1) < 256
  X = [interpft(X(1:end/2,:),96);interpft(X(end/2+1:end,:),96)];
end

[x,y] = o.getXY(X);
N = size(x,1);
[Dx,Dy] = o.getDXY(X);
length = sum(sqrt(Dx.^2 + Dy.^2))*2*pi/N;
area = sum(x.*Dy - y.*Dx)*pi/N;

reducedArea = 4*pi*area./length.^2;

end % geomProp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,nv] = initConfig(o,N,varargin)       
% [X,nv] = initConfig(n,varargin) returns N coordinates of boundary
% points.
% X = BOUNDARY(N,OPTIONS) can be used to call for different
% configuration.  Available OPTIONS are
%
%   'nv'          - (followed by) number of vesicles (in 
%                   linear configuration),
%   'angle'       - inclination angle of the vesicle(s) form the vertical
%                   position, 
%   'curly'       - returns a somehow wiggly vesicle.
%   'couette'     - the vesicles inside the default geometry of 
%                   couette flow (confined). 
%   'scale'       - multiplies the size of the output boundary
%   'choke'       - returns a choked domain.  Usually a solid boundary
%
%   'couette'     - returns a domain for a couette apparatus.
%   'figureEight' - returns a domain that resembles a pinched ellipse
%   'tube'        - returns a domain that is an ellongated ellipse

% EXAMPLE: X = boundary(64,'nv',3,'theta',pi/6);
%

options = varargin;
X = [];

% number of vesicles
if(any(strcmp(options,'nv')))
  nv = options{find(strcmp(options,'nv'))+1};
else
  nv = 1;
end

% rotate the vesicle by a certain angle
if(any(strcmp(options,'angle')))
  theta = options{find(strcmp(options,'angle'))+1};
else
  theta = zeros(nv,1);
end

% pick the center of the vesicles
if(any(strcmp(options,'center')))
  cen = options{find(strcmp(options,'center'))+1};
else
  cen = [(0:nv-1);zeros(1,nv)];
end

% desired reduced area
if(any(strcmp(options,'reducedArea')))
  ra = options{find(strcmp(options,'reducedArea'))+1};
else 
  ra = [];
end

% scale the size of the vesicle (reduced area is invariant under
% scale)
if(any(strcmp(options,'scale')))
  scale = options{find(strcmp(options,'scale'))+1};
else
  scale = 1;
end


t = (0:N-1)'*2*pi/N;
% Discretization in parameter space

if any(strcmp(options,'curly'))
  a = 1; b = 3*a; c = 0.85; 
  r = 0.5*sqrt( (a*cos(t-theta)).^2 + (b*sin(t-theta)).^2) + ...
      .07*cos(12*(t-theta));
  x = scale*c*r.*cos(t);
  y = scale*r.*sin(t);
%   tNew = o.arcLengthParamter(x,y);
%   r = 0.5*sqrt( (a*cos(tNew-theta)).^2 + (b*sin(tNew-theta)).^2) + ...
%       .07*cos(12*(tNew-theta));
%   x = scale*c*r.*cos(tNew);
%   y = scale*r.*sin(tNew);
  X = [x;y];
  % radius of curly vesicle

elseif any(strcmp(options,'cShape'))
  x = [-(1.5+sin(t)).*(cos(0.99*pi*cos(t)))];
  y = [+(1.5+sin(t)).*(sin(0.99*pi*cos(t)))];
  X0 = [x;y];

elseif any(strcmp(options,'crazy'))
  mollify = exp(1./((t/pi-1).^2-1));
  mollify(1) = 0;
  mollify(end) = 0;
  % mollifier supported on [0,2*pi] and centered at pi

  index = ceil(N/8);
  mollify1 = [mollify(index+1:end);mollify(1:index)];
  % shift the center of the mollifier

  r1 = .8*cos(16*t) + .16*sin(20*t);
  r1 = r1.*mollify1;
  % radius with additional high frequencies 'wiggles'

  index = ceil(N/2);
  mollify2 = [mollify(index+1:end);mollify(1:index)];
  % shift the center of the mollifier
  r2 = 2.4*sin(7*t) + .64*cos(10*t);
  r2 = r2.*mollify2;
  % radius with additional high frequencies 'wiggles'

  X = [8*cos(t) + r1.*cos(t) + r2.*cos(t);...
      sin(t) + r1.*sin(t) + r2.*sin(t)];
  % A somewhat complicated geometry

elseif any(strcmp(options,'star'))
  if any(strcmp(options,'folds'))
    folds = options{find(strcmp(options,'folds'))+1};
  else
    folds = 4;
  end
  radius = 1 + 0.3*cos(folds*t);
  X = [radius.*cos(t);radius.*sin(t)];
  % a star that comes very close to intersecting itself
  % at the origin

elseif any(strcmp(options,'openStar'))
 if any(strcmp(options,'folds'))
    folds = options{find(strcmp(options,'folds'))+1};
  else
    folds = 4;
  end
  if any(strcmp(options,'amplitude'))
    amp = options{find(strcmp(options,'amplitude'))+1};
  else
    amp = 5e-2;
  end
  radius = 1 + 0.2*cos(folds*t) + amp*cos(30*t);
  x = radius.*cos(t);
  y = radius.*sin(t);
  tNew = o.arcLengthParameter(x,y);
  % interpolate inverse function to find new parameter spacing

  radius = 1 + 0.2*cos(folds*tNew) + amp*cos(30*tNew);
  x = radius.*cos(tNew);
  y = radius.*sin(tNew);
  % re-define geometry equi-spaced in arclength

  X = [x;y];
  % a star that comes very close to intersecting itself at the origin

elseif any(strcmp(options,'wigglyStar'))
%  radius = 1 + .8*cos(2*t);
%  X = [radius.*cos(t);radius.*sin(t)] + [0.01*cos(80*t).*cos(t);0.01*sin(80*t).*sin(t)];
  % a star with only two folds, but a high frequency
  % added on top
  radius = 1 + 0.5*cos(3*t) + 0.05*cos(30*t);
  X = [radius.*cos(t);radius.*sin(t)];

elseif any(strcmp(options,'spikes'))
  % This shape is not implemented very well.  The curve
  % is too crazy
  N = numel(t);
  uprate = 5;
  t1 = t;
  t2 = (uprate*N-1:-1:0)'*2*pi/(uprate*N);
  % want to upsample the top part of the vesicle since it
  % has much sharper features

  x = [t1;t2];
  x = x(1:(uprate+1):end);
  y = [0.3*cos(10*t1)+0.7;ones(size(t2))+0.1*(cos(2*t2)+1.1)];
  y = y(1:(uprate+1):end);

  width = 0.2;
  cen = [0.5 2 3 3.8 5];
  height = [5 9 13 9 5];
  for k = 1:numel(cen)
    s = find(x > cen(k)-width & x < cen(k)+width & y > 1.01);
    y(s) = y(s) + height(k)*cos(pi/2/width*(x(s)-cen(k)));
  end

  z = x+1i*y;
  z = z - 1/2*(max(x) + min(x)) - 1i/2*(max(y) + min(y));
  z1 = z;
  zh = fftshift(fft(z));
  modes = (-N/2+1:N/2)';
  zh(abs(modes) > 100) = 0;
  z = ifft(ifftshift(zh));
  x = real(z); y = imag(z);
  X = [x;y];

elseif any(strcmp(options,'choke'))
  a = 10; b = 3; c = 0.6; order = 8;
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
%  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

elseif any(strcmp(options,'choke2'))
  a = 40; b = 3; c = 0.6; order = 8;
  % parameters for the boundary
%  Nsides = ceil(0.5*b/(2*a+2*b)*N);
%  Ntop = (N-4*Nsides)/2;
%  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
%  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
%  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
%  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
%  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
%  t = [t1;t2;t3;t4;t5];
  t = (0:N-1)'*2*pi/N;
  % Parameterize t so that geometry is closer to equi-spaced in arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);

  tNew = arcLengthParameter(o,x,y);
  r = (cos(tNew).^order + sin(tNew).^order).^(-1/order);
  x = a*r.*cos(tNew); y = b*r.*sin(tNew);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);

  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

elseif any(strcmp(options,'doublechoke'))
  a = 10; b = 3; c = 0.6; order = 8;
  shift = pi/2 + 0.1;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x-shift) < pi/2;
  y(ind) = y(ind).*(1-c*cos(2*(x(ind)-shift)))/(1+c);
  ind = abs(x+shift) < pi/2;
  y(ind) = y(ind).*(1-c*cos(2*(x(ind)+shift)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity.  shift controls the distance between the two
  % regions where the domain is restricted

elseif any(strcmp(options,'couette'))
  x = [20*cos(t)+cen(1,1) 10*cos(-t)+cen(1,2)];
  y = [20*sin(t)+cen(2,1) 10*sin(-t)+cen(2,2)];
  X = o.setXY(x,y);
  % annular domain
  
elseif any(strcmp(options,'couetteOuter'))
  x = [20*cos(t)+cen(1,1) 10*cos(-t)+cen(1,2)];
  y = [20*sin(t)+cen(2,1) 10*sin(-t)+cen(2,2)];
  X = o.setXY(x,y);
  % annular domain  

elseif any(strcmp(options,'doubleCouette'))
  x = [20*cos(t)+cen(1,1) 5*cos(-t)+cen(1,2) 5*cos(-t)+cen(1,3)];
  y = [20*sin(t)+cen(2,1) 5*sin(-t)+cen(2,2) 5*sin(-t)+cen(2,3)];
  X = o.setXY(x,y);
  % annular domain of genus 2

elseif any(strcmp(options,'quadCouette'))
  x = [20*cos(t)+cen(1,1) 5*cos(-t)+cen(1,2) 5*cos(-t)+cen(1,3) ...
                          5*cos(-t)+cen(1,4) 5*cos(-t)+cen(1,5)];
  y = [20*sin(t)+cen(2,1) 5*sin(-t)+cen(2,2) 5*sin(-t)+cen(2,3) ...
                          5*sin(-t)+cen(2,4) 5*sin(-t)+cen(2,5)];
  X = o.setXY(x,y);
  % annular domain of genus 4

elseif any(strcmp(options,'doubleFlower'))
  r = 17 + 2*cos(7*t);
  x = [r.*cos(t)+cen(1,1) 5*cos(-t)+cen(1,2) 5*cos(-t)+cen(1,3)];
  y = [r.*sin(t)+cen(2,1) 5*sin(-t)+cen(2,2) 5*sin(-t)+cen(2,3)];
  X = o.setXY(x,y);
  % annular domain of genus 2

elseif any(strcmp(options,'cylinder'))
  x = [20*scale*cos(t)+cen(1,1)];
  y = [20*scale*sin(t)+cen(2,1)];
  X = o.setXY(x,y);
  % single cylinder

elseif any(strcmp(options,'figureEight'))
  r = 1 + .7*cos(2*t);
  x = r.*cos(t); x = 2*x/max(x);
  y = r.*sin(t); y = y/max(y); 
  X = o.setXY(x,y);
  % pinched elliptical cylinder

elseif any(strcmp(options,'tube'))
  a = 4; b = 0.5; order = 10;
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

elseif any(strcmp(options,'diffuser'))
  a = 10; b = 1; c = 4; order = 8;
  % parameters for the boundary
%  Nsides = ceil(0.5*b/(2*a+2*b)*N);
%  Ntop = (N-4*Nsides)/2;
%  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
%  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
%  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
%  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
%  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
%  t = [t1;t2;t3;t4;t5];
  t = (0:N-1)'*2*pi/N;
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
%  ind = abs(x) < pi;
%  y(ind) = y(ind) + c/2*sign(y(ind)).*exp(-c./(pi^2 - x(ind).^2));
  ind = (abs(x) < 0.5 & y < 0);
  c = 0.5;
  y(ind) = y(ind) + 10*c/2.*exp(-c./(0.5^2 - x(ind).^2));

  X = [x;y];
  % opposite of chocked domain.  a and b control the length and
  % height.  c controls the width of the gap, and order controls the
  % regularity

elseif any(strcmp(options,'microfluidic'));
  t = (0:N-1)'*2*pi/N;

  a = 2; b = 2; order = 12;
  % parameters for the outer boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  Xout = [a*r.*cos(t);b*r.*sin(t)];

  a = 0.85; b = 0.85; order = 12;
  % parameters for the inner boundary
  r = (cos(-t).^order + sin(-t).^order).^(-1/order);
  Xin1 = [a*r.*cos(-t)-0.95;b*r.*sin(-t)-0.95];
  Xin2 = [a*r.*cos(-t)-0.95;b*r.*sin(-t)+0.95];
  Xin3 = [a*r.*cos(-t)+0.95;b*r.*sin(-t)-0.95];
  Xin4 = [a*r.*cos(-t)+0.95;b*r.*sin(-t)+0.95];

  X = [Xout Xin1 Xin2 Xin3 Xin4];
  
elseif any(strcmp(options,'createVesforDLD'))
  % generate vesicles for DLD simulation
  oc = curve;
  %   build reference elliptical vesicle  
  if(any(strcmp(options,'refVesicle')))
    X0 = options{find(strcmp(options,'refVesicle'))+1};
  else
    X0 = o.ellipse(N,ra);  
    [~,area,~] = oc.geomProp(X0);
    scale = 2*scale/sqrt(area/pi);
    X0 = scale*X0;
  end

  % Randomly place vesicles in the given DLD domain
  if any(strcmp(options,'randomlyPlaceVes'))
    ind = find(strcmp(options,'randomlyPlaceVes'));
    XwallsExt = options{ind+1};
    XwallsInt = options{ind+2};
    xrange = options{ind+3};
    yrange = options{ind+4};
    ttObj = options{ind+5};
    domainType = options{ind+6};
    if strcmp(domainType,'DLDWhole')
      prams = options{ind+7};
      Dpost = prams.Dpost; Dy = prams.Dy; epsilon = prams.epsilon; 
      ncol = prams.ncol; xfreeze = prams.xfreeze;
      volFrac = options{ind+8};
    else
      Dpost = []; Dy = []; epsilon = []; ncol = []; xfreeze = [];
      volFrac = [];
    end
    
    [X,~] = oc.fillDomain(X0,'domain',domainType,'numVes',nv,...
        'wallsExt',XwallsExt,'wallsInt',XwallsInt,'volFrac',volFrac,...
        'tstepObj',ttObj,'xrange',xrange,'yrange',yrange,'Dpost',Dpost,...
        'Dy',Dy,'epsilon',epsilon,'ncol',ncol,'xfreeze',xfreeze);
    
  else
    scale = 2*scale/sqrt(area/pi);
    X0 = scale*X0;    
  end
  % end of building reference vesicles.  Only need to rotate
  % and shift them as desired
  
else
  if ~isempty(ra)
    X0 = o.ellipse(N,ra);
    % build a vesicle of reduced area ra with N points
  else
%     a = 1; b = 4.8*a; c = 1; 
%     r = 0.5*sqrt( (a*cos(t)).^2 + (b*sin(t)).^2);
%     x = c*r.*cos(t); y = r.*sin(t);
    
    % UNCOMMENT IF YOU WANT TO START W/ POINTS EQUI-DIST IN ARC-LENGTH
    %tNew = o.arcLengthParameter(x,y);
    %r = 0.5*sqrt( (a*cos(tNew)).^2 + (b*sin(tNew)).^2);
    %x = c*r.*cos(tNew); y = r.*sin(tNew);
    
%     X0 = [x;y];
    % this shape has reduced area very close to 0.65
    X0 = o.ellipse(N,0.65);
  end
  oc = curve;
  
  [ra,area,length] = oc.geomProp(X0);
  
  scale = scale/sqrt(area/pi);
  X0 = scale*X0;
%   build reference elliptical vesicle

  if any(strcmp(options,'volFrac'))
    ind = find(strcmp(options,'volFrac'));
    volFrac = options{ind+1};
    Xwalls = options{ind+2};
    ttObj = options{ind+3};
    X = oc.fillDomain(X0,'domain','couette','volFrac',volFrac,'walls',Xwalls,'tstepObj',ttObj);
    nv = size(X,2);
  else
    scale = 2*scale/sqrt(area/pi);
    X0 = scale*X0;
  end

  if size(cen,2) ~= nv
    b = 3*max(X0(1:N));
    cen = [b*(0:nv-1);zeros(1,nv)];
  end
  % centers if there are multiple vesicles.
end 
% end of building reference vesicles.  Only need to rotate
% and shift them as desired

if isempty(X)
  % if X has not been defined, we only have X0 which is 
  % the reference vesicle which not needs to be
  % rotated and translated
  X = zeros(2*N,nv);
  for k=1:nv
    X(1:N,k) = cos(theta(k)) * X0(1:N) - ...
      sin(theta(k)) * X0(N+1:2*N) + cen(1,k);
    X(N+1:2*N,k) = sin(theta(k)) * X0(1:N) +  ...
      cos(theta(k)) * X0(N+1:2*N) + cen(2,k);
  end
  % Rotate vesicles as requested 
end

end % initConfig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X0 = cElegans(o,N,a,b,px,py)
u = linspace(-pi,pi,N)';
R = [a, b];
P = [px, py];
M = 2./P;


x = R(1).*sign(cos(u)).*abs(cos(u)).^M(1);
y = R(2).*sign(sin(u)).*abs(sin(u)).^M(2);

X0 = [x(1:end-1); y(1:end-1)];
    
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X0 = ellipse(o,N,ra)
% X0 = o.ellipse(N,ra) finds the ellipse (a*cos(theta),sin(theta)) so
% that the reduced area is ra.  Uses N points.  Parameter a is found 
% by using bisection method

t = (0:N-1)'*2*pi/N;
a = (1 - sqrt(1-ra^2))/ra;
% initial guess using approximation length = sqrt(2)*pi*sqrt(a^2+1)

X0 = [a*cos(t);sin(t)];

[raNew,~,~] = o.geomProp(X0);

cond = abs(raNew - ra)/ra < 1e-4;
maxiter = 10;
iter = 0;
while (~cond && iter < maxiter);
  iter = iter + 1;
  if (raNew > ra)
    a = 0.9 * a;
  else
    a = 1.05 * a;
  end
  % update the major axis
  X0 = [cos(t);a*sin(t)];
  % Compute new possible configuration
  [raNew,~,~] = o.geomProp(X0);
  % compute new reduced area
  cond = abs(raNew - ra) < 1e-2;
  % check if the residual is small enough
end
% iteration quits if reduced area is achieved within 1% or 
% maxiter iterations have been performed

% UNCOMMENT IF YOU WANT TO START W/ POINTS EQUI-DIST IN ARC-LENGTH
% tNew = o.arcLengthParameter(X0(1:N),X0(N+1:2*N));
% X0 = [cos(tNew);a*sin(tNew)];


end % ellipse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XwallsInt,XwallsExt,L,H,x_freeze,radiusRight,radiusUp,...
        Dpostx,Dposty] = ...
        initConfigDLD(o,itype,Nint,Next,Xpost,Dpost,epsilon,...
        periodN,Dx,Dy,nrow,ncol)       
% [XwallsInt,XwallsExt,L,H,x_freeze] = ...
% initConfigDLDOptim(o,Nint,Next,Xpost,varargin) builds the initial
% configuration of the DLD device with the given post geometry

if strcmp(itype,'wholeDLD')
  % If Xpost is not given, then suppose it is a circle.
  if isempty(Xpost)
    Xpost = [Dpost/2*cos(-(0:Nint-1)'/Nint*2*pi);...
        Dpost/2*sin(-(0:Nint-1)'/Nint*2*pi)];  
    Dpostx = Dpost;
    Dposty = Dpost;
    xup = interpft(Xpost(1:end/2),256); 
    yup = interpft(Xpost(end/2+1:end),256);
    radiusRight = Dpost/2;
    radiusUp = Dpost/2;
    
    maxDist2Left = Dpost/2;
  else
    % if we fit the post in a rectangle, then we put the posts as if we 
    % put the corresponding rectangular posts in DLD
    xup = interpft(Xpost(1:end/2),512); 
    yup = interpft(Xpost(end/2+1:end),512);
      
    % distances between top and bottom, left and right of the posts
    Dpostx = max(xup)-min(xup);
    Dposty = max(yup)-min(yup);
      
    % Move shape such that the rectangle around the shape is centered
    x = Xpost(1:end/2); y = Xpost(end/2+1:end);
    x2 = x-(max(xup)-Dpostx/2);
    y2 = y-(max(yup)-Dposty/2);
    Xpost = [x2;y2];
    % radiusLeft=radiusRight, and same for the top-bottom, since the
    % post's and rectangle's centers coincides
%       radiusLeft = mean(xup)-min(xup); radiusRight = max(xup)-mean(xup);
%       radiusUp = max(yup)-mean(yup); radiusDown = mean(yup)-min(yup);
    radiusLeft = Dpostx/2; radiusRight = radiusLeft;
    radiusUp = Dposty/2; radiusDown = radiusUp;
    % this resembles, Ranjan,  Zeming - Zhang, LoC (2014) paper in which 
    % they studied separation for several post shapes
    
    % Find where the smallest gap size is, place vesicles there and scale the
    x2up = interpft(x2,512);
    y2up = interpft(y2,512);

    % look at the points in the horizontal direction and find the largest size
    % across the shape in y-direction
    maxGsize = -1E+4;
    sort_x2 = sort(x2up);
    for ij = 1 : numel(x2up)
      ids = find(abs(x2up-sort_x2(ij))<=1e-2);
      ys = y2up(ids);
      gapSize = max(ys)-min(ys);

      if gapSize>=maxGsize 
        maxGsize = gapSize;
        [~,idx] = max(y2up(ids));
        maxDist2Left = x2up(ids(idx))+Dpostx/2;  
      end
    end

  end
  % suppose same separation in both directions if Dx is not given
  if isempty(Dx)
    Dx = Dy;    
  end
  
  % Height and length of the exterior wall
  H = (nrow+1)*Dy + nrow*Dposty; % necessary for large periodN
  L = (ncol+1)*(Dx+Dpostx);

  % parameters for the exterior wall
  a = L/2;
  b = H/2; order = 20;


  % a and b control length and height of exterior wall.
  % Order controls the regularity.
  t = (0:Next-1)'*2*pi/Next;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  XwallsExt = [x;y];

  % Row shift  
  delLat = (Dy+Dposty)*epsilon;

  % Generate more rows than asked and remove ones outside exterior wall
  pnrow = ceil(4*nrow);
  pHeight = (pnrow-1)*Dy + (pnrow-1)*Dposty;

  % Place horizontally the obstacles in the first column
  centy1stCol = linspace(-pHeight/2,pHeight/2,pnrow)';
  centx = linspace(-L/2+1.5*Dpostx-maxDist2Left+Dx,L/2-Dpostx/2-maxDist2Left-Dx,ncol);

  % when to freeze vesicles (after slightly after one period)
  if numel(centx)>=ceil(periodN)+2
    x_freeze = centx(1*ceil(periodN)+2)+Dpostx/2;
  else
    x_freeze = [];
  end

  % Place vertically all the obstacles
  centx = centx(ones(pnrow,1),:);

  % Shift centers of the obstacles
  delLatVect = [0:ncol-1]*delLat;
  centy = delLatVect(ones(pnrow,1),:) + centy1stCol(:,ones(1,ncol));

  % Place the obstacles
  Xobs = zeros(2*Nint,pnrow*ncol);
  for iwall = 1 : pnrow*ncol
    Xobs(:,iwall) = [Xpost(1:end/2)+centx(iwall);...
        Xpost(end/2+1:end)+centy(iwall)];
  end


  % Remove obstacles outside the geometry, replace the ones intersecting the
  % geometry with ellipses
  tobs = [0:Nint-1]'*2*pi/Nint;
  jdx = 1;
  for iwall = 1 : pnrow*ncol
    xwall = Xobs(1:end/2,iwall); ywall = Xobs(end/2+1:end,iwall);
    idx1 = isempty(find(abs(xwall) >= 0.98*L,1));
    idx2 = isempty(find(abs(ywall) >= 0.98*H/2,1));
    if idx1 && idx2
      XwallsInt(:,jdx) = [xwall;ywall];
      jdx = jdx + 1;
    else
      centx = mean(xwall);
      if isempty(find(ywall<0,1))
        minyloc = min(ywall);
        smallDiam = 0.98*H/2-minyloc;
        if smallDiam >= 0.25
%           widthX = sqrt((Dposty/2)^2-(Dposty/2-smallDiam/2)^2);    
          widthX = Dpostx/2;
          centy = minyloc + smallDiam/2;
          sxwalls = widthX*cos(-tobs)+centx;
          sywalls = smallDiam/2*sin(-tobs)+centy;
          sidx = isempty(find(sywalls >= 0.99*H/2,1));
          if sidx
            XwallsInt(:,jdx) = [sxwalls;sywalls];
            jdx = jdx + 1;
          end
        end % smallDiam
      else
        maxyloc = max(ywall);
        smallDiam = maxyloc+0.98*H/2;
        if smallDiam >= 0.25
%           widthX = sqrt((Dposty/2)^2-(Dposty/2-smallDiam/2)^2);    
          widthX = Dpostx/2;
          centy = maxyloc - smallDiam/2;
          sxwalls = widthX*cos(-tobs)+centx;
          sywalls = smallDiam/2*sin(-tobs)+centy;
          sidx = isempty(find(sywalls <= -0.99*H/2,1));
          if sidx
            XwallsInt(:,jdx) = [sxwalls;sywalls];
            jdx = jdx + 1;
          end % sidx
        end % smallDiam
      end % find(ywall<0)
    end %idx1&&idx2
  end % iwall  
    
elseif strcmp(itype,'optimDLD')
  % Smaller DLD device for optimization (nrow x ncol array of a DLD with 
  % a period p)  
  % If Xpost is not given, then suppose it is a circle.
  if isempty(Xpost)
    Xpost = [Dpost/2*cos(-(0:Nint-1)'/Nint*2*pi);...
        Dpost/2*sin(-(0:Nint-1)'/Nint*2*pi)];  
    Dpostx = Dpost;
    Dposty = Dpost;
    xup = interpft(Xpost(1:end/2),256); 
    yup = interpft(Xpost(end/2+1:end),256);
    centx = mean(xup); centy = mean(yup); 
    radiusRight = Dpost/2;
    radiusUp = Dpost/2;
  else
    if 0
      % Get the distance between the center of the post and the point at x = 0,
      % y<0 and x=0, y>0
      xup = interpft(Xpost(1:end/2),512); 
      yup = interpft(Xpost(end/2+1:end),512);
      centx = mean(xup); centy = mean(yup);
      idsAtCentx = abs(xup-centx)<0.1; idsUp = yup>centy; idsDown = yup<centy;
      idsAtCentxUp = idsAtCentx&idsUp; idsAtCentxDown = idsAtCentx&idsDown;
      radiusDown = sqrt((xup(idsAtCentxDown)-centx).^2+...
          (yup(idsAtCentxDown)-centy).^2);
      radiusUp = sqrt((xup(idsAtCentxUp)-centx).^2+...
          (yup(idsAtCentxUp)-centy).^2);

      radiusDown = max(radiusDown); radiusUp = max(radiusUp);

      idsLeft = xup<centx;
      idsAtCenty = abs(yup-centy)<0.1;
      idsAtCentyLeft = idsAtCenty&idsLeft;
      radiusLeft = sqrt((xup(idsAtCentyLeft)-centx).^2+...
        (yup(idsAtCentyLeft)-centy).^2);
      radiusLeft = max(radiusLeft);

      idsRight = xup>centx;
      idsAtCentyRight = idsAtCenty&idsRight;
      radiusRight = sqrt((xup(idsAtCentyRight)-centx).^2+...
        (yup(idsAtCentyRight)-centy).^2);
      radiusRight = max(radiusRight);
      % radiusUp is the vertical distance between top of the device (at x = centx) and the
      % center of the device, likewise, radiusDown is the distance to the bottom
      % of the device. This is essential to keep the gap size same as the shape
      % of the device changes.
      Dpostx = radiusLeft + radiusRight;
      Dposty = radiusDown + radiusUp;
    else
      % if we fit the post in a rectangle, then we put the posts as if we 
      % put the corresponding rectangular posts in DLD
      xup = interpft(Xpost(1:end/2),512); 
      yup = interpft(Xpost(end/2+1:end),512);
      
      % distances between top and bottom, left and right of the posts
      Dpostx = max(xup)-min(xup);
      Dposty = max(yup)-min(yup);
      
      % Move shape such that the rectangle around the shape is centered
      x = Xpost(1:end/2); y = Xpost(end/2+1:end);
      x2 = x-(max(xup)-Dpostx/2);
      y2 = y-(max(yup)-Dposty/2);
      Xpost = [x2;y2];
      % radiusLeft=radiusRight, and same for the top-bottom, since the
      % post's and rectangle's centers coincides
%       radiusLeft = mean(xup)-min(xup); radiusRight = max(xup)-mean(xup);
%       radiusUp = max(yup)-mean(yup); radiusDown = mean(yup)-min(yup);
      radiusLeft = Dpostx/2; radiusRight = radiusLeft;
      radiusUp = Dposty/2; radiusDown = radiusUp;
      % this resembles, Ranjan,  Zeming - Zhang, LoC (2014) paper in which 
      % they studied separation for several post shapes
    end

  end % if isempty(Xpost)

  % suppose same separation in both directions
  if isempty(Dx)
    Dx = Dy;    
  end

  % Height and length of the exterior wall
  H = (nrow-1)*Dy + nrow*Dposty; % necessary for large periodN
  L = (ncol+1)*(Dx+Dpostx);
  
  % Exterior long tube
  a = L/2; b = H/2; order = 20;
  t = (0:Next-1)'*2*pi/Next;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  XwallsExt = [x;y];
  

  % Row shift  
  delLat = (Dy+Dposty)*epsilon;

  % Generate more rows than asked and remove ones outside exterior wall
  pnrow = ceil(nrow^2); 
  pHeight = (pnrow-1)*Dy + (pnrow-1)*Dposty;

  % Place horizontally the obstacles in the first column
  centy1stCol = linspace(-pHeight/2,pHeight/2,pnrow)';
  centx = linspace(-L/2+Dpostx+Dx,L/2-Dpostx-Dx,ncol);

  % when to freeze vesicles (above the third column)
  % this is used for optimization
  x_freeze = centx(3)+1.5*radiusRight;
  
  % Place vertically all the obstacles
  centx = centx(ones(pnrow,1),:);

  % Shift centers of the obstacles
  delLatVect = [0:ncol-1]*delLat;
  centy = delLatVect(ones(pnrow,1),:) + centy1stCol(:,ones(1,ncol));

  % Place the obstacles
  Xobs = zeros(2*Nint,pnrow*ncol);
  for iwall = 1 : pnrow*ncol
    Xobs(:,iwall) = [Xpost(1:end/2)+centx(iwall);...
        Xpost(end/2+1:end)+centy(iwall)];
  end


  % Remove obstacles outside the geometry, replace the ones intersecting the
  % geometry with ellipses
  tobs = [0:Nint-1]'*2*pi/Nint;
  jdx = 1;
  for iwall = 1 : pnrow*ncol
    xwall = Xobs(1:end/2,iwall); ywall = Xobs(end/2+1:end,iwall);
    idx1 = isempty(find(abs(xwall) >= 0.98*L,1));
    idx2 = isempty(find(abs(ywall) >= 0.98*H/2,1));
    if idx1 && idx2
      XwallsInt(:,jdx) = [xwall;ywall];
      jdx = jdx + 1;
    else
      centx = mean(xwall);
      if isempty(find(ywall<0,1))
        minyloc = min(ywall);
        smallDiam = 0.98*H/2-minyloc;
        if smallDiam >= 0.1*Dposty
%           widthX = sqrt((Dposty/2)^2-(Dposty/2-smallDiam/2)^2);    
          widthX = Dpostx/2;
          centy = minyloc + smallDiam/2;
          sxwalls = widthX*cos(-tobs)+centx;
          sywalls = smallDiam/2*sin(-tobs)+centy;
          sidx = isempty(find(sywalls >= 0.99*H/2,1));
          if sidx
            XwallsInt(:,jdx) = [sxwalls;sywalls];
            jdx = jdx + 1;
          end
        end % smallDiam
      else
        maxyloc = max(ywall);
        smallDiam = maxyloc+0.98*H/2;
        if smallDiam >= 0.1*Dposty
%           widthX = sqrt((Dposty/2)^2-(Dposty/2-smallDiam/2)^2);    
          widthX = Dpostx/2;
          centy = maxyloc - smallDiam/2;
          sxwalls = widthX*cos(-tobs)+centx;
          sywalls = smallDiam/2*sin(-tobs)+centy;
          sidx = isempty(find(sywalls <= -0.99*H/2,1));
          if sidx
            XwallsInt(:,jdx) = [sxwalls;sywalls];
            jdx = jdx + 1;
          end % sidx
        end % smallDiam
      end % find(ywall<0)
    end %idx1&&idx2
  end % iwall  
  
end % if itype
end % initConfigDLD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XwallsInt,XwallsExt,L,H,x_freeze,radiusRight,radiusUp,...
        Dpostx1,Dposty1,Dpostx2,Dposty2] = ...
        initConfigDLD_multiPillar(o,Nint,Next,Xpost1,Xpost2,...
        epsilon,periodN,Dx,Dy1,nrow,ncol)       
% ASSUMING BOTH PILLARS HAVE THE SAME SIZE
% if we fit the post in a rectangle, then we put the posts as if we 
% put the corresponding rectangular posts in DLD
xup = interpft(Xpost1(1:end/2),512); 
yup = interpft(Xpost1(end/2+1:end),512);

% distances between top and bottom, left and right of the posts
Dpostx1 = max(xup)-min(xup);
Dposty1 = max(yup)-min(yup);

% Move shape such that the rectangle around the shape is centered
x = Xpost1(1:end/2); y = Xpost1(end/2+1:end);
x2 = x-(max(xup)-Dpostx1/2);
y2 = y-(max(yup)-Dposty1/2);
Xpost1 = [x2;y2];


% if we fit the post in a rectangle, then we put the posts as if we 
% put the corresponding rectangular posts in DLD
xup = interpft(Xpost2(1:end/2),512); 
yup = interpft(Xpost2(end/2+1:end),512);

% distances between top and bottom, left and right of the posts
Dpostx2 = max(xup)-min(xup);
Dposty2 = max(yup)-min(yup);

% Move shape such that the rectangle around the shape is centered
x = Xpost2(1:end/2); y = Xpost2(end/2+1:end);
x2 = x-(max(xup)-Dpostx2/2);
y2 = y-(max(yup)-Dposty2/2);
Xpost2 = [x2;y2];

% radiusLeft=radiusRight, and same for the top-bottom, since the
% post's and rectangle's centers coincides
%       radiusLeft = mean(xup)-min(xup); radiusRight = max(xup)-mean(xup);
%       radiusUp = max(yup)-mean(yup); radiusDown = mean(yup)-min(yup);
radiusLeft = Dpostx1/2; radiusRight = radiusLeft;
radiusUp = Dposty1/2; radiusDown = radiusUp;
% this resembles, Ranjan,  Zeming - Zhang, LoC (2014) paper in which 
% they studied separation for several post shapes

% Find where the smallest gap size is, place vesicles there and scale the
x1up = interpft(Xpost1(1:end/2),512);
y1up = interpft(Xpost1(end/2+1:end),512);

x2up = interpft(Xpost2(1:end/2),512);
y2up = interpft(Xpost2(end/2+1:end),512);

% look at the points in the horizontal direction and find the largest size
% across the shape in y-direction
maxGsize1 = -1E+4;
sort_x1 = sort(x1up);

maxGsize2 = -1E+4;
sort_x2 = sort(x2up);

for ij = 1 : numel(x2up)
  ids = find(abs(x1up-sort_x1(ij))<=1e-2);
  ys = y1up(ids);
  gapSize1 = max(ys)-min(ys);

  if gapSize1>=maxGsize1 
    maxGsize1 = gapSize1;
    [~,idx] = max(y1up(ids));
    maxDist2Left1 = x1up(ids(idx))+Dpostx1/2;  
  end  
    
    
  ids = find(abs(x2up-sort_x2(ij))<=1e-2);
  ys = y2up(ids);
  gapSize2 = max(ys)-min(ys);

  if gapSize2>=maxGsize2 
    maxGsize2 = gapSize2;
    [~,idx] = max(y2up(ids));
    maxDist2Left2 = x2up(ids(idx))+Dpostx2/2;  
  end
end

% Height and length of the exterior wall
H = (nrow+1)*Dy1 + nrow*Dposty1; % necessary for large periodN
L = (ncol+1)*Dx + (ncol+1)/2 * (Dpostx1+Dpostx2);
% parameters for the exterior wall
a = L/2;
b = H/2; order = 20;


% a and b control length and height of exterior wall.
% Order controls the regularity.
t = (0:Next-1)'*2*pi/Next;
r = (cos(t).^order + sin(t).^order).^(-1/order);
x = a*r.*cos(t); y = b*r.*sin(t);
XwallsExt = [x;y];

% Row shift  
delLat = (Dy1+Dposty1)*epsilon;

% Generate more rows than asked and remove ones outside exterior wall
pnrow = ceil(4*nrow);
pHeight = (pnrow-1)*Dy1 + (pnrow-1)*Dposty1;

% Place horizontally the obstacles in the first column
centy1stCol = linspace(-pHeight/2,pHeight/2,pnrow)';
centx = linspace(-L/2+Dpostx1+Dpostx2/2-maxDist2Left1+Dx,L/2-Dpostx2/2-maxDist2Left1-Dx,ncol);

% when to freeze vesicles (after slightly after one period)
if numel(centx)>=ceil(periodN)+2
x_freeze = centx(1*ceil(periodN)+2)+Dpostx1/2;
else
x_freeze = [];
end

% Place vertically all the obstacles
centx = centx(ones(pnrow,1),:);

% Shift centers of the obstacles
delLatVect = [0:ncol-1]*delLat;
centy = delLatVect(ones(pnrow,1),:) + centy1stCol(:,ones(1,ncol));

% Place the obstacles
Xobs = zeros(2*Nint,pnrow*ncol);
iwall = 1;
for irow = 1 : pnrow
  for icol = 1 : ncol
    
    if rem(icol,2)==0
      Xobs(:,iwall) = [Xpost1(1:end/2)+centx(irow,icol);...
          Xpost1(end/2+1:end)+centy(irow,icol)];
    else
      Xobs(:,iwall) = [Xpost2(1:end/2)+centx(irow,icol);...
          Xpost2(end/2+1:end)+centy(irow,icol)];  
    end
    
    iwall = iwall + 1;
  end
end


% Remove obstacles outside the geometry, replace the ones intersecting the
% geometry with ellipses
tobs = [0:Nint-1]'*2*pi/Nint;
jdx = 1;
for iwall = 1 : pnrow*ncol
xwall = Xobs(1:end/2,iwall); ywall = Xobs(end/2+1:end,iwall);
idx1 = isempty(find(abs(xwall) >= 0.98*L,1));
idx2 = isempty(find(abs(ywall) >= 0.98*H/2,1));
if idx1 && idx2
  XwallsInt(:,jdx) = [xwall;ywall];
  jdx = jdx + 1;
else
  centx = mean(xwall);
  if isempty(find(ywall<0,1))
    minyloc = min(ywall);
    smallDiam = 0.98*H/2-minyloc;
    if smallDiam >= 0.25
%           widthX = sqrt((Dposty/2)^2-(Dposty/2-smallDiam/2)^2);    
      Dpostx = max(xwall)-min(xwall);
      widthX = Dpostx/2;
      centy = minyloc + smallDiam/2;
      sxwalls = widthX*cos(-tobs)+centx;
      sywalls = smallDiam/2*sin(-tobs)+centy;
      sidx = isempty(find(sywalls >= 0.99*H/2,1));
      if sidx
        XwallsInt(:,jdx) = [sxwalls;sywalls];
        jdx = jdx + 1;
      end
    end % smallDiam
  else
    maxyloc = max(ywall);
    smallDiam = maxyloc+0.98*H/2;
    if smallDiam >= 0.25
%           widthX = sqrt((Dposty/2)^2-(Dposty/2-smallDiam/2)^2);
      Dpostx = max(xwall)-min(xwall);
      widthX = Dpostx/2;
      centy = maxyloc - smallDiam/2;
      sxwalls = widthX*cos(-tobs)+centx;
      sywalls = smallDiam/2*sin(-tobs)+centy;
      sidx = isempty(find(sywalls <= -0.99*H/2,1));
      if sidx
        XwallsInt(:,jdx) = [sxwalls;sywalls];
        jdx = jdx + 1;
      end % sidx
    end % smallDiam
  end % find(ywall<0)
end %idx1&&idx2
end % iwall  

end % initConfigDLD_multiPillar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xpost,flgFits] = samplePostSymm(o,Npost,DpostBounds,fewPoints)
% Xpost = samplePostSym(o,Npost,fewPoints,Dpost) builds a post shape using
% 4 control points Pi assuming a symmetry in the 
% post shape and 4th order b-splines.

% fewPoints(1): x coordinate of point 1 (x>0, y=0)
% fewPoints(2): x coordinate of point 2 (x>0, y>0)
% fewPoints(3): y coordinate of point 2 (x>0, y>0)
% fewPoints(4): y coordinate of point 3 (x=0, y>0)
% fewPoints(5): x coordinate of point 4 (x<0, y>0)
% fewPoints(6): y coordinate of point 4 (x<0, y>0)
% fewPoints(7): x coordinate of point 5 (x<0, y=0)
% fewPoints(8): angle


% build the symmetric points with respect to the x axis
points(1,:) = [fewPoints(1) 0];
points(2,:) = [fewPoints(2) fewPoints(3)];
points(3,:) = [0 fewPoints(4)];
points(4,:) = [-fewPoints(5) fewPoints(6)];
points(5,:) = [-fewPoints(7) 0];
points(6,:) = [-fewPoints(5) -fewPoints(6)];
points(7,:) = [0 -fewPoints(4)];
points(8,:) =[fewPoints(2) -fewPoints(3)];

alpha = fewPoints(8); % angle

% build a closed loop using 4th order b-splines
c = spcrv([points(end,:);points;points(1,:);points(2,:);points(3,:)].',5);


% reach the desired resolution Npost, remove the last three points which
% overlap with the first three points

values = [c(:,end-9:end) c(:,5:end-10)];

x = interpft(values(1,:),Npost)';
y = interpft(values(2,:),Npost)';

x = flipud(x); y = flipud(y);
Xpost = [x;y];

% Center it to the origin
Xpost(1:end/2) = Xpost(1:end/2)-mean(Xpost(1:end/2));
Xpost(end/2+1:end) = Xpost(end/2+1:end)-mean(Xpost(end/2+1:end));

XpostRotate = zeros(numel(Xpost),1);
XpostRotate(1:end/2,1) = cos(alpha)*Xpost(1:end/2)-sin(alpha)*Xpost(end/2+1:end);
XpostRotate(end/2+1:end,1) = sin(alpha)*Xpost(1:end/2)+cos(alpha)*Xpost(end/2+1:end);
Xpost = XpostRotate;

% redistribute the disc. points equally in arclength
for it = 1 : 100
[Xpost,~,~] = o.redistributeArcLength(Xpost);
end

% Center it to the origin
Xpost(1:end/2) = Xpost(1:end/2)-mean(Xpost(1:end/2));
Xpost(end/2+1:end) = Xpost(end/2+1:end)-mean(Xpost(end/2+1:end));

% Fit the post into a square with a length of Dpost
Dx = max(Xpost(1:end/2))-min(Xpost(1:end/2));
Dy = max(Xpost(end/2+1:end))-min(Xpost(end/2+1:end));
maxSize = abs(max(Dx,Dy));

% if larger than maxSize or smaller than minSize, penalize.
% DpostBounds(1) = upperBound; DpostBounds(2) = lowerBound;
flgFits = true; % flag to show if the post shape fits
if maxSize > DpostBounds(1) || maxSize < DpostBounds(2)
  flgFits = false;
end

end %samplePostSymm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xpost,Dx,Dy] = samplePost(o,Npost,points)
% Xpost = samplePost(o,Npost,points) builds a post shape using
% 8 control points Pi (each 1x2) and 4th order b-splines.
    
% build a closed loop using 4th order b-splines
%c = spcrv([points(8,:);points;points(1,:);points(2,:)].',4);
c = spcrv([points(end,:);points;points(1,:);points(2,:);points(3,:)].',5);
%CS = cscvn([points;points(1,:)].');
%[c,t] = fnplt(CS);

% reach the desired resolution Npost, remove the last three points which
% overlap with the first three points

values = [c(:,end-9:end) c(:,5:end-10)];

%x = interpft(c(1,2:end-2),Npost)';
%y = interpft(c(2,2:end-2),Npost)';

x = interpft(values(1,:),Npost)';
y = interpft(values(2,:),Npost)';

x = flipud(x); y = flipud(y);

% reparametrize to kill high frequencies

Xpost = [x;y];


% redistribute the disc. points so that corners are smooth
for it = 1 : 10
[Xpost,~] = o.reparametrize(Xpost,[],6,20);
end

% Center it to the origin
Xpost(1:end/2) = Xpost(1:end/2)-mean(Xpost(1:end/2));
Xpost(end/2+1:end) = Xpost(end/2+1:end)-mean(Xpost(end/2+1:end));

% Fit the post into a rectangle of Dx by Dy
Dx = max(Xpost(1:end/2))-min(Xpost(1:end/2));
Dy = max(Xpost(end/2+1:end))-min(Xpost(end/2+1:end));

end %samplePost

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XwallsInt,XwallsExt,L,H,radiusRight,radiusUp,...
        Dpostx,Dposty] = ...
        initConfigDLDRot(o,itype,Nint,Next,Xpost,Dpost,epsilon,Dx,Dy,nrow,ncol)       

% Smaller DLD device for optimization (nrow x ncol array of a DLD with 
% a period p)  
% If Xpost is not given, then suppose it is a circle.
if isempty(Xpost)
  Xpost = [Dpost/2*cos(-(0:Nint-1)'/Nint*2*pi);...
      Dpost/2*sin(-(0:Nint-1)'/Nint*2*pi)];  
  Dpostx = Dpost;
  Dposty = Dpost;
  xup = interpft(Xpost(1:end/2),256); 
  yup = interpft(Xpost(end/2+1:end),256);
  centx = mean(xup); centy = mean(yup); 
  radiusRight = Dpost/2;
  radiusUp = Dpost/2;
  
  maxDist2Left = Dpost/2;
else

  % if we fit the post in a rectangle, then we put the posts as if we 
  % put the corresponding rectangular posts in DLD
  xup = interpft(Xpost(1:end/2),512); 
  yup = interpft(Xpost(end/2+1:end),512);

  % distances between top and bottom, left and right of the posts
  Dpostx = max(xup)-min(xup);
  Dposty = max(yup)-min(yup);

  % Move shape such that the rectangle around the shape is centered
  x = Xpost(1:end/2); y = Xpost(end/2+1:end);
  x2 = x-(max(xup)-Dpostx/2);
  y2 = y-(max(yup)-Dposty/2);
  Xpost = [x2;y2];
  % radiusLeft=radiusRight, and same for the top-bottom, since the
  % post's and rectangle's centers coincides
  
  radiusLeft = Dpostx/2; radiusRight = radiusLeft;
  radiusUp = Dposty/2; radiusDown = radiusUp;
  % this resembles, Ranjan,  Zeming - Zhang, LoC (2014) paper in which 
  % they studied separation for several post shapes

  % Find where the smallest gap size is, place vesicles there and scale the
    x2up = interpft(x2,512);
    y2up = interpft(y2,512);

    % look at the points in the horizontal direction and find the largest size
    % across the shape in y-direction
    maxGsize = -1E+4;
    sort_x2 = sort(x2up);
    for ij = 1 : numel(x2up)
      ids = find(abs(x2up-sort_x2(ij))<=1e-2);
      ys = y2up(ids);
      gapSize = max(ys)-min(ys);

      if gapSize>=maxGsize 
        maxGsize = gapSize;
        [~,idx] = max(y2up(ids));
        maxDist2Left = x2up(ids(idx))+Dpostx/2;  
      end
    end
    
end % if isempty(Xpost)

% suppose same separation in both directions if Dx is not given
if isempty(Dx)
  Dx = Dy;    
end

if itype == 1 % exterior wall crosses centers of imaginary pillars
  % Height and length of the exterior wall
  H = (nrow+1)*(Dy+Dposty); 
  L = (ncol+1)*(Dx+Dpostx);
elseif itype == 2 % exterior wall is in the middle of the gap
  H = nrow*(Dy+Dposty);
  L = ncol*(Dx+Dpostx);
end

% Exterior long tube
a = L/2; b = H/2; order = 20;
t = (0:Next-1)'*2*pi/Next;
r = (cos(t).^order + sin(t).^order).^(-1/order);
x = a*r.*cos(t); y = b*r.*sin(t);

% ROTATE
pillarGrad = epsilon*(Dposty+Dy)/(Dpostx+Dx);
if itype == 1
  y = y+pillarGrad*(x+L/2-Dx-Dpostx);
else
  y = y+pillarGrad*(x+L/2-Dx/2-maxDist2Left);
end

XwallsExt = [x;y];

% Row shift  
delLat = (Dy+Dposty)*epsilon;


if itype == 2
  % Place horizontally the obstacles in the first column  
  centy1stCol = linspace(-H/2+Dposty/2+Dy/2,H/2-Dposty/2-Dy/2,nrow)';
  centx = linspace(-L/2+Dpostx/2+Dx/2,L/2-Dpostx/2-Dx/2,ncol);

  % this is used for optimization when we have multiple vesicles in sim.
  % x_freeze is prepared for 8 vesicles
  % first 4 vesicles are frozen when they reach the 2nd column
  % last 4 vesicles are frozen when they reach the 4th column
 
  % Place vertically all the obstacles
  centx = centx(ones(nrow,1),:);

  % Shift centers of the obstacles
  delLatVect = [0:ncol-1]*delLat;
  centy = delLatVect(ones(nrow,1),:) + centy1stCol(:,ones(1,ncol));

  % Place the obstacles
  XwallsInt = zeros(2*Nint,nrow*ncol);
  for iwall = 1 : nrow*ncol
    XwallsInt(:,iwall) = [Xpost(1:end/2)+centx(iwall);...
      Xpost(end/2+1:end)+centy(iwall)];
  end
else
  
  % to restrict the flow in type 1 device add ellipses as two more
  % rows at the top and at the bottom
  centy1stCol = linspace(-H/2+Dposty+Dy,H/2-Dposty-Dy,nrow)';
  centx = linspace(-L/2+Dpostx+Dx,L/2-Dpostx-Dx,ncol);


  % Place vertically all the obstacles
  centx = centx(ones(nrow,1),:);

  % Shift centers of the obstacles
  delLatVect = [0:ncol-1]*delLat;
  centy = delLatVect(ones(nrow,1),:) + centy1stCol(:,ones(1,ncol));

  % Place the obstacles
  XwallsInt = zeros(2*Nint,(nrow+2)*ncol);
  for iwall = 1:nrow*ncol
    XwallsInt(:,iwall+ncol) = [Xpost(1:end/2)+centx(iwall);...
      Xpost(end/2+1:end)+centy(iwall)];
  end
  
  centyEllipTop = H/2;
  centyEllipBot = -H/2;
  centxEllip = linspace(-L/2+Dpostx+Dx,L/2-Dpostx-Dx,ncol);
  
  centyTop = delLatVect + centyEllipTop;
  centyBot = delLatVect + centyEllipBot;
  
  % Place ellipses, first bottom then top
  height = 0.8*Dposty/4; % half height
  width = Dpostx/2; % half widht
  tobs = [0:Nint-1]'*2*pi/Nint;
  for j = 1 : ncol
    centyEll = centyBot(j)+Dposty/2-height;
    XwallsInt(:,j) = [width*cos(-tobs)+centxEllip(j); ...
        height*sin(-tobs)+centyEll];
  end
    
  for j = 1 : ncol
    centyEll = centyTop(j)-Dposty/2+height;
    XwallsInt(:,(nrow+1)*ncol+j) = [width*cos(-tobs)+centxEllip(j); ...
        height*sin(-tobs)+centyEll];
  end
      
end
  
end % initConfigDLDRot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XwallsInt,XwallsExt,L,H,radiusRight,radiusUp,...
        Dpostx1,Dposty1,Dpostx2,Dposty2] = ...
        initConfigDLDRot_multiPillar(o,Nint,Next,Xpost1,Xpost2,...
        epsilon,Dx,Dy1,nrow,ncol)       

% Xpost1: triangle, Xpost2: circle
% ASSUMING Xpost1 and Xpost2 CAN FIT IN A RECTANGLE OF THE SAME SIZE

% Smaller DLD device for optimization (nrow x ncol array of a DLD with 
% a period p)  

% if we fit the post in a rectangle, then we put the posts as if we 
% put the corresponding rectangular posts in DLD


xup = interpft(Xpost1(1:end/2),512); 
yup = interpft(Xpost1(end/2+1:end),512);
% distances between top and bottom, left and right of the posts
Dpostx1 = max(xup)-min(xup);
Dposty1 = max(yup)-min(yup);

% Move shape such that the rectangle around the shape is centered
x = Xpost1(1:end/2); y = Xpost1(end/2+1:end);
x2 = x-(max(xup)-Dpostx1/2);
y2 = y-(max(yup)-Dposty1/2);
Xpost1 = [x2;y2];

xup = interpft(Xpost2(1:end/2),512); 
yup = interpft(Xpost2(end/2+1:end),512);
Dpostx2 = max(xup)-min(xup);
Dposty2 = max(yup)-min(yup);
x = Xpost2(1:end/2); y = Xpost2(end/2+1:end);
x2 = x-(max(xup)-Dpostx2/2);
y2 = y-(max(yup)-Dposty2/2);
Xpost2 = [x2;y2];

% radiusLeft=radiusRight, and same for the top-bottom, since the
% post's and rectangle's centers coincides
radiusLeft = Dpostx1/2; radiusRight = radiusLeft;
radiusUp = Dposty1/2; radiusDown = radiusUp;

% Find where the smallest gap size is, place vesicles there and scale the
x1up = interpft(Xpost1(1:end/2),512);
y1up = interpft(Xpost1(end/2+1:end),512);

x2up = interpft(Xpost2(1:end/2),512);
y2up = interpft(Xpost2(end/2+1:end),512);

% look at the points in the horizontal direction and find the largest size
% across the shape in y-direction
maxGsize1 = -1E+4; maxGsize2 = -1E+4;
sort_x1 = sort(x1up); sort_x2 = sort(x2up);
for ij = 1 : numel(x2up)
  
  ids = find(abs(x1up-sort_x1(ij))<=1e-2);
  ys = y1up(ids);
  gapSize_1 = max(ys)-min(ys);
    
  if gapSize_1>=maxGsize1
    maxGsize1 = gapSize_1;
    [~,idx] = max(y1up(ids));
    maxDist2Left1 = x1up(ids(idx))+Dpostx1/2;  
  end
  
  ids = find(abs(x2up-sort_x2(ij))<=1e-2);
  ys = y2up(ids);
  gapSize_2 = max(ys)-min(ys);

  if gapSize_2>=maxGsize2
    maxGsize2 = gapSize_2;
    [~,idx] = max(y2up(ids));
    maxDist2Left2 = x2up(ids(idx))+Dpostx2/2;  
  end
end
    

% suppose same separation in both directions if Dx is not given
if isempty(Dx)
  Dx = Dy1;    
end


H = nrow*(Dy1+Dposty1);
L = ncol*Dx + (ncol+1)/2*Dpostx1 + ((ncol+1)/2-1)*Dpostx2;

% Exterior long tube
a = L/2; b = H/2; order = 20;
t = (0:Next-1)'*2*pi/Next;
r = (cos(t).^order + sin(t).^order).^(-1/order);
x = a*r.*cos(t); y = b*r.*sin(t);

% ROTATE
pillarGrad = epsilon*(Dposty1+Dy1)/(Dpostx1/2+Dpostx2/2+Dx);
y = y+pillarGrad*(x+L/2-Dx/2-maxDist2Left1);


XwallsExt = [x;y];

% Row shift  
delLat = (Dy1+Dposty1)*epsilon;



% Place horizontally the obstacles in the first column  
centy1stCol = linspace(-H/2+Dposty1/2+Dy1/2,H/2-Dposty1/2-Dy1/2,nrow)';
centx = linspace(-L/2+Dpostx1/2+Dx/2,L/2-Dpostx1/2-Dx/2,ncol);

% this is used for optimization when we have multiple vesicles in sim.
% x_freeze is prepared for 8 vesicles
% first 4 vesicles are frozen when they reach the 2nd column
% last 4 vesicles are frozen when they reach the 4th column

% Place vertically all the obstacles
centx = centx(ones(nrow,1),:);

% Shift centers of the obstacles
delLatVect = [0:ncol-1]*delLat;
centy = delLatVect(ones(nrow,1),:) + centy1stCol(:,ones(1,ncol));

% Place the obstacles
XwallsInt = zeros(2*Nint,nrow*ncol);
iwall = 1;
for icol = 1 : ncol
  for irow = 1 : nrow
    if rem(icol,2) ~= 0 
      XwallsInt(:,iwall) = [Xpost1(1:end/2)+centx(irow,icol);...
          Xpost1(end/2+1:end)+centy(irow,icol)];
    else
      XwallsInt(:,iwall) = [Xpost2(1:end/2)+centx(irow,icol);...
          Xpost2(end/2+1:end)+centy(irow,icol)];  
    end
    iwall = iwall + 1; 
  end
end

  
end % initConfigDLDRot_multiPillar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = fillCouetteArea(o,Xref,Xwalls,volFrac,fmm,op)
N = numel(Xref)/2;
x0 = Xref(1:end/2); y0 = Xref(end/2+1:end);
xwalls = Xwalls(1:end/2,:); ywalls = Xwalls(end/2+1:end,:);
% build object for walls
walls = capsules(Xwalls,[],[],0,0,true);
% need radx and rady to decide where to randomly
% sample the geometry.  This way we do not pick 
% a center that is  outside of the solid walls 
% on a regular basis
radx = 1/2*(max(xwalls(:,1)) - min(xwalls(:,1)));
rady = 1/2*(max(ywalls(:,1)) - min(ywalls(:,1)));
% Total area occupied by the physical domain
[~,area,~] = o.geomProp(Xwalls);
areaGeom = sum(area);

% shuffle random number generator
rng('shuffle');

% DONE WITH THE PARAMETERS

% area of a reference vesicle
[~,areaVes,~] = o.geomProp(Xref);

% number of vesicles to initialize
nv = ceil(volFrac*areaGeom/areaVes);
X = zeros(2*N,nv);
XLarge = X;

% counter for the number of successfully placed vesicles
k = 1;
while k <= nv 
   cx = 2*radx*(rand-1/2);
   cy = 2*rady*(rand-1/2);
  
   % center
   phi = 2*pi*rand;
    
   % potential vesicle location
   xpot = cx + x0*cos(phi) - y0*sin(phi);
   ypot = cy + x0*sin(phi) + y0*cos(phi);
  
   xpotLarge = cx + 1.05*(x0*cos(phi) - y0*sin(phi));
   ypotLarge = cy + 1.05*(x0*sin(phi) + y0*cos(phi));

   accept = true; % tentatively accept the vesicle

   % create capsule with the potential vesicle
   vesicle = capsules([xpotLarge;ypotLarge],[],[],[],[],true);
    
   [~,NearV2W] = vesicle.getZone(walls,3);
   [~,icollisionWall] = vesicle.collision(walls,[],NearV2W,fmm,op);
   % at least one of the vesicles's points is outside of one
   % of the solid wall components
   if icollisionWall
     accept = false;
   end
  
   % reject vesicle if it is outside the domain
   if 1.1*sqrt(mean(xpot)^2 + mean(ypot)^2)>=max(xwalls(:,1)) || ...
           0.9*sqrt(mean(xpot)^2 + mean(ypot)^2)<=max(xwalls(:,2)) 
     accept = false;
   end
  
  
   % if vesicle is not outside of physical walls, accept it as
   % a potential new vesicle.  It will be kept as long as it intersects
   % no other vesicles
   if accept 
     X(:,k) = [xpot;ypot];
     XLarge(:,k) = [xpotLarge;ypotLarge]; 
     % create an object with the current configuration of
     % vesicles
     vesicle = capsules(XLarge(:,1:k),[],[],0,0,true);  
      
     % see if vesicles have crossed using collision detection code
     % Can be used with or without the fmm
     NearV2V = vesicle.getZone([],1);
     icollisionVes = vesicle.collision(...
         [],NearV2V,[],fmm,op);
    
     if icollisionVes
       X(:,k) = 0;
       % if they've crossed, reject it
     else
       k = k + 1;
       fprintf('%d Vesicles left to fill the domain to the desired volume fraction\n', nv-k+1)
       % if they haven't crossed, increase the total number of vesicles
     end
   end % if accept

end % while k <= nv  
    


end % fillCouetteArea

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, volFrac] = fillDomain(o,Xref,varargin)
    
options = varargin;

% GET THE PARAMETERS

% Get the domain type (couette, DLD, rectangular like a tube or choke)
if(any(strcmp(options,'domain')))
  ind = find(strcmp(options,'domain'));
  itype = options{ind+1};
else
  itype = [];
end

% Get volume fraction
if(any(strcmp(options,'volFrac')))
  ind = find(strcmp(options,'volFrac'));
  volFrac = options{ind+1};
else
  volFrac = [];
end

% Get number of vesicles
if(any(strcmp(options,'numVes')))
  ind = find(strcmp(options,'numVes'));
  nv = options{ind+1};
else
  nv = [];
end

% Get walls (all of which have the same disc.)
if(any(strcmp(options,'walls')))
  ind = find(strcmp(options,'walls'));
  Xwalls = options{ind+1};
else
  Xwalls = [];
end

% Get interior walls (when there are two sets of walls)
if(any(strcmp(options,'wallsInt')))
  ind = find(strcmp(options,'wallsInt'));
  XwallsInt = options{ind+1};
else
  XwallsInt = [];
end

% Get exterior walls (when there are two sets of walls)
if(any(strcmp(options,'wallsExt')))
  ind = find(strcmp(options,'wallsExt'));
  XwallsExt = options{ind+1};
else
  XwallsExt = [];
end

% Get the tstep object
if(any(strcmp(options,'tstepObj')))
  ind = find(strcmp(options,'tstepObj'));
  tt = options{ind+1};
else
  tt = [];
end

% Get xrange and yrange
if(any(strcmp(options,'xrange')))
  ind = find(strcmp(options,'xrange'));
  xrange = options{ind+1};
else
  xrange = [];
end
if(any(strcmp(options,'yrange')))
  ind = find(strcmp(options,'yrange'));
  yrange = options{ind+1};
else
  yrange = [];
end

% Get the coordinates of left bottom point
if(any(strcmp(options,'leftBotX')))
  ind = find(strcmp(options,'leftBotX'));
  leftBotX = options{ind+1};
else
  leftBotX = [];
end
if(any(strcmp(options,'leftBotY')))
  ind = find(strcmp(options,'leftBotY'));
  leftBotY = options{ind+1};
else
  leftBotY = [];
end

% Get the length and height of the domain to fill
if(any(strcmp(options,'length')))
  ind = find(strcmp(options,'length'));
  L = options{ind+1};
else
  L = [];
end
if(any(strcmp(options,'height')))
  ind = find(strcmp(options,'height'));
  H = options{ind+1};
else
  H = [];
end

% Get the DLD domain information
if(any(strcmp(options,'Dpost')))
  ind = find(strcmp(options,'Dpost'));
  Dpost = options{ind+1};
else
  Dpost = [];
end

if(any(strcmp(options,'Dx')))
  ind = find(strcmp(options,'Dx'));
  Dx = options{ind+1};
else
  Dx = [];
end

if(any(strcmp(options,'Dy')))
  ind = find(strcmp(options,'Dy'));
  Dy = options{ind+1};
else
  Dy = [];
end

if(any(strcmp(options,'epsilon')))
  ind = find(strcmp(options,'epsilon'));
  epsilon = options{ind+1};
else
  epsilon = [];
end
if(any(strcmp(options,'ncol')))
  ind = find(strcmp(options,'ncol'));
  ncol = options{ind+1};
else
  ncol = [];
end
if(any(strcmp(options,'xfreeze')))
  ind = find(strcmp(options,'xfreeze'));
  xfreeze = options{ind+1};
else
  xfreeze = [];
end

% DONE WITH THE PARAMETERS

% reference vesicle
[x0,y0] = o.getXY(Xref);
N = numel(x0);

% area of a reference vesicle
[~,areaVes,~] = o.geomProp(Xref);

if strcmp(itype,'couette')
  % solid walls      
  [xwalls,ywalls] = o.getXY(Xwalls);
  nvbd = 2;
  % Total area occupied by the physical domain
  [~,area,~] = o.geomProp(Xwalls);
  areaGeom = sum(area);
  
  if isempty(nv)
    % total number of vesicles to achieve volume fraction
    nv = ceil(volFrac*areaGeom/areaVes);
  else
    volFrac = nv*areaVes/areaGeom;  
  end

  % build object for walls
  walls = capsules(Xwalls,[],[],0,0,true);
  % need radx and rady to decide where to randomly
  % sample the geometry.  This way we do not pick 
  % a center that is  outside of the solid walls 
  % on a regular basis
  radx = 1/2*(max(xwalls(:,1)) - min(xwalls(:,1)));
  rady = 1/2*(max(ywalls(:,1)) - min(ywalls(:,1)));
  
  X = zeros(2*N,nv);
  % counter for the number of successfully placed vesicles
  k = 1;

  while k <= nv 
    cx = 2*radx*(rand-1/2);
    cy = 2*rady*(rand-1/2);
  
    % center
    phi = 2*pi*rand;
    
    % potential vesicle location
    xpot = cx + x0*cos(phi) - y0*sin(phi);
    ypot = cy + x0*sin(phi) + y0*cos(phi);
  

    accept = true; % tentatively accept the vesicle

    % create capsule with the potential vesicle
    vesicle = capsules([xpot;ypot],[],[],[],[],true);
    
    [~,NearV2W] = vesicle.getZone(walls,3);
    [~,icollisionWall] = vesicle.collision(walls,[],NearV2W,tt.fmm,tt.op);
    % at least one of the vesicles's points is outside of one
    % of the solid wall components
    if icollisionWall
      accept = false;
    end
  
    % reject vesicle if it is outside the domain
    if 1.1*sqrt(mean(xpot)^2 + mean(ypot)^2)>=max(xwalls(:,1)) || ...
            0.9*sqrt(mean(xpot)^2 + mean(ypot)^2)<=max(xwalls(:,2)) 
      accept = false;
    end
  
  
    % if vesicle is not outside of physical walls, accept it as
    % a potential new vesicle.  It will be kept as long as it intersects
    % no other vesicles
    if accept 
      X(:,k) = [xpot;ypot];
      
      % create an object with the current configuration of
      % vesicles
      vesicle = capsules(X(:,1:k),[],[],0,0,true);
      
      % see if vesicles have crossed using collision detection code
      % Can be used with or without the fmm
      NearV2V = vesicle.getZone([],1);
      icollisionVes = vesicle.collision(...
          [],NearV2V,[],tt.fmm,tt.op);
    
      if icollisionVes
        X(:,k) = 0;
        % if they've crossed, reject it
      else
        k = k + 1;
        fprintf('%d Vesicles left to fill the domain to the desired volume fraction\n', nv-k+1)
        % if they haven't crossed, increase the total number of vesicles
      end
    end % if accept

  end % while k <= nv  
    
elseif strcmp(itype,'rectangular')
  % Total area occupied by the physical domain  
  areaGeom = L*H;
  
  if isempty(nv)
    % total number of vesicles to achieve volume fraction
    nv = ceil(volFrac*areaGeom/areaVes);
  else
    volFrac = nv*areaVes/areaGeom;  
  end
  
  % maximum and minimum x and y coordinates for which vesicles can reside
  leftX = leftBotX; rightX = leftX+L;
  botY = leftBotY; topY = botY+H;

  X = zeros(2*N,nv);
  XLarge = X;
  
  % counter for the number of successfully placed vesicles
  k = 1;

  while k <= nv 
    % center  
    cx = leftX + L*rand;
    cy = botY + H*rand;
  
    phi = 2*pi*rand;
    
    xpot = cx + x0*cos(phi) - y0*sin(phi);
    ypot = cy + x0*sin(phi) + y0*cos(phi);
  
    % potential vesicle location (make it larger to avoid ves-ves
    % collision)
    xpotLarge = cx + 1.1*(x0*cos(phi) - y0*sin(phi));
    ypotLarge = cy + 1.1*(x0*sin(phi) + y0*cos(phi));
  
    accept = true; % tentatively accept the vesicle
    
    % reject vesicle if it is outside the domain
    if 1.1*mean(xpot)>=rightX || 0.9*mean(xpot)<=leftX 
      accept = false;
    end
    if 1.1*mean(ypot)>=topY || 0.9*mean(xpot)<=botY
      accept = false;
    end

    if accept 
      % if vesicle is not outside of physical walls, accept it as
      % a potential new vesicle.  It will be kept as long as it intersects
      % no other vesicles  
      X(:,k) = [xpot;ypot];
      XLarge(:,k) = [xpotLarge;ypotLarge];
      
      % create an object with the current configuration of
      % vesicles
      vesicle = capsules(XLarge(:,1:k),[],[],0,0,true);
      
      % see if vesicles have crossed using collision detection code
      % Can be used with or without the fmm
      NearV2V = vesicle.getZone([],1);
      icollisionVes = vesicle.collision(...
          [],NearV2V,[],tt.fmm,tt.op);
    
      if icollisionVes
        X(:,k) = 0;
        XLarge(:,k) = 0;
        % if they've crossed, reject it
      else
        k = k + 1;
        fprintf('%d Vesicles left to fill the domain to the desired volume fraction\n', nv-k+1)
        % if they haven't crossed, increase the total number of vesicles
      end % if icollisionVes
    end % if accept

  end % while k <= nv
  
elseif strcmp(itype,'DLDEntrance')
  om = tt.om;
  % initializes vesicles only in the entrance  
    
  % build object for walls      
  wallsExt = capsules(XwallsExt,[],[],0,0,true);
  wallsInt = capsules(XwallsInt,[],[],0,0,true);

  % hematocrit, volume fraction of vesicles
  areaGeom = abs((yrange(2)-yrange(1))*(xrange(2)-xrange(1)));
  volFrac = areaVes*nv/areaGeom;

  X = zeros(2*N,nv);
  XLarge = X;

  % Counter
  k = 1;
  while k <= nv 
 
    % randomly choose vesicle's center
    cx = xrange(1) + (xrange(2)-xrange(1))*rand;
    cy = yrange(1) + (yrange(2)-yrange(1))*rand;
  
  
    phi = -pi/6 + pi/3*rand;
    % potential vesicle (rotated a bit)
    xpot = cx + x0*cos(phi) - y0*sin(phi);
    ypot = cy + x0*sin(phi) + y0*cos(phi);
  
    xpotLarge = cx + 1.1*(x0*cos(phi) - y0*sin(phi));
    ypotLarge = cy + 1.1*(x0*sin(phi) + y0*cos(phi));

    accept = true; % tentatively accept the vesicle

    vesicle = capsules([xpotLarge;ypotLarge],[],[],[],[],true);
  
    fprintf('Obtaining vesicle-wall getZone...\n');
    [~,NearV2WInt] = vesicle.getZone(wallsInt,3);
    [~,NearV2WExt] = vesicle.getZone(wallsExt,3);
    fprintf('DONE.\n')
  
  
    % create capsule with the potential vesicle
    [~,icollisionWallExt] = vesicle.collision(wallsExt,[],NearV2WExt,...
        tt.fmm,tt.op);
    [~,icollisionWallInt] = vesicle.collision(wallsInt,[],NearV2WInt,...
        tt.fmm,tt.op);
    if icollisionWallExt || icollisionWallInt
      accept = false;
      
      message = ['Vesicle crossed the outer wall.'];
      om.writeMessage(message,'%s\n')
      % at least one of the vesicles's points is outside of one
      % of the solid wall components
    end
  
    if ~(cx <= xrange(2) && cx>= xrange(1) && cy <= yrange(2) ...
            && cy >= yrange(1))
        accept = false;
      
        message = ['Vesicle was outside the domain.'];
        om.writeMessage(message,'%s\n')
    end
    % reject vesicle if it is outside the domain
  
    if accept 
      % if vesicle is not outside of physical walls, accept it as
      % a potential new vesicle.  It will be kept as long as it intersects
      % no other vesicles  
      X(:,k) = [xpot;ypot];
      XLarge(:,k) = [xpotLarge;ypotLarge];
    
      % create an object with the current configuration of vesicles
      vesicle = capsules(XLarge(:,1:k),[],[],0,0,true);  
    
    
      % see if vesicles have crossed using collision detection code
      % Can be used with or without the fmm    
      fprintf('Obtaining vesicle-vesicle getZone...\n');
      NearV2V = vesicle.getZone([],1);
      icollisionVes = vesicle.collision(...
          [],NearV2V,[],tt.fmm,tt.op);
      fprintf('DONE.\n');
    
    
      if icollisionVes
        X(:,k) = 0;
        XLarge(:,k) = 0;
        message = ['Vesicles crossed.'];
        om.writeMessage(message,'%s\n')
        % if they've crossed, reject it
      else
        k = k + 1;
        message = [num2str(nv-k+1,'%d') ' vesicles left to fill the domain'];
        om.writeMessage(message,'%s\n')
        % if they haven't crossed, increase the total number of vesicles
      end
    end
  end % while k <= nv
  
elseif strcmp(itype,'DLDWhole')
  % fill in the lanes in the middle, one above, and one below for a given 
  % volume fraction
  
  om = tt.om;
  
  % Left and right x-coordinate of the lanes (at the very left point 
  % of the first pillar) + some buffer
  xLeft = min(XwallsInt(1:end/2,1));
  xRight = xfreeze;
  
  % range for vesicles
  xrange = [xLeft;xRight]-0.05*abs([xLeft;xRight]);
  
  % build objects for walls
  wallsExt = capsules(XwallsExt,[],[],0,0,true);
  wallsInt = capsules(XwallsInt,[],[],0,0,true);
  
  % fill in 3 lanes (each is a parallelogram)
  areaGeom = 3*(xrange(2)-xrange(1))*Dy;
  
  % number of vesicles to fill 3 lanes
  nv = ceil(volFrac*areaGeom/areaVes); 
    
  X = zeros(2*N,nv);
  XLarge = X;
  
  % Region info (1: above, 2: middle, 3: below)
  regionyLB(1) = Dy/2+Dpost;
  regionyLB(2) = -Dy/2;
  regionyLB(3) = -1.5*Dy-Dpost;
  nvInRegs = zeros(3,1);
  
  % counter
  k = 1;
  
  while k <= nv
    
    % randomly pick x-location for a center of vesicle
    cx = xrange(1) + (xrange(2)-xrange(1))*rand;
    
    % we need to pick cy now, but decide on the lane, so
    % find the region which has the minimum number of vesicles 
    [~,regID] = min(nvInRegs);
    
    yBot = (cx-xLeft)*epsilon+regionyLB(regID);
    yTop = yBot+Dy; 
    yrange = [yBot+0.05*abs(yBot); yTop-0.05*abs(yTop)];
    
    cy = yrange(1) + (yrange(2)-yrange(1))*rand;
    
    phi = pi/8*rand;
    % potential vesicle
    
    xpot = cx + x0*cos(phi) + y0*sin(phi);
    ypot = cy - x0*sin(phi) + y0*cos(phi);
  
    xpotLarge = cx + 1.1*(x0*cos(phi) + y0*sin(phi));
    ypotLarge = cy - 1.1*(x0*sin(phi) - y0*cos(phi));

    accept = true; % tentatively accept the vesicle

    % create capsule with the potential vesicle
    vesicle = capsules([xpotLarge;ypotLarge],[],[],[],[],true);
  
    fprintf('Obtaining vesicle-wall getZone...\n');
    [~,NearV2WInt] = vesicle.getZone(wallsInt,3);
    [~,NearV2WExt] = vesicle.getZone(wallsExt,3);
    fprintf('DONE.\n')
    
     
    [~,icollisionWallExt] = vesicle.collision(wallsExt,[],NearV2WExt,...
        tt.fmm,tt.op);
    [~,icollisionWallInt] = vesicle.collision(wallsInt,[],NearV2WInt,...
        tt.fmm,tt.op);
    if icollisionWallExt || icollisionWallInt
      accept = false;
      
      message = ['Vesicle crossed the outer wall.'];
      om.writeMessage(message,'%s\n')
      % at least one of the vesicles's points is outside of one
      % of the solid wall components
    end
    
    if ~(cx <= xrange(2) && cx>= xrange(1) && cy <= yrange(2) ...
            && cy >= yrange(1))
        accept = false;
      
        message = ['Vesicle was outside the domain.'];
        om.writeMessage(message,'%s\n')
    end
    % reject vesicle if it is outside the domain
    
    if accept 
      % if vesicle is not outside of physical walls, accept it as
      % a potential new vesicle.  It will be kept as long as it intersects
      % no other vesicles  
      X(:,k) = [xpot;ypot];
      XLarge(:,k) = [xpotLarge;ypotLarge];
    
      % create an object with the current configuration of vesicles
      vesicle = capsules(XLarge(:,1:k),[],[],0,0,true);  
    
    
      % see if vesicles have crossed using collision detection code
      % Can be used with or without the fmm    
      fprintf('Obtaining vesicle-vesicle getZone...\n');
      NearV2V = vesicle.getZone([],1);
      icollisionVes = vesicle.collision(...
          [],NearV2V,[],tt.fmm,tt.op);
      fprintf('DONE.\n');
    
    
      if icollisionVes
        X(:,k) = 0;
        XLarge(:,k) = 0;
        message = ['Vesicles crossed.'];
        om.writeMessage(message,'%s\n')
        % if they've crossed, reject it
      else
        nvInRegs(regID) = nvInRegs(regID) + 1;  
        k = k + 1;
        message = [num2str(nv-k+1,'%d') ' vesicles left to fill the domain'];
        om.writeMessage(message,'%s\n')
        % if they haven't crossed, increase the total number of vesicles
      end
    end % if accept
    
  end % while k <= 3*nv 
      
      
end % end if

end % fillDomain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = seedVesicles(o,seedRate,Xsamp,Xact,xrange,yrange,tt)

k = 1;
% Counter
nFails = 0;
% number of failed attempts to initialize vesicles
while k<=seedRate && nFails<=50
 
  kSample = 1 + floor((size(Xsamp,2)-1)*rand);
  
  x00 = Xsamp(1:end/2,kSample)-mean(Xsamp(1:end/2,kSample));
  y00 = Xsamp(end/2+1:end,kSample)-mean(Xsamp(end/2+1:end,kSample));
  
  % rotate to 0 angle
  IA = o.getIncAngle([x00;y00]);
  
  x0 = x00*cos(-IA) - y00*sin(-IA);
  y0 = x00*sin(-IA) + y00*cos(-IA);
  % choose a relaxed sample vesicle and move it to (0,0)
  
  cx = xrange(1) + (xrange(2)-xrange(1))*rand;
  cy = yrange(1) + (yrange(2)-yrange(1))*rand;
  % randomly choose vesicle's center
  
  % randomly rotate it
  phi = -pi/6 + pi/3*rand;
  xpot = cx + x0*cos(phi) - y0*sin(phi);
  ypot = cy + x0*sin(phi) + y0*cos(phi);
  % potential vesicle 
  xpotLarge = cx + 1.2*(x0*cos(phi) - y0*sin(phi));
  ypotLarge = cy + 1.2*(x0*sin(phi) + y0*cos(phi));
    
  accept = true;
  if ~(cx <= xrange(2) && cx>= xrange(1) && cy <= yrange(2) ...
          && cy >= yrange(1))
      accept = false;
  end
  % reject vesicle if it is outside the domain (not very common)
  
  if accept
    Xnew(:,k) = [xpot;ypot];
    XLarge(:,k) = [xpotLarge;ypotLarge];
    % if vesicle is not outside of physical walls, accept it as
    % a potential new vesicle.  It will be kept as long as it intersects
    % no other vesicles

    vesicle = capsules([XLarge(:,1:k) Xact],[],[],0,0,true);  
    % create an object with the current configuration of
    % vesicles
    
    NearV2V = vesicle.getZone([],1);
    icollisionVes = vesicle.collision(...
      [],NearV2V,[],tt.fmm,tt.op);
  
    
    % see if vesicles have crossed using collision detection code
    % Can be used with or without the fmm
    if icollisionVes
      Xnew(:,k) = 0;
      XLarge(:,k) = 0;

      nFails = nFails+1;
      % if it tries more than 50 times to initialize the vesicle and 
      % fails then do not continue initializing. Probably not enough 
      % space to initialize more vesicles.

    else
      k = k + 1;
      nFails = 0;
    end
  end

end % while k <= nv

if norm(Xnew(:,end))==0
  Xnew = Xnew(:,1:end-1);  
end    
% remove the last vesicle if it is all zeros, that's the case 
% where it failed to generate a vesicle
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = addAndRemove(o,X,walls,near,fmm)
% X = addAndRemove(X,walls,near,fmm) is used to fudge in periodicity.
% When a vesicle passes a certain threshold (hard-coded into this
% function), it randomly places it at the start of the computational
% domain with the same shape.  If a collision with another vesicle
% results, it tries again.  If too many attempts are made, it simply
% continues with its current position and tries again at the next
% iteration

near = false;
[x,y] = o.getXY(X);
[N,nv] = size(x);

ind = find(max(x) > 8);
if numel(ind) > 0
  [~,s] = max(max(x(:,ind)));
  % find the index of the point farthest to the right
  x = x(:,[1:ind(s)-1 (ind(s)+1:nv)]);
  y = y(:,[1:ind(s)-1 (ind(s)+1:nv)]);
  Xminus = o.setXY(x,y);
  % all points except the one that is farthest to the right

  icollisionVes = true;
  maxTrys = 10;
  % maximum number of attempts before attempt to build in peridocity is
  % postponed
  ntrys = 0;
  while (icollisionVes && ntrys < maxTrys)
    ntrys = ntrys + 1;
%    center = [-8.5;0.2*rand(1)];
    center = [-8;0.7*2*(rand(1)-0.5)];
    % random center
    Xadd = X(:,ind(s));
    Xadd(1:end/2) = Xadd(1:end/2) - mean(Xadd(1:end/2));
    Xadd(1+end/2:end) = Xadd(1+end/2:end) - mean(Xadd(1+end/2:end));
    Xadd(1:end/2) = Xadd(1:end/2) + center(1);
    Xadd(1+end/2:end) = Xadd(1+end/2:end) + center(2);
    if max(abs(Xadd(1+end/2:end,:)) > 9.5e-1)
      icollisionWall = true;
    else
      icollisionWall = false;
    end
    % solid wall is at about 0.95 at the point -8.5
    if ~icollisionWall
      % parameterization of new vesicle
      vesicle = capsules([Xminus Xadd],[],[],0,0,o.antiAlias);
      % class corresponding to the new configuration
      if near
        [NearV2V,NearV2W] = vesicle.getZone(walls,3);
      else
        NearV2V = [];
        NearV2W = [];
      end
      icollisionVes = o.collision(vesicle,walls,...
          NearV2V,NearV2W,fmm,near);
      % check if the vesicles have crossed
    end
  end

  if ntrys < maxTrys
    X = vesicle.X;
  end
  % if it's tried to many times, simply move on and try again at the
  % next iteration
end

end % addAndRemove
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,u,Xstore,sigStore,uStore,rem_ids,genSoFar_ids,prams,...
        om,nInstant] = freezeAndStream(o,om,tt,Xstore,sigStore,uStore,X,...
        sigma,u,prams,XwallsExt,rem_ids,genSoFar_ids,istreaming,...
        allVesViscConts,time,nInstant)
%  [X,sigma,u,Xstore,sigStore,uStore,rem_ids,genSoFar_ids,prams,...
%         om] = freezeAndStream(o,om,tt,Xstore,sigStore,uStore,X,...
%         sigma,u,prams,XwallsExt,rem_ids,genSoFar_ids,istreaming)    
% decides which vesicles are going to be frozen and streams new vesicles.
% This is written only for DLD examples. Freezing for other examples
% requires another subroutine.

% ids of the vesicles to be frozen and to be remained
will_freeze=[]; will_remain=[];

% length of the exterior wall
Lext = max(XwallsExt(1:end/2))-min(XwallsExt(1:end/2));
delLambda = (prams.Dposty+prams.Dy)*prams.epsilon;

% yZigZag is the lateral position of the center of the post above which the
% vesicle is initialized. That center moves up by delLambda after every
% column. 

% Loop over all the active vesicles
for iv = 1:size(X,2)
    
  % Which column the vesicle has just left behind
  atCol = floor((mean(X(1:end/2,iv))+Lext/2)/(prams.Dx+prams.Dpostx));
    
  % if the centx of a vesicle reaches x_freeze or it is 
  % close to the top and bottom of the DLD device, then freeze it 
%   if numel(prams.xfreeze) == 1
  if prams.freezeIfZigZag
    if mean(X(1:end/2,iv)) >= prams.xfreeze || ...
            mean(X(end/2+1:end,iv))<=prams.yZigZag(iv)+(atCol-1)*delLambda

      % tag the ones we need to freeze
      will_freeze = [will_freeze;iv];

      message = ['Vesicle ' num2str(rem_ids(iv),'%d') ' is frozen'];
      om.writeMessage(message,'%s\n')
    else
      % tag the ones we need to keep
      will_remain = [will_remain;iv];
    end % if mean
  else
    if mean(X(1:end/2,iv)) >= prams.xfreeze

      % tag the ones we need to freeze
      will_freeze = [will_freeze;iv];

      message = ['Vesicle ' num2str(rem_ids(iv),'%d') ' is frozen'];
      om.writeMessage(message,'%s\n')
    else
      % tag the ones we need to keep
      will_remain = [will_remain;iv];
    end % if mean
  end % if prams.freezeIfZigZag
    
%   else % then different vesicles have different xfreeze points
%     if mean(X(1:end/2,iv)) > prams.xfreeze(iv) || ...
%             max(X(1:end/2,iv)) > 0.98*max(XwallsExt(1:end/2))
%         
%       % tag the ones we need to freeze
%       will_freeze(idx) = iv;
%       idx = idx+1;
% 
%       message = ['Vesicle ' num2str(rem_ids(iv),'%d') ' is frozen'];
%       om.writeMessage(message,'%s\n')
%     else
%       % tag the ones we need to keep
%       will_remain(idx2) = iv;
%       idx2 = idx2+1;        
%     end % if mean  
%       
%       
%   end % if numel(prams.xfreeze)==1
end % for iv

% save the frozen ones and remained ones
Xstore(:,rem_ids) = X;
sigStore(:,rem_ids) = sigma;
uStore(:,rem_ids) = u;

if~isempty(will_freeze)
  % then there are vesicles to be frozen  
  if istreaming
    % then save them, b/c we may want to plot them  
    Xfrozen = X(:,will_freeze);
    frozenIds = rem_ids(will_freeze);
  end

  X = X(:,will_remain);
  sigma = sigma(:,will_remain);
  u = u(:,will_remain);
  om.area = om.area(will_remain);
  om.length = om.length(will_remain);
  rem_ids = rem_ids(will_remain);
  prams.viscCont = allVesViscConts(rem_ids);
  if prams.freezeIfZigZag
    prams.yZigZag = prams.yZigZag(will_remain);
  end
  prams.nv = numel(rem_ids);
else
  Xfrozen = [];
  frozenIds = [];
end

if istreaming && numel(genSoFar_ids)<prams.totnv
  % this sets the streaming rate by extending the area where we check if
  % there is any vesicle, and if not, then seed. 
  checkArea = 1/prams.streamRate;
  
  % if all the vesicles are not seeded, keep seeding  
  
  % do not seed more than the predetermined value
  nSeed = prams.nSeed; % number of cells to stream 
  seedRate = min(nSeed,prams.totnv-numel(genSoFar_ids));  

  % check if there are vesicles where we seed them, if yes, do not
  % add
  anyVesInlet = false;
  for iv = 1:size(X,2)
    if min(X(1:end/2,iv))<=prams.xrange(1)+checkArea*0.5*(prams.xrange(2)-...
      prams.xrange(1))
      anyVesInlet = true;
    end
  end
    
  % if there is no vesicle in the inlet, add vesicles
  if ~anyVesInlet
    Xnew = o.seedVesicles(seedRate,...
        Xstore(:,1:max(genSoFar_ids)),X,...
        [prams.xrange(1) prams.xrange(1)+0.5*...
        (prams.xrange(2)-prams.xrange(1))],prams.yrange,tt);

    % number of recently seeded vesicles
    nSeeded = size(Xnew,2);

    % keep the ids of vesicles generated so far
    genSoFar_ids = [genSoFar_ids;genSoFar_ids(end)+[1:nSeeded]'];

    % add them to the remaining ids as wall
    rem_ids = [rem_ids;genSoFar_ids(end-nSeeded+1:end)];

    % save the new ones to the matrices keeping info about all vesicles
    Xstore(:,genSoFar_ids(end-nSeeded+1:end)) = Xnew;
    sigStore(:,genSoFar_ids(end-nSeeded+1:end)) = zeros(prams.N,nSeeded);
    uStore(:,genSoFar_ids(end-nSeeded+1:end)) = zeros(2*prams.N,nSeeded);

    prams.viscCont = allVesViscConts(rem_ids);
    X = [X Xnew];
    sigma = [sigma zeros(prams.N,nSeeded)];
    u = [u zeros(2*prams.N,nSeeded)];

    [~,areaNew,lengthNew] = o.geomProp(Xnew);
    om.area = [om.area areaNew];
    om.length = [om.length lengthNew];
    om.areaAll = [om.areaAll areaNew];
    om.lengthAll = [om.lengthAll lengthNew];

    prams.nv = prams.nv + nSeeded;
  end % ~anyVesInlet

end % istreaming && numel(genSoFar_ids)<prams.totnv

if istreaming
  % now save data, what we need is global ids of active vesicles (including 
  % recently streamed ones), current
  % time, their positions, and also the global ids and positions of the
  % vesicles frozen in this time step
  nInstant = nInstant + 1;
  fileName = [prams.folderName prams.runName '_DataInst_' num2str(nInstant) '.mat'];
  save(fileName,'X','rem_ids','time','Xfrozen','frozenIds');       
end
end % freezeAndStream

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,u,Xstore,sigStore,uStore,rem_ids,nPutBacks,prams,...
        om] = freezeAndPutBack(o,om,Xstore,sigStore,uStore,X,...
        sigma,u,prams,rem_ids,allVesViscConts,nPutBacks,nPutBackMax)
    
% if vesicle is beyond xfreeze, then move it back to
% pillarGrad = prams.epsilon*(prams.Dposty+prams.Dy)/(prams.Dpostx+prams.Dx);
pillarGrad = prams.epsilon*(prams.Dposty1+prams.Dy1)/(prams.Dpostx1/2+prams.Dpostx2/2+prams.Dx);
for k = 1 : size(X,2)
  centxn = mean(interpft(X(1:end/2,k),256));  
  if centxn>=prams.xfreeze
    centyn = mean(interpft(X(end/2+1:end,k),256));
    centxNew = centxn-(prams.Dx+prams.Dpostx1/2+prams.Dpostx2/2);
    centyNew = centyn-pillarGrad*(prams.Dx+prams.Dpostx1/2+prams.Dpostx2/2);
    X(:,k) = [X(1:end/2,k)-centxn+centxNew;X(end/2+1:end,k)-centyn+centyNew];
    nPutBacks(k) = nPutBacks(k)+1;
  end
end
      
idx = 1; idx2 = 1;

% ids of the vesicles to be frozen and to be remained
will_freeze=[]; will_remain=[];

% Loop over all the active vesicles
for iv = 1:size(X,2)
    
  
  if nPutBackMax > ceil(prams.periodN)+1 % then do not freeze zig-zagging ones
    if mean(X(end/2+1:end,iv))<=prams.yZigZag(iv)
      % if zig-zags move it up again  
      X(end/2+1:end,iv) = X(end/2+1:end,iv)+prams.Dposty+prams.Dy;
    end
      
    if nPutBacks(iv)>=nPutBackMax
      % tag the ones we need to freeze
      will_freeze(idx) = iv;
      idx = idx+1;

      message = ['Vesicle ' num2str(rem_ids(iv),'%d') ' is frozen'];
      om.writeMessage(message,'%s\n')    
    else
      % tag the ones we need to keep
      will_remain(idx2) = iv;
      idx2 = idx2+1;        
    end
      
  else % we freeze the zig-zagging ones
    % Freeze if nPutBacks exceeds nPutBackMax or vesicle has just zig-zagged  
    if  nPutBacks(iv)>=nPutBackMax || mean(X(end/2+1:end,iv))<=prams.yZigZag(iv)

      % tag the ones we need to freeze
      will_freeze(idx) = iv;
      idx = idx+1;

      message = ['Vesicle ' num2str(rem_ids(iv),'%d') ' is frozen'];
      om.writeMessage(message,'%s\n')
    else
      % tag the ones we need to keep
      will_remain(idx2) = iv;
      idx2 = idx2+1;        
    end  % if mean    
      
  end % nPutBackMax > prams.periodN+1
  
end % for iv

% save the frozen ones and remained ones
Xstore(:,rem_ids) = X;
sigStore(:,rem_ids) = sigma;
uStore(:,rem_ids) = u;

if~isempty(will_freeze)
  % then there are vesicles to be frozen  
  X = X(:,will_remain);
  sigma = sigma(:,will_remain);
  u = u(:,will_remain);
  om.area = om.area(will_remain);
  om.length = om.length(will_remain);
  rem_ids = rem_ids(will_remain);
  prams.viscCont = allVesViscConts(rem_ids);
  nPutBacks = nPutBacks(will_remain);
  prams.yZigZag = prams.yZigZag(will_remain);
  prams.nv = numel(rem_ids);
else
  Xfrozen = [];
  frozenIds = [];
end


end % freezeAndPutBack
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = checkCollisionForLRCA(o,Xcorr,X,walls,fmm,op,om,uprate)
% X = checkCollisionForLRCA(o,Xcorr,X,walls,fmm,op,om)
% checks if LRCA (reparam. and correction lead to collision) if so then we
% undo the correction.

% Build vesicle structure for corrected shape
% antialiasing is on
if uprate == 1
  vesicle = capsules(Xcorr,[],[],[],[],1);
else
  Nup = uprate*size(Xcorr,1);  
  Xup = [interpft(X(1:end/2,:),Nup);interpft(Xcorr(end/2+1:end,:),Nup)];  
  vesicle = capsules(Xup,[],[],[],[],1);
end

% Get the near structure
if ~isempty(walls)
  [NearV2V,NearV2W] = vesicle.getZone(walls,3);
else
  [NearV2V,~] = vesicle.getZone([],1);
  NearV2W = [];
end

% check if there is a collision
[icollisionVes,icollisionWall] = vesicle.collision(walls,NearV2V,...
  NearV2W,fmm,op);

% if collision, do not correct, if not correct
if icollisionVes || icollisionWall
  message = ['LRCA RESULTS IN CROSSING, SO UNDO CORRECTION'];
  om.writeMessage(message,'%s\n')
else
  X = Xcorr;
end

end % checkCollisionForLRCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = fixCollisionsZinc(o,X,Xwalls)

idebug = false;
N = size(X,1)/2; nv = size(X,2); nvbd = size(Xwalls,2); Nbd = size(Xwalls,1)/2;
confined = ~isempty(Xwalls);

% near-zone defined as tol*spacing
[~,~,len] = o.geomProp(X);
beta = 1.1; % small buffer
htol = max(len)/N;
[~,~,lenW] = o.geomProp(Xwalls);
htolWall = min(lenW/Nbd);
% htolWall = htol;

% build vesicle class
vesicle = capsules(X,[],[],[],[],1);
vesicle.setUpRate();

% get near structure
NearV2V = vesicle.getZone([],1);
if confined
  walls = capsules(Xwalls,[],[],[],[],1);
  walls.setUpRate();  
  [~,NearW2V] = walls.getZone(vesicle,2);
  zoneW2V = NearW2V.zone;
  nearestW2V = NearW2V.nearest;
  icpW2V = NearW2V.icp;
end
zoneV2V = NearV2V.zone;
nearestV2V = NearV2V.nearest;
icpV2V = NearV2V.icp;

% Smooth function to move points
interpOrder = 6; % include 1/8 if N points to move (std = 1/5 of order)
% So, the region of contact is L/8
xfilt = linspace(-interpOrder/2,interpOrder/2,interpOrder+1);
yfilt = exp(-xfilt'.^2/(2*1.5));
% normalize the function so that the middle point is 1
yfilt = yfilt/yfilt(interpOrder/2+1); 

% Go over each vesicle, check its interaction with other vesicles
Xnew = X;
for k1 = 1 : nv
  
  % vesicle-vesicle interaction first
  K = [(1:k1-1) (k1+1:nv)];
  X0 = Xnew(:,k1); % save the original shape
  for k2 = K
    J = find(zoneV2V{k2}(:,k1) == 1);
    % set of points on vesicle k2 close to vesicle k1   
    if (numel(J) ~= 0)
      for i = 1:numel(J)
        % distance vector from a point on k1 to its projection on k2
        distx = (Xnew(J(i),k1) - nearestV2V{k2}(J(i),k1));
        disty = (Xnew(J(i)+N,k1) - nearestV2V{k2}(J(i)+N,k1));
        dist = sqrt(distx^2+disty^2);
        if dist <= htol
          if idebug
          figure(1);clf;hold on;
          title('Detected')
          plot(Xnew(1:end/2,:),Xnew(end/2+1:end,:),'r','linewidth',2)
          plot(Xwalls(1:end/2,:),Xwalls(end/2+1:end,:),'k','linewidth',2)
          plot(Xnew(J(i),k1),Xnew(J(i)+N,k1),'bo','markerfacecolor','b','markersize',8);
          plot(Xnew(1:N,k2),Xnew(N+1:2*N,k2),'ko','markerfacecolor','k','markersize',6);
          plot(nearestx,nearesty,'go','markerfacecolor','g','markersize',8)
          axis equal
          pause
          end
        
          % normal direction between the point on k2 and the projection of k2
          % on vesicle k1
          nx = distx/dist; ny = disty/dist;
          % Move only 1/2 of htol*beta(buffer zone)
          % solve (distx/2+kmove*nx)^2+(disty/2+kmove*ny)^2 =
          % (beta*htol/2)^2
          quadB = (distx*nx+disty*ny); quadC = dist^2/4-beta^2*htol^2/4;
          kmove = max(roots([1 quadB quadC]));
          if kmove <= 0 || ~isreal(kmove)
            kmove = beta*htol/2-dist/2; 
          end
        
          % move the points in the normal direction by kmove
          if 1 % move several points using a smooth function
            % find nearby points' idcs (the middle point is J(i))
            idcs = mod((J(i)-interpOrder/2:J(i)+interpOrder/2)'-1,N)+1;
            % distribute the motion, make sure pth point move kmove*nx
            Xnew(idcs,k1) = Xnew(idcs,k1) + yfilt*kmove*nx;
            % for y-component
            % distribute the motion, make sure pth point move kmove*nx
            Xnew(idcs+N,k1) = Xnew(idcs+N,k1) + yfilt*kmove*ny;
          else % move only one point
            Xnew(J(i),k1) = Xnew(J(i),k1) + kmove*nx; 
            Xnew(J(i)+N,k1) = Xnew(J(i)+N,k1) + kmove*ny; 
          end
        
          if idebug
          figure(1);
          title('Moved')
          plot(Xnew(J(i),k1),Xnew(J(i)+N,k1),'mo','markerfacecolor','m','markersize',8);
          axis equal
          pause
          end
        end % if sqrt(distx^2+disty^2) <= htol/2
        
      end % end for i = 1: numel(J)
      % Check if k1 is self-intersecting. If so, undo jiggling
      [xIntersect,~,~] = o.selfintersect(Xnew(:,k1));
      if ~isempty(xIntersect)
        disp('Jiggling causes self-intersection, so undo it')
        Xnew(:,k1) = X0;
      end
    end % end if(numel(j)~=0)
  end %end for k2 = K
  
  if confined
  % Go over each wall, check its interaction with vesicles
  % handle wall-vesicle interaction
  X0 = Xnew(:,k1); % save the original shape
  for k2 = 1 : nvbd
    J = find(zoneW2V{k2}(:,k1) == 1);
    % set of points on vesicle k1 close to wall k2
    if (numel(J) ~= 0)
      for i = 1 : numel(J)
        % distance vector from k1 to projection of the point J(i) on wall k2
        distx = (Xnew(J(i),k1) - nearestW2V{k2}(J(i),k1));
        disty = (Xnew(J(i)+N,k1)- nearestW2V{k2}(J(i)+N,k1));
        dist = sqrt(distx^2+disty^2);
        if dist <= htolWall
          if idebug
          figure(1);clf;hold on;
          title('Detected')
          plot(Xnew(1:end/2,:),Xnew(end/2+1:end,:),'r','linewidth',2)
          plot(Xwalls(1:end/2,:),Xwalls(end/2+1:end,:),'k-o','markerfacecolor','k','linewidth',2)
          plot(Xnew(J(i),k1),Xnew(J(i)+N,k1),'bo','markerfacecolor','b','markersize',8);
          plot(nearestW2V{k2}(J(i),k1),nearestW2V{k2}(J(i)+N,k1),'go','markerfacecolor','g','markersize',8)
          axis equal
          pause
          end
        
          % normal vector
          nx = distx/dist;
          ny = disty/dist;
          % solve (distx+kmove*nx)^2+(disty+kmove*ny)^2 = (beta*htol)^2
          % do not divide by 2 since wall does not move
          quadB = 2*(distx*nx+disty*ny); 
          quadC = distx^2+disty^2-beta^2*htolWall^2;
          kmove = max(roots([1 quadB quadC]));
          if kmove <= 0 || ~isreal(kmove)
            kmove = beta*htolWall-dist; 
          end
          
          % move the points in the normal direction by kmove
          if 1 % move several points using a smooth function
            % find nearby points' idcs (the middle point is J(i))
            idcs = mod((J(i)-interpOrder/2:J(i)+interpOrder/2)'-1,N)+1;
            % distribute the motion, make sure pth point move kmove*nx
            Xnew(idcs,k1) = Xnew(idcs,k1) + yfilt*kmove*nx;
            % for y-component
            % distribute the motion, make sure pth point move kmove*nx
            Xnew(idcs+N,k1) = Xnew(idcs+N,k1) + yfilt*kmove*ny;
          else % move only one point
            Xnew(J(i),k1) = Xnew(J(i),k1) + kmove*nx;
            Xnew(J(i)+N,k1) = Xnew(J(i)+N,k1) + kmove*ny; 
          end
        
          if idebug
          figure(1);
          title('Moved')
          plot(Xnew(J(i),k1),Xnew(J(i)+N,k1),'mo','markerfacecolor','m','markersize',8);
          axis equal
          pause
          end
        end % if sqrt(distx^2+disty^2) <= htolWall(k2)/2
        
      end % end for i = 1 : numel(J)
      % Check if k2 is self-intersecting. If so, undo jiggling
      [xIntersect,~,~] = o.selfintersect(Xnew(:,k1));
      if ~isempty(xIntersect)
        disp('Jiggling causes self-intersection, so undo it')
        Xnew(:,k1) = X0;
      end
    end % if numel(J) ~= 0
  end % end for k2 = 1 : 2
  end % if o.confined
      
end % end for k1 = 1 : nv
end % fixCollisionZinc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = fixCollisionsDLD(o,X,XwallsInt)

Xorig = X;
enScale = 1.1;
N = size(X,1)/2; nv = size(X,2); 
% for k = 1 : nv
%   X(:,k) = [enScale*(X(1:end/2,k)-mean(X(1:end/2,k))); enScale*(X(end/2+1:end,k)-mean(X(end/2+1:end,k)))];
% end

nvbdInt = size(XwallsInt,2); NbdInt = size(XwallsInt,1)/2;
% do not expect vesicle to get close to exterior wall (so, won't jiggle for
% that)

% near-zone defined as tol*spacing
[~,~,len] = o.geomProp(X);
beta = 1.1; % small buffer
htol = max(len)/N;
[~,~,lenW] = o.geomProp(XwallsInt);
htolWallInt = min(lenW/NbdInt);
htolWallInt = htol;

% build vesicle class
vesicle = capsules(X,[],[],[],[],1);
vesicle.setUpRate();

% get near structure
NearV2V = vesicle.getZone([],1);

wallsInt = capsules(XwallsInt,[],[],[],[],1);
wallsInt.setUpRate();  
[~,NearWint2V] = wallsInt.getZone(vesicle,2);
zoneWint2V = NearWint2V.zone;
nearestWint2V = NearWint2V.nearest;
icpWint2V = NearWint2V.icp;

zoneV2V = NearV2V.zone;
nearestV2V = NearV2V.nearest;
icpV2V = NearV2V.icp;

% Smooth function to move points
interpOrder = 6; % include 1/8 if N points to move (std = 1/5 of order)
% So, the region of contact is L/8
xfilt = linspace(-interpOrder/2,interpOrder/2,interpOrder+1);
yfilt = exp(-xfilt'.^2/(2*1.5));
% normalize the function so that the middle point is 1
yfilt = yfilt/yfilt(interpOrder/2+1); 

% Go over each vesicle, check its interaction with other vesicles
Xnew = X;
for k1 = 1 : nv
  
  % vesicle-vesicle interaction first
  K = [(1:k1-1) (k1+1:nv)];
  X0 = Xnew(:,k1); % save the original shape
  for k2 = K
    J = find(zoneV2V{k2}(:,k1) == 1);
    % set of points on vesicle k2 close to vesicle k1   
    if (numel(J) ~= 0)
      for i = 1:numel(J)
        % distance vector from a point on k1 to its projection on k2
        distx = (Xnew(J(i),k1) - nearestV2V{k2}(J(i),k1));
        disty = (Xnew(J(i)+N,k1) - nearestV2V{k2}(J(i)+N,k1));
        dist = sqrt(distx^2+disty^2);
        if dist <= htol
          % normal direction between the point on k2 and the projection of k2
          % on vesicle k1
          nx = distx/dist; ny = disty/dist;
          % Move only 1/2 of htol*beta(buffer zone)
          % solve (distx/2+kmove*nx)^2+(disty/2+kmove*ny)^2 =
          % (beta*htol/2)^2
          quadB = (distx*nx+disty*ny); quadC = dist^2/4-beta^2*htol^2/4;
          kmove = max(roots([1 quadB quadC]));
          if kmove <= 0 || ~isreal(kmove)
            kmove = beta*htol/2-dist/2; 
          end
        
          % move the points in the normal direction by kmove
          if 1 % move several points using a smooth function
            % find nearby points' idcs (the middle point is J(i))
            idcs = mod((J(i)-interpOrder/2:J(i)+interpOrder/2)'-1,N)+1;
            % distribute the motion, make sure pth point move kmove*nx
            Xnew(idcs,k1) = Xnew(idcs,k1) + yfilt*kmove*nx;
            % for y-component
            % distribute the motion, make sure pth point move kmove*nx
            Xnew(idcs+N,k1) = Xnew(idcs+N,k1) + yfilt*kmove*ny;
          else % move only one point
            Xnew(J(i),k1) = Xnew(J(i),k1) + kmove*nx; 
            Xnew(J(i)+N,k1) = Xnew(J(i)+N,k1) + kmove*ny; 
          end
        end % if sqrt(distx^2+disty^2) <= htol/2
        
      end % end for i = 1: numel(J)
      % Check if k1 is self-intersecting. If so, undo jiggling
      [xIntersect,~,~] = o.selfintersect(Xnew(:,k1));
      if ~isempty(xIntersect)
        disp('Jiggling causes self-intersection, so undo it')
        Xnew(:,k1) = X0;
      end
    end % end if(numel(j)~=0)
  end %end for k2 = K
  
  % Go over each wall, check its interaction with vesicles
  % handle wall-vesicle interaction
  X0 = Xnew(:,k1); % save the original shape
  for k2 = 1 : nvbdInt
    J = find(zoneWint2V{k2}(:,k1) == 1);
    % set of points on vesicle k1 close to wall k2
    if (numel(J) ~= 0)
      for i = 1 : numel(J)
        % distance vector from k1 to projection of the point J(i) on wall k2
        distx = (Xnew(J(i),k1) - nearestWint2V{k2}(J(i),k1));
        disty = (Xnew(J(i)+N,k1)- nearestWint2V{k2}(J(i)+N,k1));
        dist = sqrt(distx^2+disty^2);
        if dist <= htolWallInt
          % normal vector
          nx = distx/dist;
          ny = disty/dist;
          % solve (distx+kmove*nx)^2+(disty+kmove*ny)^2 = (beta*htol)^2
          % do not divide by 2 since wall does not move
          quadB = 2*(distx*nx+disty*ny); 
          quadC = distx^2+disty^2-beta^2*htolWallInt^2;
          kmove = max(roots([1 quadB quadC]));
          if kmove <= 0 || ~isreal(kmove)
            kmove = beta*htolWallInt-dist; 
          end
          
          % move the points in the normal direction by kmove
          if 1 % move several points using a smooth function
            % find nearby points' idcs (the middle point is J(i))
            idcs = mod((J(i)-interpOrder/2:J(i)+interpOrder/2)'-1,N)+1;
            % distribute the motion, make sure pth point move kmove*nx
            Xnew(idcs,k1) = Xnew(idcs,k1) + yfilt*kmove*nx;
            % for y-component
            % distribute the motion, make sure pth point move kmove*nx
            Xnew(idcs+N,k1) = Xnew(idcs+N,k1) + yfilt*kmove*ny;
          else % move only one point
            Xnew(J(i),k1) = Xnew(J(i),k1) + kmove*nx;
            Xnew(J(i)+N,k1) = Xnew(J(i)+N,k1) + kmove*ny; 
          end
        end % if sqrt(distx^2+disty^2) <= htolWall(k2)/2
        
      end % end for i = 1 : numel(J)
      % Check if k2 is self-intersecting. If so, undo jiggling
      [xIntersect,~,~] = o.selfintersect(Xnew(:,k1));
      if ~isempty(xIntersect)
        disp('Jiggling causes self-intersection, so undo it')
        Xnew(:,k1) = X0;
      end
    end % if numel(J) ~= 0
  end % end for k2 = 1 : 2
      
end % end for k1 = 1 : nv
% for k = 1 : nv
%   Xnew(:,k) = [1/enScale*(Xnew(1:end/2,k)-mean(Xnew(1:end/2,k))); 1/enScale*(Xnew(end/2+1:end,k)-mean(Xnew(end/2+1:end,k)))];
% end
end % fixCollisionDLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew] = correctAreaAndLength(o,X,timeTolerance,om)
% Xnew = correctAreaAndLength(X,a0,l0) changes the shape of the vesicle
% by finding the shape Xnew that is closest to X in the L2 sense and
% has the same area and length as the original shape

% tolConstraint (which controls area and length) comes from the area-length
% tolerance for time adaptivity.

a0 = om.area;
l0 = om.length;
[~,at,lt] = o.geomProp(X);
eAt = abs(at-a0)./a0;
eLt = abs(lt-l0)./l0;

N  = size(X,1)/2;
  
tolConstraint = 1e-2; % 1 percent error in constraints
% tolConstraint = timeTolerance;
tolFunctional = 1e-2; % Allowed to change shape by 1 percent

options = optimset('Algorithm','sqp','TolCon',tolConstraint,...
    'TolFun',tolFunctional,'display','off','MaxFunEvals',3000);

% Options for Algorithm are:
% 'active-set', 'interior-point', 'interior-point-convex' 'sqp'

Xnew = zeros(size(X));

for k = 1:size(Xnew,2)
    
  minFun = @(z) 1/N*min(sum((z - X(:,k)).^2));
  [Xnew(:,k),~,iflag] = fmincon(minFun,X(:,k),[],[],[],[],[],[],...
      @(z) o.nonlcon(z,a0(k),l0(k)),options);
  if iflag~=1 && iflag~=2
    message = ['Correction scheme failed, do not correct at this step'];
    om.writeMessage(message,'%s\n')
    Xnew(:,k) = X(:,k);
  end
  % if fmincon fails, keep the current iterate for this time step.
  % Hopefully it'll be corrected at a later step.
  
end

% Looping over vesicles, correct the area and length of each vesicle

end % correctAreaAndLength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,iFailCorrection] = correctAreaAndLength2(o,X,a0,l0)
% Xnew = correctAreaAndLength(X,a0,l0) changes the shape of the vesicle
% by finding the shape Xnew that is closest to X in the L2 sense and
% has the same area and length as the original shape

% tolConstraint (which controls area and length) comes from the area-length
% tolerance for time adaptivity.

[~,at,lt] = o.geomProp(X);

N  = size(X,1)/2;

tolConstraint = 1e-2; % 1 percent error in constraints
% tolConstraint = timeTolerance;
tolFunctional = 1e-2; % Allowed to change shape by 1 percent

options = optimset('Algorithm','sqp','TolCon',tolConstraint,...
    'TolFun',tolFunctional,'display','off','MaxFunEvals',3000);

% Options for Algorithm are:
% 'active-set', 'interior-point', 'interior-point-convex' 'sqp'

Xnew = zeros(size(X));

% flag for failing correction
iFailCorrection = false;
for k = 1:size(Xnew,2)
    
  minFun = @(z) 1/N*min(sum((z - X(:,k)).^2));
  [Xnew(:,k),~,iflag] = fmincon(minFun,X(:,k),[],[],[],[],[],[],...
      @(z) o.nonlcon(z,a0(k),l0(k)),options);
  if iflag~=1 && iflag~=2
    iFailCorrection = true;
    message = ['Correction scheme failed, do not correct at this step'];
    disp(message)
    Xnew(:,k) = X(:,k);
  end
  % if fmincon fails, keep the current iterate for this time step.
  % Hopefully it'll be corrected at a later step.
  
end
% Looping over vesicles, correct the area and length of each vesicle

end % correctAreaAndLength2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,iFailCorrection] = correctAreaAndLength3(o,X,a0,l0,downTo)
% Xnew = correctAreaAndLength(X,a0,l0) changes the shape of the vesicle
% by finding the shape Xnew that is closest to X in the L2 sense and
% has the same area and length as the original shape

% tolConstraint (which controls area and length) comes from the area-length
% tolerance for time adaptivity.

[~,at,lt] = o.geomProp(X);

N  = size(X,1)/2;

tolConstraint = 1e-2; % 1 percent error in constraints
% tolConstraint = timeTolerance;
tolFunctional = 1e-2; % Allowed to change shape by 1 percent

options = optimset('Algorithm','sqp','TolCon',tolConstraint,...
    'TolFun',tolFunctional,'display','off','MaxFunEvals',3000);

% Options for Algorithm are:
% 'active-set', 'interior-point', 'interior-point-convex' 'sqp'

Xnew = zeros(size(X));

% flag for failing correction
iFailCorrection = false;
for k = 1:size(Xnew,2)
  Xd = [interpft(X(1:end/2,k),downTo);interpft(X(end/2+1:end,k),downTo)];  
  minFun = @(z) 1/N*min(sum((z - Xd).^2));
  [XnewD,~,iflag] = fmincon(minFun,Xd,[],[],[],[],[],[],...
      @(z) o.nonlcon(z,a0(k),l0(k)),options);
  if iflag~=1 && iflag~=2
    iFailCorrection = true;
    message = ['VesID: ' num2str(k) 'Correction scheme failed, do not correct at this step'];
    disp(message)
    Xnew(:,k) = X(:,k);
  else
    Xnew(:,k) = [interpft(XnewD(1:end/2),N);interpft(XnewD(end/2+1:end),N)];  
  end
  % if fmincon fails, keep the current iterate for this time step.
  % Hopefully it'll be corrected at a later step.
  
end
% Looping over vesicles, correct the area and length of each vesicle

end % correctAreaAndLength2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cIn,cEx] = nonlcon(o,X,a0,l0)
% [cIn,cEx] = nonlcon(X,a0,l0) is the non-linear constraints required
% by fmincon

[~,a,l] = o.geomProp(X);

cIn = [];
% new inequalities in the constraint
cEx = [(a-a0)/a0 (l-l0)/l0];
% want to keep the area and length the same

end % nonlcon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew] = alignCenterAngle(o,Xorg,X)
% Xnew = alignCenterAngle(o,Xorg,X) uses
% rigid body translation and rotation to match X having the corrected area 
% and length but wrong center and inclination angle with Xorg having the 
% right center and IA but wrong area and length. So that, Xnew has the
% correct area,length,center and inclination angle.
N = size(X,1)/2;
nv = size(X,2);
Xnew = zeros(size(X));

for k = 1 : nv
initMean = [mean(Xorg(1:end/2,k)); mean(Xorg(end/2+1:end,k))];
newMean = [mean(X(1:end/2,k)); mean(X(end/2+1:end,k))];

initAngle = o.getIncAngle2(Xorg(:,k));
newAngle = o.getIncAngle2(X(:,k));
if newAngle > pi
  newAngle2 = newAngle-pi;
else
  newAngle2 = newAngle+pi;
end
newAngles = [newAngle;newAngle2];
diffAngles = abs(initAngle-newAngles); [~,id] = min(diffAngles);
newAngle = newAngles(id);

% move to (0,0) new shape
Xp = [X(1:end/2,k)-newMean(1); X(end/2+1:end,k)-newMean(2)];

% tilt it to the original angle
XpNew = zeros(size(Xp)); thet = 0;% -newAngle+initAngle;
XpNew(1:end/2) = Xp(1:end/2)*cos(thet)-Xp(end/2+1:end)*sin(thet);
XpNew(end/2+1:end) = Xp(1:end/2)*sin(thet)+Xp(end/2+1:end)*cos(thet);

% move to original center
Xnew(:,k) = [XpNew(1:end/2)+initMean(1); XpNew(end/2+1:end)+initMean(2)];
end

% OLD VERSION (the one in 2018 JCP paper about low-resolutions)
% N = size(X,1)/2;
% nv = size(X,2);
% Nup = 96;
% 
% Xnew = zeros(size(X));
% 
% %tolerance for the error 
% % tolFun = 1e-6;
% 
% for k = 1 : nv
% 
%     % Find the centers
%     xmean_org = mean(interpft(Xorg(1:end/2,k),Nup));
%     ymean_org = mean(interpft(Xorg(end/2+1:end,k),Nup));
%     
%     xmean = mean(interpft(X(1:end/2,k),Nup));
%     ymean = mean(interpft(X(end/2+1:end,k),Nup));
%     
%     % Decenter and form H
%     xmean = xmean(1,ones(N,1)); ymean = ymean(1,ones(N,1));
%     xmean_org = xmean_org(1,ones(N,1)); ymean_org = ymean_org(1,ones(N,1));
%     H = ([X(1:N,k)';X(N+1:2*N,k)']-[xmean;ymean])*...
%         ([Xorg(1:N,k)';Xorg(N+1:2*N,k)']-[xmean_org;ymean_org])';
%     
%     [U,~,V] = svd(H);
%     
%     % Rotation matrix
%     R = V*U';
%     if det(R)<0
%         R(:,2) = R(:,2)*(-1);
%     end
%     
%     % Find translation matrix    
%     t = -R*[xmean(1);ymean(1)] + [xmean_org(1);ymean_org(1)];
% 
%     % Align center and angle    
%     t = t(:,ones(N,1));
%     xnew = R*[X(1:N,k)';X(N+1:2*N,k)'] + t;
%     Xnew(:,k) = [xnew(1,:)';xnew(2,:)'];
%     
% end


end % alignCenterAngle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = filterShape(o,X)
% delete high frequencies from the vesicle shape
N = size(X,1)/2;
nv = size(X,2);
modes = [(0:N/2-1) (-N/2:-1)];

[x,y] = o.getXY(X);

for k = 1:nv
  z = x(:,k) + 1i*y(:,k);
  z = fft(z);
%   z(abs(modes) > 2/3*N/2) = 0;
  z(abs(modes) > 1/2*N/2) = 0;
  z = ifft(z);
  X(:,k) = o.setXY(real(z),imag(z));
end

end % filterShape


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,filtered,ratio] = adapFilterShape(o,X)
% delete high frequencies from the vesicle shape
N = size(X,1)/2;
nv = size(X,2);

modes = [(0:N/2-1) (-N/2:-1)];
ratio = zeros(nv,1);

if N <= 32
  nModes = 12;
  tol = 1e-2;
elseif N == 64 
  nModes = 12; % 16
  tol = 8e-3; % 5e-3
else
  nModes = 16;
  tol = 1e-3;
end

[x,y] = o.getXY(X);

% to keep list of filtered shapes
filtered = zeros(nv,1);

for k = 1:nv 
  z = x(:,k) + 1i*y(:,k);
  z = fft(z);
  

  lowEnergy = sqrt(sum(abs(z(abs(modes) < nModes & modes ~= 0)).^2));
  highEnergy = sqrt(sum(abs(z(abs(modes) >= nModes)).^2));
  ratio(k) = highEnergy/lowEnergy;
  
  if ratio(k) > tol 
    filtered(k) = 1;
    z(abs(modes) >= nModes) = 0;
    z = ifft(z);
    X(:,k) = o.setXY(real(z),imag(z));
    
    % this filtering can cause self-intersection, if so, undo it
    [xIntersect,~,~] = o.selfintersect(X(:,k));
    if ~isempty(xIntersect)
      X(:,k) = [x(:,k);y(:,k)];
    end
  end % if ratio(k) >= 2e-2
end

end % adapFilterShape

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [derivA,derivL] = computeAreaLengthDerivs(o,X,Xold,dt)
% [derivA,derivL] = computeAreaLengthDerivs(o,X,Xold,dt) computes dA/dt 
% and dL/dt based on the area, length computation through parameterization
N = size(X,1)/2;
nv = size(X,2);
Nup = 96;

[x,y] = o.getXY(X);
[xOld,yOld] = o.getXY(Xold);

% Upsample
xUp = interpft(x,Nup); yUp = interpft(y,Nup);
xOldUp = interpft(xOld,Nup); yOldUp = interpft(yOld,Nup);

% Compute velocity on upsampled grid using finite differences
vel = ([xUp;yUp]-[xOldUp;yOldUp])/dt;

derivA = zeros(nv,1);
derivL = zeros(nv,1);

for k = 1 : nv  
    % Compute dx/da of vesicle k
    [DxOld,DyOld] = o.getDXY([xOldUp(:,k);yOldUp(:,k)]);
    % Compute dv/da of vesicle k
    [DvelX,DvelY] = o.getDXY(vel(:,k));
    % Get x and y components of velocity of vesicle k
    [velx,vely] = o.getXY(vel(:,k));
    
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
function [X,u,sigma] = redistributeArcLength(o,X,u,sigma)
% [X,u,sigma] = resdistributeArcLength(o,X,u,sigma) resdistributes
% the vesicle shape eqiuspaced in arclength and adjusts the tension and
% velocity according to the new parameterization

N = size(X,1)/2;
nv = size(X,2);
modes = [(0:N/2-1) (-N/2:-1)];
jac = o.diffProp(X);
jac1 = jac;
tol = 1e-10;

u = [];
sigma = [];

%figure(1); clf;
%plot(X(1:end/2),X(end/2+1:end),'r-o')
%hold on
%axis equal
%axis(1.1*[-3 3 -3 3])
for k = 1:nv
  if norm(jac(:,k) - mean(jac(:,k)),inf) > tol*mean(jac(:,k))
    theta = o.arcLengthParameter(X(1:end/2,k),...
        X(end/2+1:end,k));
    zX = X(1:end/2,k) + 1i*X(end/2+1:end,k);
    zXh = fft(zX)/N;
    zX = zeros(N,1);
    for j = 1:N
      zX = zX + zXh(j)*exp(1i*modes(j)*theta);
    end
    X(:,k) = o.setXY(real(zX),imag(zX));
    if nargin > 2
      zu = u(1:end/2,k) + 1i*u(end/2+1:end,k);
      zuh = fft(zu)/N;
      sigmah = fft(sigma(:,k))/N;
      zu = zeros(N,1);
      sigma(:,k) = zeros(N,1);
      for j = 1:N
        zu = zu + zuh(j)*exp(1i*modes(j)*theta);
        sigma(:,k) = sigma(:,k) + sigmah(j)*exp(1i*modes(j)*theta);
      end
      sigma = real(sigma);
      u(:,k) = o.setXY(real(zu),imag(zu));
    else
      u = [];
      sigma = [];
    end
    % redistribute the vesicle positions and tension so that it is
    % equispaced in arclength
  end
end
%plot(X(1:end/2),X(end/2+1:end),'b*')
%jac2 = o.diffProp(X);
%figure(2); clf;
%semilogy(abs(jac1 - mean(jac1)),'r')
%hold on
%semilogy(abs(jac2 - mean(jac2)),'b')
%norm(jac1 - mean(jac1),inf)
%norm(jac2 - mean(jac2),inf)
%ylim([1e-3 1e1])
%pause

end % redistributeArcLength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta,arcLength] = arcLengthParameter(o,x,y)
% theta = arcLengthParamter(o,x,y) finds a discretization of parameter
% space theta so that the resulting geometry will be equispaced in
% arclength

uprate = 1;
N = numel(x);
Nup = uprate*N;
t = (0:Nup-1)'*2*pi/Nup; % this is not correct when you iterate
x = interpft(x,Nup);
y = interpft(y,Nup);
[~,~,len] = o.geomProp([x;y]);
% find total perimeter
[Dx,Dy] = o.getDXY([x;y]);
% find derivative
arc = sqrt(Dx.^2 + Dy.^2);
arch = fft(arc);
modes = -1i./[(0:Nup/2-1) 0 (-Nup/2+1:-1)]';
modes(1) = 0;
modes(Nup/2+1) = 0;
arcLength = real(ifft(modes.*arch) - sum(modes.*arch/Nup) + ...
    arch(1)*t/Nup);
% arclength from 0 to t(n) in parameter space
%clf
%semilogy(t,abs(arcLength -len/2/pi*t),'b-o')
%ylim([1e-20 1])
%pause

z1 = [arcLength(end-6:end)-len;arcLength;arcLength(1:7)+len];
z2 = [t(end-6:end)-2*pi;t;t(1:7)+2*pi];
% put in some overlap to account for periodicity

theta = [interp1(z1,z2,(0:N-1)'*len/N,'spline')];

end % arcLengthParamter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,niter] = reparametrize(o,X,dX,uprate,maxIter)
% [X,niter] = reparametrize applies the reparametrization with
% minimizing the energy in the high frequencies (Veerapaneni et al. 2011, 
% doi: 10.1016/j.jcp.2011.03.045, Section 6).

% Plot how the coefficients decay and the points move
plotError = false;

% Decay power (k^pow)
pow = 4;

N     = size(X,1)/2; % # of points per vesicle
Nup   = N*uprate;
nv    = size(X,2)  ; % # of vesicles

niter = ones(nv,1); % store # of iterations per vesicle
tolg = 1e-3;
if isempty(dX)
  [~,~,len] = o.geomProp(X);    
  dX = len/N;
  toly = 1e-5*dX;
else
  normDx = sqrt(dX(1:end/2,:).^2+dX(end/2+1:end,:).^2);
  toly = 1e-3*min(normDx(:));  
end


beta = 0.1;
dtauOld = 0.05;

for k = 1:nv
    
    % Get initial coordinates of kth vesicle (upsample if necessary)
    x0 = interpft(X(1:end/2,k),Nup);
    y0 = interpft(X(end/2+1:end,k),Nup);
    
    % Compute initial projected gradient energy
    g0 = o.computeProjectedGradEnergy(x0,y0,pow);  
    x = x0; y = y0; g = g0;
    
    % Explicit reparametrization
    while niter(k) <= maxIter
        dtau = dtauOld;
        xn = x - g(1:end/2)*dtau; yn = y - g(end/2+1:end)*dtau;
        gn = o.computeProjectedGradEnergy(xn,yn,pow);
        while norm(gn) > norm(g)
            dtau = dtau*beta;
            xn = x - g(1:end/2)*dtau; yn = y - g(end/2+1:end)*dtau;
            gn = o.computeProjectedGradEnergy(xn,yn,pow);
        end
        dtauOld = dtau*1/beta;
%         xn = x - g(1:end/2)*dtau;
%         yn = y - g(end/2+1:end)*dtau;
%         gn = o.computeProjectedGradEnergy(xn,yn,pow);
        
        if norm(gn) < max(toly/dtau,tolg*norm(g0))
            break
        end
        
        x = xn; y = yn; g = gn;
        niter(k) = niter(k)+1;
        
        if plotError
            % Check if it decays with desired rate
            modes = (0:Nup/2-1);
            zX = xn + 1i*yn;
            zXh = abs(fft(zX)/Nup);
            subplot(1,2,1)
            plot([x0;x0(1)],[y0;y0(1)],'r-o')
            hold on
            plot([xn;xn(1)],[yn;yn(1)],'b-o')
            title(['# of iter.: ' num2str(niter(k))])
            hold off
            axis equal
            subplot(1,2,2)
            plot(modes,zXh(1:Nup/2),'b-o','linewidth',2)
            hold on
            plot(modes,zXh(2)*1./modes.^pow,'r-o','linewidth',2)
            hold off
            axis square
            title('termination condition')
            pause(0.01)
         end
    end
    % Combine x and y
%     xnh = fft(xn)/uprate;
%     ynh = fft(yn)/uprate;
%     X(:,k) = [real(ifft([xnh(1:N/2);xnh(Nup-N/2+1:Nup)]));...
%         real(ifft([ynh(1:N/2);ynh(Nup-N/2+1:Nup)]))];
    X(:,k) = [interpft(xn,N);interpft(yn,N)];
end

end % end reparamEnergyMin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g] = computeProjectedGradEnergy(o,x,y,pow)
% g = computeProjectedGradEnergy(o,x,y) computes the projected gradient of
% the energy of the surface. We use this in reparamEnergyMin(o,X). For the
% formulation see (Veerapaneni et al. 2011 doi: 10.1016/j.jcp.2011.03.045,
% Section 6)

N = numel(x);

% to be used in computing gradE 
modes = [(0:N/2-1) (-N/2:-1)]'; 

% get tangent vector at each point (tang_x;tang_y) 
[~,tang] = o.diffProp([x;y]);
% get x and y components of normal vector at each point
nx = tang(N+1:2*N);
ny = -tang(1:N);

% Compute gradE
% first, get Fourier coefficients
zX = x + 1i*y;
zXh = fft(zX)/N;
% second, compute zX with a_k = k^pow
zX = ifft(N*zXh.*abs(modes).^pow);

% Compute Energy
gradE = [real(zX);imag(zX)]; % [gradE_x;gradE_y]

% ######################################################################
% % Compute projected gradient g (there might be a way to avoid for loop?)
% g = zeros(2*N,1); 
% for j = 1 : N
%     % compute x and y components of projected gradient at each point
%     gloc = (eye(2) - [nx(j);ny(j)]*[nx(j);ny(j)]')* ...
%         [gradE(j);gradE(j+N)];
%     g(j) = gloc(1); 
%     g(j+N) = gloc(2);
% end
% ######################################################################
% A dyadic product property (a (ban) a)b = a(a.b) can be used to avoid the
% for loop as follows
normals = [nx;ny];
% do the dot product n.gradE
prod = normals.*gradE;
dotProd = prod(1:N)+prod(N+1:2*N);
% do (I-(n ban n))gradE = gradE - n(n.gradE) for each point
g = gradE - normals.*[dotProd;dotProd];
% norm(g-gAltern)
% pause
end % end computeProjectedGradEnergy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,niter,u,sigma] = equiArcLengthDist(o,X,u,sigma)
% function [X,niter,u,sigma] = equiArcLengthDist(o,X,u,sigma) is another
% implementation of Bryan's preliminary implementation of redistributing
% the points to have equidistributed arc length. 

N     = size(X,1)/2; % # of points per vesicle
nv    = size(X,2)  ; % # of vesicles

tol     = 1e-6; % tolarence for difference btw arcLength distribution and
                % theoretical equidistributed arcLength
maxIter = 1000; 

% Upsampling
uprate = 2;
Nup    = uprate*N;
theta0 = (0:Nup-1)'*2*pi/Nup; 
modes = [(0:Nup/2-1) (-Nup/2:-1)];

% Store # of iterations
niter = zeros(nv,1);  

for k = 1:nv
  theta = theta0;
  Xup = [interpft(X(1:end/2,k),Nup);interpft(X(end/2+1:end,k),Nup)];
  [arcLength,len] = o.computeArcLengthArray ...
          (Xup(1:end/2),Xup(end/2+1:end),theta);    
    
  while (norm(arcLength-(0:Nup-1)'*len/Nup)/...
      norm((0:Nup-1)'*len/Nup)>tol && niter(k,:)<=maxIter)
    z1 = [arcLength(end-6:end)-len;arcLength;arcLength(1:7)+len];
    z2 = [theta(end-6:end)-2*pi;theta;theta(1:7)+2*pi];
    % put in some overlap to account for periodicity
       
    theta = interp1(z1,z2,(0:Nup-1)'*len/Nup,'linear');
       
    zX = Xup(1:end/2) + 1i*Xup(end/2+1:end);
    zXh = fft(zX)/Nup;
    zX = zeros(Nup,1);
    for j = 1:Nup
      zX = zX + zXh(j)*exp(1i*modes(j)*theta);
    end
    Xup = o.setXY(real(zX),imag(zX));
        
    if uprate ~= 1
      X(1:end/2,k) = interpft(Xup(1:end/2),N);
      X(end/2+1:end,k) = interpft(Xup(end/2+1:end),N);
    else 
      X = Xup;
    end
        
    [arcLength,len] = o.computeArcLengthArray ...
        (Xup(1:end/2,k),Xup(end/2+1:end,k),theta);
    niter(k,1) = niter(k,1) + 1;
  end
  if nargin > 2
    zu = u(1:end/2,k) + 1i*u(end/2+1:end,k);
    zuh = fft(zu)/N;
    sigmah = fft(sigma(:,k))/N;
    zu = zeros(N,1);
    sigma(:,k) = zeros(N,1);
    for j = 1:N
      zu = zu + zuh(j)*exp(1i*modes(j)*theta);
      sigma(:,k) = sigma(:,k) + sigmah(j)*exp(1i*modes(j)*theta);
    end
    sigma = real(sigma);
    u(:,k) = o.setXY(real(zu),imag(zu));
  end
  % redistribute the vesicle positions and tension so that it is
  % equispaced in arclength
end

end % equiArcLengthDist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [arcLength,len] = computeArcLengthArray(o,x,y,t)
% function [arcLength,len] = computeArcLengthArray(o,x,y,t) computes
% arcLength s(t(j)) as a function of parametrization t(j) and s(2*pi).

% # of points per vesicle
N = numel(x);
% Find total perimeter
[~,~,len] = o.geomProp([x;y]);
% Find derivative of x and y w.r.t. parametrization t
[Dx,Dy] = o.getDXY([x;y]);
% Compute x'
arc = sqrt(Dx.^2 + Dy.^2);
% Compute Fourier coefficients of x' to integrate
arch = fft(arc);
modes = -1i./[(0:N/2-1) 0 (-N/2+1:-1)]';
modes(1) = 0;
modes(N/2+1) = 0;
arcLength = real(ifft(modes.*arch) - sum(modes.*arch/N) + ...
    arch(1)*t/N);

end % computeArcLengthArray
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x0,y0,segments]=selfintersect(o,X)

%SELFINTERSECT Self-intersections of a curve.
%
%    [X0,Y0,SEGMENTS] = SELFINTERSECT(X,Y) computes the locations where
%    a curve self-intersects in a fast and robust way.
%    The curve can be broken with NaNs or have vertical segments.
%    Segments of the curve involved in each of the self-interesections are
%    also provided.
%
%    Vectors X and Y are equal-length vectors of at least four points defining
%    the curve.
%    X0 and Y0 are column vectors with the x- and y- coordinates, respectively
%    of the N self-intersections.
%    SEGMENTS is an N x 2 matrix containing the pairs of segments involved in
%    each self-intersection.
%
%    This program uses the theory of operation of the file "Fast and Robust Curve
%    Intersections" submitted by Douglas M. Schwartz (intersections.m, F.Id: 11837).
%
%    Example of use
% 	 N=201;
% 	 th=linspace(-3*pi,4*pi,N);
% 	 R=1;
% 	 x=R*cos(th)+linspace(0,6,N);
% 	 y=R*sin(th)+linspace(0,1,N);
%    t0=clock;
%    [x0,y0,segments]=selfintersect(x,y)
% 	 etime(clock,t0)
%    plot(x,y,'b',x0,y0,'.r');
% 	 axis ('equal'); grid

x = X(1:end/2,1); y = X(end/2+1:end,1);

% Input checks.
error(nargchk(2,2,nargin))
% x and y must be vectors with same number of points (at least 4 for self-intersection).
if sum(size(x) > 3) ~= 1 || sum(size(y) > 3) ~= 1 || ...
		length(x) ~= length(y)
	error('X and Y must be equal-length vectors of at least 4 points.')
end

x0=[];
y0=[];
segments=[];

% Two similar curves are firstly created.
x1=x; x2=x;
y1=y; y2=y;

x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;

dxy1 = diff([x1 y1]);
dxy2 = diff([x2 y2]);

% Determine the combinations of i and j where the rectangle enclosing the
% i'th line segment of curve 1 overlaps with the rectangle enclosing the
% j'th line segment of curve 2.
[i,j] = find(repmat(min(x1(1:end-1),x1(2:end)),1,n2) <= ...
	repmat(max(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(max(x1(1:end-1),x1(2:end)),1,n2) >= ...
	repmat(min(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(min(y1(1:end-1),y1(2:end)),1,n2) <= ...
	repmat(max(y2(1:end-1),y2(2:end)).',n1,1) & ...
	repmat(max(y1(1:end-1),y1(2:end)),1,n2) >= ...
	repmat(min(y2(1:end-1),y2(2:end)).',n1,1));

% Removing coincident and adjacent segments.
remove=find(abs(i-j)<2);
i(remove)=[];
j(remove)=[];

% Removing duplicate combinations of segments.
remove=[];
for ii=1:size(i,1)
	ind=find((i(ii)==j(ii:end))&(j(ii)==i(ii:end)));
	remove=[remove;ii-1+ind];
end
i(remove)=[];
j(remove)=[];

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
i(remove) = [];
j(remove) = [];

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
i(remove) = [];
j(remove) = [];

% Initialize matrices.  We'll put the T's and B's in matrices and use them
% one column at a time.  For some reason, the \ operation below is faster
% on my machine when A is sparse so we'll initialize a sparse matrix with
% the fixed values and then assign the changing values in the loop.
n = length(i);
T = zeros(4,n);
A = sparse([1 2 3 4],[3 3 4 4],-1,4,4,8);
B = -[x1(i) x2(j) y1(i) y2(j)].';
index_dxy1 = [1 3];  %  A(1) = A(1,1), A(3) = A(3,1)
index_dxy2 = [6 8];  %  A(6) = A(2,2), A(8) = A(4,2)

% Loop through possibilities.  Set warning not to trigger for anomalous
% results (i.e., when A is singular).
warning_state = warning('off','MATLAB:singularMatrix');
try
	for k = 1:n
		A(index_dxy1) = dxy1(i(k),:);
		A(index_dxy2) = dxy2(j(k),:);
		T(:,k) = A\B(:,k);
	end
	warning(warning_state)
catch
	warning(warning_state)
	rethrow(lasterror)
end

% Find where t1 and t2 are between 0 and 1 and return the corresponding x0
% and y0 values.  Anomalous segment pairs can be segment pairs that are
% colinear (overlap) or the result of segments that are degenerate (end
% points the same).  The algorithm will return an intersection point that
% is at the center of the overlapping region.  Because of the finite
% precision of floating point arithmetic it is difficult to predict when
% two line segments will be considered to overlap exactly or even intersect
% at an end point.  For this algorithm, an anomaly is detected when any
% element of the solution (a single column of T) is a NaN.

in_range = T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1;
anomalous = any(isnan(T));
if any(anomalous)
	ia = i(anomalous);
	ja = j(anomalous);
	% set x0 and y0 to middle of overlapping region.
	T(3,anomalous) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
		min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1))))/2;
	T(4,anomalous) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
		min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1))))/2;
	x0 = T(3,in_range | anomalous).';
	y0 = T(4,in_range | anomalous).';
	i=i(in_range | anomalous);
	j=j(in_range | anomalous);
else
	x0 = T(3,in_range).';
	y0 = T(4,in_range).';
	i=i(in_range);
	j=j(in_range);
end

segments=sort([i,j],2);

end %selfintersect

end % methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = arcDeriv(f,m,sa,IK)
% f = arcDeriv(f,m,s,IK,col) is the arclength derivative of order m.
% f is a matrix of scalar functions (each function is a column)
% f is assumed to have an arbitrary parametrization
% sa = d s/ d a, where a is the aribtrary parameterization
% IK is the fourier modes which is saved and used to accelerate 
% this routine

for j=1:m
  f = sa.*ifft(IK.*fft(f));
end
f = real(f);

end % arcDeriv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zCoarse = restrict(z,Nfine,Ncoarse)
% zCoarse = restrict(z,Nfine,Ncoarse) restricts the periodic function
% z.  z has a column vector containing arbitrarly many periodic
% functions of size Nfine.  The output, zCoarse, has the same number of
% periodic copies at a grid with Ncoarse points.  NOTE: Nfine/Ncoarse
% must be a power of 2 if method == 'local'

method = 'spectral';
%method = 'local';

if Nfine == Ncoarse
  zCoarse = z;
  % if coarse and fine grids are the same, nothing to do
else
  nSecs = size(z,1)/Nfine;
  % number of functions that need to be restricted.  This is the number
  % of periodic functions that are stacked in each row

  if strcmp(method,'spectral')
    zFine = zeros(Nfine*nSecs,size(z,2));
    for j = 1:nSecs
      istart = (j-1)*Nfine + 1;
      iend = istart + Nfine - 1;
      istart2 = (j-1)*Ncoarse + 1;
      iend2 = istart2 + Ncoarse - 1;
      zh = fft(z(istart:iend,:));
      % take fft of a matrix of periodic columns
      zh = [zh(1:Ncoarse/2,:);zh(end-Ncoarse/2+1:end,:)]*Ncoarse/Nfine;
      % zero all the high frequencies and scale appropriatly
      zCoarse(istart2:iend2,:) = real(ifft(zh));
      % move to physical space and take real part as we are assuming
      % input is real
    end
  end
  % spectral restriction

  if strcmp(method,'local')
    nRestrict = log2(Nfine/Ncoarse);
    % number of halvings to do
    for k = 1:nRestrict
      Nfine = Nfine/2;
      zCoarse = zeros(size(z,1)/2,size(z,2));
      % create grid half the size
      for j = 1:nSecs
        istart = (j-1)*Nfine + 1;
        iend = istart + Nfine - 1;
        zCoarse(istart:iend,:) = 1/2*z(2*(istart:iend)-1,:);
        zCoarse(istart:iend,:) = zCoarse(istart:iend,:) + ...
          1/4*z(2*(istart:iend),:);
        zCoarse(istart+1:iend,:) = zCoarse(istart+1:iend,:) + ...
          1/4*z(2*(istart+1:iend)-2,:);
        zCoarse(istart,:) = zCoarse(istart,:) + 1/4*z(2*iend,:);
        % usual (1/4,1/2,1/4) restriction with periodicity built in
      end
      z = zCoarse;
    end
  end
  % local restriction
end

end % restrict

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zFine = prolong(z,Ncoarse,Nfine)
% zFine = prolong(z,Ncoarse,Nfine) prolongs the periodic function z.  z
% has a column vector containing arbitrarly many periodic functions of
% size Ncoarse.  The output, zfine, has the same number of periodic
% copies at a grid with Nfine points.  NOTE: Nfine/Ncoarse must be a
% power of 2

method = 'spectral';
%method = 'local';

if Nfine == Ncoarse
  zFine = z;
% if coarse and fine grids are the same, nothing to do
else
  nSecs = size(z,1)/Ncoarse;
  % number of functions that need to be prolonged.  This is the number of
  % periodic functions that are stacked in each row

  if strcmp(method,'spectral')
    zFine = zeros(Nfine*nSecs,size(z,2));
    for j = 1:nSecs
      istart = (j-1)*Ncoarse + 1;
      iend = istart + Ncoarse - 1;
      istart2 = (j-1)*Nfine + 1;
      iend2 = istart2 + Nfine - 1;

      zFine(istart2:iend2,:) = interpft(z(istart:iend,:),Nfine);
    end
  end
  % spectral prolongation

  if strcmp(method,'local')
    nProlong = log2(Nfine/Ncoarse);
    % number of doublings to do

    for k = 1:nProlong
      Ncoarse = Ncoarse*2;
      zFine = zeros(size(z,1)*2,size(z,2));
      for j = 1:nSecs
        istart = (j-1)*Ncoarse + 1;
        iend = istart + Ncoarse - 1;
        zFine(istart:2:iend,:) = z((istart+1)/2:iend/2,:);
        zFine(istart+1:2:iend,:) = 1/2*z((istart+1)/2:iend/2,:);
        zFine(istart+1:2:iend-2,:) = zFine(istart+1:2:iend-2,:) + ...
            1/2*z(((istart+1)/2:(iend-2)/2)+1,:);
        zFine(iend,:) = zFine(iend,:) + 1/2*z((istart+1)/2,:);
        % usual (1/2,1,1/2) prolongation with periodicity built in
      end
      z = zFine;
    end
  end
  % local prolongation
end

end % prolong

end % methods (Static)

end % classdef
