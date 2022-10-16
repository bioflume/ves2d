% script to check jumps for  stress of 2D Stokes single layer

clear all; clear globals; clf;
H=1e-3;   % distance from the boundary

fx=@(s) ones(size(s));  % hydro density, x -component
fy=@(s) ones(size(s));  % hydro density, y -component

% can also use variable density, jumps still get computed correctly.
if 0
fx=@(s) cos(2*s)+sin(3*s); fy=@(s) -0.2*cos(s)+sin(4*s); 
end

% shorthand for numerical integration of a function in 0-2*pi
int2pi = @(f) quad(f,0,2*pi,1e-6);

% symbolic representation of boundary
syms sS    % 0-2*pi parameterization
sR = 1+0.2*cos(2*sS)+0.3*sin(8*sS);  % radius
%sR = 1;
sX = sR.*cos(sS); %symbolic X and Y positions
sY = sR.*sin(sS); 
sDX=diff(sX,sS);  %symbolic derivatives
sDY=diff(sY,sS);  

% convert symbolic expressions to numerical functions
gx = matlabFunction(sX);    % position
gy = matlabFunction(sY);
tx = matlabFunction(sDX);   % unnormailized tangent
ty = matlabFunction(sDY);

% plot shape
figure(1),plot( gx([0:128]*2*pi/128), gy([0:128]*2*pi/128) ); axis equal; pause(0.05);

% given some target point x, get $r$ and other stokes kernel quantities.
rx=@(s,x) x(1)-gx(s);
ry=@(s,x) x(2)-gy(s);
nx=@(s) ty(s);               
ny=@(s) -tx(s);
rho2=@(s,x) rx(s,x).^2 + ry(s,x).^2;   % rho^2
rho =@(s,x) sqrt(rho2(s,x));           % rho
jac =@(s) sqrt(tx(s).^2 + ty(s).^2);   % jacobian

% here we compute  1/pi * (r \otimes r) / \rho^2
dx=@(s,x) rx(s,x)./rho(s,x);
dy=@(s,x) ry(s,x)./rho(s,x);
D11=@(s,x) dx(s,x).*dx(s,x)/pi;
D12=@(s,x) dx(s,x).*dy(s,x)/pi;
D22=@(s,x) dy(s,x).*dy(s,x)/pi;


% given density f,  compute (f dot r)/ (r.r)
fr=@(s,x) (rx(s,x).*fx(s) + ry(s,x).*fy(s))./rho2(s,x);

% functions for STRESS  of single layer STOKES kernel at target point $x$
T11 = @(x) int2pi( @(s) D11(s,x).*fr(s,x).*jac(s) );
T12 = @(x) int2pi( @(s) D12(s,x).*fr(s,x).*jac(s) );
T22 = @(x) int2pi( @(s) D22(s,x).*fr(s,x).*jac(s) );

% Constructing target points
%st = pi/3  + 0*rand(1)*2*pi; % pick a parametrization point
st = pi/3;
xst = [gx(st); gy(st)];   % compute corresponding point
nst = [nx(st); ny(st)]/jac(st); % normal
tst = [tx(st); ty(st)]/jac(st); % tangent

xo = xst + H * nst;  % target outside domain
xi = xst - H * nst;  % target inside  domain

% rotation matrix to for n-tau coordinate system
% used to correct stress
Rot = [ nst';tst'];  

% HERE WE COMPUTE STRESS JUMPS  (exterior-interior)
if 1
J11=T11(xo)-T11(xi);
J12=T12(xo)-T12(xi);
J22=T22(xo)-T22(xi);
numericalJ=[[J11,J12];[J12,J22]]
end

% Here we decompose the Stress tensor to normal and tangent
%     we need to first define the normal and tangent laplace double layers:
% laplace tangent  (r dot tau) / rho^2
lapt=@(s,x) (rx(s,x).*tx(s) + ry(s,x).*ty(s))./rho2(s,x);
% laplace normal   (r dot n) / rho^2
lapn=@(s,x) (rx(s,x).*nx(s) + ry(s,x).*ny(s))./rho2(s,x);

% tangent stress tensor
Tt11 = @(x) int2pi( @(s) D11(s,x).*lapt(s,x));
Tt12 = @(x) int2pi( @(s) D12(s,x).*lapt(s,x));
Tt22 = @(x) int2pi( @(s) D22(s,x).*lapt(s,x));

% normal stress tensor
Tn11 = @(x) int2pi( @(s) D11(s,x).*lapn(s,x));
Tn12 = @(x) int2pi( @(s) D12(s,x).*lapn(s,x));
Tn22 = @(x) int2pi( @(s) D22(s,x).*lapn(s,x));

% compute jumps of tangent stress tensor
if 1
Jt11=Tt11(xo) - Tt11(xi);
Jt12=Tt12(xo) - Tt12(xi);
Jt22=Tt22(xo) - Tt22(xi);
Jt = [[Jt11,Jt12];[Jt12,Jt22]];
% rotate to normal-tangent coord system. In this way, on the circle
% we obtain an angle-invariant jump (due to the Dij terms)
RJt=Rot*Jt*Rot'
Jtexact = [2*prod(tst) tst(2)^2-tst(1)^2; tst(2)^2-tst(1)^2 -2*prod(tst)];
end
% this jump should have ones in off-diagonal and zeros in diagonal

% compute jumps in normal stress tensor
if 1
Jn11=Tn11(xo) - Tn11(xi);
Jn12=Tn12(xo) - Tn12(xi);
Jn22=Tn22(xo) - Tn22(xi);
Jn = [[Jn11,Jn12];[Jn12,Jn22]];
RJn=Rot*Jn*Rot'
end
% this jump should by the identity matrix
 
% Given known normal and tangent jumps, compute the jump 
% for a given density analytically:
% tangential component
t = @(s) (fx(s).*tx(s) + fy(s).*ty(s))./jac(s);
% normal component
fn = @(s) (fx(s).*nx(s) + fy(s).*ny(s))./jac(s);
analyticJ = fn(st)*[[1,0];[0,1]] + ft(st)*Rot'*[[0,1];[1,0]]*Rot;

% this is error should be linear on $H$.
jump_error_norm = norm(numericalJ - analyticJ)


% plot stress function.
if 0
n=2^7; 
stv = linspace(st-5*pi*H,st+5*pi*H);
vals = D11(stv,xo).*fr(stv,xo).*jac(stv); vals=vals(:);
figure(2),plot(vals,'-bo');
end


% unused code for accurate integration, matlab does a good enough job.
%int2pi1 = @(f) quad(f,         0, st-pi*5*H ,1e-11);
%int2pi2 = @(f) quad(f,st-pi*5*H, st+pi*5*H ,1e-11);
%int2pi3 = @(f) quad(f,st+pi*5*H, 2*pi       ,1e-11);
%int2pi  = @(f)  int2pi1(f) + int2pi2(f) + int2pi3(f);
