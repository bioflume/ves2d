%check the jumps for the pressure of the 2D stokes DLP

clear all; clear globals; clf
H = 1e-3; %Distance from the boundary

st = 20*pi/256;
%syms sf
%fx = @(s) cos(0*s);
%fy = @(s) 0*cos(0*s);
fx = @(s) cos(s);
fy = @(s) cos(s);
%traction jump
Dfx = @(s) -sin(s);
Dfy = @(s) -sin(s);
if 0
  fx = cos(2*sf) + sin(3*sf);
  fy = 0.2*cos(sf) + sin(4*sf);
end
%Some random periodic traction jump
%Dfx = diff(fx,sf);
%Dfy = diff(fy,sf);
%Need the derivative of the traction jump to compute
%the jump in pressure

%fx = matlabFunction(fx);
%fy = matlabFunction(fy);
%Dfx = matlabFunction(Dfx);
%Dfy = matlabFunction(Dfy);

int2pi = @(f) quad(f,0,2*pi,1e-6);
%short-hand for numerical integration over [0,2*pi]

syms sS
%Define a symbolic variable
sR = 1 + 3*0.1*sin(8*sS);
%sR = 1 + 0.2*cos(2*sS) + 0.3*sin(8*sS);
%sR = 1;
%radius of the boundary
sX = sR.*cos(sS);
sY = sR.*sin(sS);
%Symbolic boundary of domain
sDX = diff(sX,sS);
sDY = diff(sY,sS);
%symbolic derivative of boundary
jacobian = sqrt(sDX.^2+sDY.^2);
%sDX = sDX./jacobian;
%sDY = sDY./jacobian;

gx = matlabFunction(sX);
gy = matlabFunction(sY);
tx = matlabFunction(sDX);
ty = matlabFunction(sDY);
%Unormalized tangent
%Convert symbolic expressions to numeric ones

N = 1024;
theta = (0:N-1)*2*pi/N;
plot(gx(theta),gy(theta))
axis equal

rx = @(s,x) x(1)-gx(s);
ry = @(s,x) x(2)-gy(s);
rho2 = @(s,x) rx(s,x).^2 + ry(s,x).^2;
rho = @(s,x) sqrt(rho2(s,x));
jac = @(s) sqrt(tx(s).^2 + ty(s).^2);
nx = @(s) ty(s)./jac(s);
ny = @(s) -tx(s)./jac(s);
%Compute handles for jacobian, normal, tangent, rx, ry

rof11 = @(s,x) rx(s,x).*fx(s);
rof12 = @(s,x) ry(s,x).*fx(s);
rof21 = @(s,x) rx(s,x).*fy(s);
rof22 = @(s,x) ry(s,x).*fy(s);
nor11 = @(s,x) nx(s).*rx(s,x);
nor12 = @(s,x) ny(s).*rx(s,x);
nor21 = @(s,x) nx(s).*ry(s,x);
nor22 = @(s,x) ny(s).*ry(s,x);
ror11 = @(s,x) rx(s,x).*rx(s,x);
ror12 = @(s,x) ry(s,x).*rx(s,x);
ror21 = @(s,x) rx(s,x).*ry(s,x);
ror22 = @(s,x) ry(s,x).*ry(s,x);
fon11 = @(s) fx(s).*nx(s);
fon12 = @(s) fy(s).*nx(s);
fon21 = @(s) fx(s).*ny(s);
fon22 = @(s) fy(s).*ny(s);
%All the outer products.  Including the 1,2 and 2,1 entries
%so that we can make sure the tensor is symmetric

ndotf = @(s) nx(s).*fx(s) + ny(s).*fy(s);
rdotn = @(s,x) rx(s,x).*nx(s) + ry(s,x).*ny(s);
rdotf = @(s,x) rx(s,x).*fx(s) + ry(s,x).*fy(s);
%All the dot products


T11 = @(s,x) 1*fon11(s)./rho2(s,x) - ...
    1*8*rdotn(s,x).*rdotf(s,x)./rho2(s,x).^3.*ror11(s,x) +  ...
    1*rdotn(s,x)./rho2(s,x).^2.*rof11(s,x) + ...
    1*rdotn(s,x).*rdotf(s,x)./rho2(s,x).^2 + ...
    1*ndotf(s)./rho2(s,x).^2.*ror11(s,x) + ...
    1*rdotf(s,x)./rho2(s,x).^2.*nor11(s,x);

T12 = @(s,x) 1*fon12(s)./rho2(s,x) - ...
    8*rdotn(s,x).*rdotf(s,x)./rho2(s,x).^3.*ror12(s,x) +  ...
    rdotn(s,x)./rho2(s,x).^2.*rof12(s,x) + ...
    ndotf(s)./rho2(s,x).^2.*ror12(s,x) + ...
    rdotf(s,x)./rho2(s,x).^2.*nor12(s,x);
% (1,2) and (2,1) components don't have the first term
% as it comes from the pressure which is multiplied by
% the identity matrix
T21 = @(s,x) 1*fon21(s)./rho2(s,x) - ...
    8*rdotn(s,x).*rdotf(s,x)./rho2(s,x).^3.*ror21(s,x) +  ...
    rdotn(s,x)./rho2(s,x).^2.*rof21(s,x) + ...
    ndotf(s)./rho2(s,x).^2.*ror21(s,x) + ...
    rdotf(s,x)./rho2(s,x).^2.*nor21(s,x);

T22 = @(s,x) 1*fon22(s)./rho2(s,x) - ...
    8*rdotn(s,x).*rdotf(s,x)./rho2(s,x).^3.*ror22(s,x) +  ...
    rdotn(s,x)./rho2(s,x).^2.*rof22(s,x) + ...
    rdotn(s,x).*rdotf(s,x)./rho2(s,x).^2 + ...
    ndotf(s)./rho2(s,x).^2.*ror22(s,x) + ...
    rdotf(s,x)./rho2(s,x).^2.*nor22(s,x);

Ten11 = @(x) int2pi(@(s) 1/pi*T11(s,x).*jac(s));
Ten12 = @(x) int2pi(@(s) 1/pi*T12(s,x).*jac(s));
Ten21 = @(x) int2pi(@(s) 1/pi*T21(s,x).*jac(s));
Ten22 = @(x) int2pi(@(s) 1/pi*T22(s,x).*jac(s));


%st = 0.09321711;
xst = [gx(st);gy(st)];
nst = [nx(st);ny(st)];

xo = xst + H*nst;
xi = xst - H*nst;

%J11 = Ten11(xo) - Ten11(xi);
%J12 = Ten12(xo) - Ten12(xi);
%J21 = Ten21(xo) - Ten21(xi);
%J22 = Ten22(xo) - Ten22(xi);

%J = [J11 J12; J21 J22];
t1 = tx(st)/jac(st); t2 = ty(st)/jac(st);

Jump = 2/jac(st)*(Dfx(st)*t1+Dfy(st)*t2) * ...
  (eye(2) + [t1^2-t2^2 2*t1*t2; ...
      2*t1*t2 -t1^2+t2^2]);

%disp(J)
%disp(Jump)




