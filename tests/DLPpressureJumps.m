%check the jumps for the pressure of the 2D stokes DLP

clear all; clear globals; clf
H = 1e-5; %Distance from the boundary

syms sf
%fx = @(s) 1;
%fy = @(s) 1;
fx = @(s) cos(s);
fy = @(s) cos(s);
Dfx = @(s) -sin(s);
Dfy = @(s) -sin(s);
%fx = 0*cos(sf)+1;
%fy = 0*cos(sf)+1;
%fx = cos(sf);
%fy = sin(sf);
%traction jump
if 0
  fx = cos(2*sf) + sin(3*sf);
  fy = 0.2*cos(sf) + sin(4*sf);
end
%Some random periodic traction jump
%Dfx = diff(fx,sf);
%Dfy = diff(fy,sf);
%Need derivative for of the traction jump to compute 
%the jump in pressure

%fx = matlabFunction(fx);
%fy = matlabFunction(fy);
%Dfx = matlabFunction(Dfx);
%Dfy = matlabFunction(Dfy);

int2pi = @(f) quad(f,0,2*pi,1e-6);
%short-hand for numerical integration over [0,2*pi]

syms sS
%Define a symbolic variable
%sR = 1 + 3*0.1*sin(8*sS);
%sR = 1 + 0.2*cos(2*sS) + 0.3*sin(8*sS);
%sR = 1;
%radius of the boundary
%sX = sR.*cos(sS);
%sY = sR.*sin(sS);
sX = 2*cos(sS);
sY = sin(sS);
%Symbolic boundary of domain
sDX = diff(sX,sS);
sDY = diff(sY,sS);
%symbolic derivative of boundary

gx = matlabFunction(sX);
gy = matlabFunction(sY);
tx = matlabFunction(sDX);
ty = matlabFunction(sDY);
%Unormalized tangent
%Convert symbolic expressions to numeric ones

N = 256;
theta = (0:N-1)*2*pi/N;
plot(gx(theta),gy(theta))
axis equal

rx = @(s,x) x(1)-gx(s);
ry = @(s,x) x(2)-gy(s);
nx = @(s) ty(s);
ny = @(s) -tx(s);
rho2 = @(s,x) rx(s,x).^2 + ry(s,x).^2;
rho = @(s,x) sqrt(rho2(s,x));
jac = @(s) sqrt(tx(s).^2 + ty(s).^2);
%Compute handles for jacobian, normal, tangent, rx, ry

st = 5.2383;
xst = [gx(st);gy(st)];
%Target point
tst = [tx(st);ty(st)]/jac(st);
nst = [nx(st);ny(st)]/jac(st);
%Unit tangent and normal

dx = @(s,x) rx(s,x)./rho(s,x);
dy = @(s,x) ry(s,x)./rho(s,x);
D11 = @(s,x) dx(s,x).*dx(s,x);
D12 = @(s,x) dx(s,x).*dy(s,x);
D22 = @(s,x) dy(s,x).*dy(s,x);

p = @(s,x) 1/pi*1./rho2(s,x).*...
    ((nx(s) - 2*D11(s,x).*nx(s) - 2*D12(s,x).*ny(s)).*fx(s) + ...
    (ny(s) - 2*D12(s,x).*nx(s) - 2*D22(s,x).*ny(s)).*fy(s)); 
p2 = @(s,x) 1/pi*1./rho2(s,x).*...
    ((nx(s) - 2*D11(s,x).*nx(s) - 2*D12(s,x).*ny(s)).*(fx(s)-fx(st)) + ...
    (ny(s) - 2*D12(s,x).*nx(s) - 2*D22(s,x).*ny(s)).*(fy(s)-fy(st))); 
%Don't need to normalize normal vector because we need to multiply by
%the jacobian when computing the integral
%pnor = @(s,x) 1/pi*1./rho2(s,x).*(fx(s).*nx(s) + fy(s).*ny(s)) - ...
%    2/pi./rho2(s,x).^2.*(rx(s,x).*nx(s)+ry(s,x).*ny(s)).^2.*...
%    (fx(s).*nx(s) + fy(s).*ny(s))./jac(s).^2;
%ptan = @(s,x) -2/pi./rho2(s,x).^2.*(rx(s,x).*nx(s)+ry(s,x).*ny(s)).*(rx(s,x).*tx(s)+ry(s,x).*ty(s)).*...
%    (fx(s).*tx(s) + fy(s).*ty(s))./jac(s).^2;
%Need to normalize two of the three normal and tangential vectors


xo = xst + H*nst;
xi = xst - H*nst;
%xi = 0*xi;
%Points slightly outside and inside the domain

press = @(x) int2pi(@(s) p(s,x));
%pressnor = @(x) int2pi(@(s) pnor(s,x));
%presstan = @(x) int2pi(@(s) ptan(s,x));

J = press(xo) - press(xi);
%Jnor = pressnor(xo) - pressnor(xi);
%Jtan = presstan(xo) - presstan(xi);

%disp(['Jump in pressure is ' num2str(J)]);
%disp(['Jump in normal component of pressure is ' num2str(Jnor)]);
%disp(['Jump in tangential component of pressure is ' num2str(Jtan)]);

Jexact = 2/jac(st)*(Dfx(st)*tst(1) + Dfy(st)*tst(2));
%err = abs(J-Jexact);
%disp(['Error in jump is ' num2str(err)]);




