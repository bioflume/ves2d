addpath ../src/
oc = curve;

N = 32;

X0 = oc.initConfig(N,'ellipse');
[~,area0,len0] = oc.geomProp(X0);
X0 = X0/len0;
scale = 1/len0;

H = max(X0(1:end/2))-min(X0(1:end/2));
W = max(X0(end/2+1:end))-min(X0(end/2+1:end)); 

l = W*125 + 124*0.06;
h = 0.2 + 0.07 * 7 + 8 * H;

h/2

cy = linspace(-h/2+0.1 + H/2, h/2-0.1-H/2,8)';
cx = linspace(W/2, l-W/2,125)';

[ccx, ccy] = meshgrid(cx,cy);
ccx = ccx + (-0.01 + 0.02 * rand(8,125));
ccy = ccy + (-0.015 + 0.03 * rand(8,125));

ccx = ccx(:)';
ccy = ccy(:)';
nv = numel(ccx); 
angle = -ones(nv,1)*pi/2;

X = oc.initConfig(N,'nv',nv,...
  'reducedArea',0.65,...
  'angle',angle,...
  'center',[ccx;ccy], 'scale',scale);

% figure(1);clf;
% plot([X(1:end/2,:);X(1,:)],[X(end/2+1:end,:);X(end/2+1,:)],'r','linewidth',2)
% axis equal
% 
% save 1000vesShapeEllips X

%%

t = (0:N-1)'*2*pi/N;
r = 1 + .7*cos(2*t);
x = r.*cos(t); x = 2*x/max(x);
y = r.*sin(t); y = y/max(y); 
X8 = [x;y];
[ra,area,len] = oc.geomProp(X8);

X8 = X8/len;
[ra,area,len] = oc.geomProp(X8);
H = max(X8(1:end/2))-min(X8(1:end/2));
W = max(X8(end/2+1:end))-min(X8(end/2+1:end)); 


l = W*125 + 124*0.06;
h = 0.2 + 0.07 * 7 + 8 * H;

h/2 

cy = linspace(-h/2+0.1 + H/2, h/2-0.1-H/2,8)';
cx = linspace(W/2, l-W/2,125)';

[ccx, ccy] = meshgrid(cx,cy);
ccx = ccx + (-0.01 + 0.02 * rand(8,125));
ccy = ccy + (-0.015 + 0.03 * rand(8,125));

X0 = oc.initConfig(N,'figureEight');
[~,area0,len0] = oc.geomProp(X0);
X0 = X0/len0;

ccx = ccx(:)';
ccy = ccy(:)';
nv = numel(ccx); 
theta = -ones(nv,1)*pi/2;

 X = zeros(2*N,nv);
  for k=1:nv
    X(1:N,k) = cos(theta(k)) * X0(1:N) - ...
      sin(theta(k)) * X0(N+1:2*N) + ccx(k);
    X(N+1:2*N,k) = sin(theta(k)) * X0(1:N) +  ...
      cos(theta(k)) * X0(N+1:2*N) + ccy(k);
  end

% figure(1);clf;
% plot([X(1:end/2,:);X(1,:)],[X(end/2+1:end,:);X(end/2+1,:)],'r','linewidth',2)
% axis equal
% 
% save 1000vesShape8 X