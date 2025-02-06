addpath ../src/
oc = curve;

chanWidth = 30;
X0 = oc.initConfig(32,'ellipse');
[~,area0,len0] = oc.geomProp(X0);
scale = 1/len0;

% Gives 50% area-fraction
% sx = [0.125:0.2:chanWidth-0.075]';
% sy = [0.225:0.5:chanWidth-0.225]';

% Gives 30% area-fraction
sx = [0.125:0.3:chanWidth-0.075]';
sy = [0.225:0.6:chanWidth-0.225]';

[cenx, ceny] = meshgrid(sx,sy);
cenx = cenx(:)';
ceny = ceny(:)';
nv = numel(cenx); 
angle = -ones(nv,1)*pi/2;

X = oc.initConfig(32,'nv',nv,...
  'reducedArea',0.65,...
  'angle',angle,...
  'center',[cenx;ceny], 'scale',scale);

Vsize = chanWidth; speed = 200;
[xx,yy] = meshgrid(linspace(-1,1.5+chanWidth,50)',linspace(-1,1.5+chanWidth,50)');
uu = speed*sin(xx/Vsize*pi).*cos(yy/Vsize*pi); 
vv = -speed*cos(xx/Vsize*pi).*sin(yy/Vsize*pi);

figure(1); clf;
plot(X(1:end/2,:),X(end/2+1:end,:),'r','linewidth',2)
hold on
axis equal
l = streamslice(xx,yy,uu,vv);
set(l,'Color',[12/255,44/255,132/255, 1])
set(l,'linewidth',3)
disp(['Number of vesicles:' num2str(nv)])
[~,area,len] = oc.geomProp(X);
areaFrac = area(1)*nv/chanWidth^2 * 100;
disp(['Area fraction:' num2str(areaFrac)])



