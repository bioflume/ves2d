clear; clc;
tube_length = 20;
tube_height = 6;

a = tube_length/2; 
b = tube_height/2;
order = 8;
N = 256;
%% Testing the height
cvals = (0:0.1:3)';

for k = 1 : numel(cvals)

c = cvals(k);

Nsides = ceil(0.5*b/(2*a+2*b)*N);
Ntop = (N-4*Nsides)/2;

t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
t = [t1;t2;t3;t4;t5];

r = (cos(t).^order + sin(t).^order).^(-1/order);
x = a*r.*cos(t); y = b*r.*sin(t);
ind = abs(x) < pi & y < 0;
y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);
X0 = [x;y];

miny = max(y(ind));
ind = abs(x) < pi & y > 0;
maxy = max(y(ind));
Hvals(k,1) = maxy-miny;

figure(1); clf;
plot(x,y,'linewidth',2)
axis equal
title(['c = ' num2str(cvals(k)) ', H = ' num2str(Hvals(k,1))])
ax = gca;
exportgraphics(ax,['~/Desktop/stenosisFigs/cval' num2str(c) '.png'],'Resolution', 300)

pause(0.1)


end

%% Testing the width of the 
wvals = (1:0.5:5)';
c = 0.5;
for k = 1 : numel(wvals)

w = wvals(k);

Nsides = ceil(0.5*b/(2*a+2*b)*N);
Ntop = (N-4*Nsides)/2;

t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
t = [t1;t2;t3;t4;t5];

r = (cos(t).^order + sin(t).^order).^(-1/order);
x = a*r.*cos(t); y = b*r.*sin(t);
ind = abs(x)/w * pi < pi & y < 0;
y(ind) = y(ind).*(1-c*cos(x(ind)/w * pi))/(1+c);
X0 = [x;y];

miny = max(y(ind));
ind = abs(x) < pi & y > 0;
maxy = max(y(ind));
Hvals(k,1) = maxy-miny;

figure(1); clf;
plot(x,y,'linewidth',2)
axis equal
title(['w = ' num2str(wvals(k)) ', H = ' num2str(Hvals(k,1))])
ax = gca;
exportgraphics(ax,['~/Desktop/stenosisFigs/wval' num2str(w) '.png'],'Resolution', 300)

pause(0.1)


end

%% With the physical inputs
% Tube is centered at (0, 0)
tube_length = 20;
tube_height = 6;
min_height = 3;
contraction_width = 2;
N = 256; % Number of discretization points

a = tube_length/2; 
b = tube_height/2;
c = tube_height/min_height - 1;
order = 8;
Nsides = ceil(0.5*b/(2*a+2*b)*N);
Ntop = (N-4*Nsides)/2;

t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
t = [t1;t2;t3;t4;t5];

r = (cos(t).^order + sin(t).^order).^(-1/order);
x = a*r.*cos(t); y = b*r.*sin(t);
ind = abs(x)/contraction_width * pi < pi & y < 0;
y(ind) = y(ind).*(1-c*cos(x(ind)/contraction_width * pi))/(1+c);
X0 = [x;y];

miny = max(y(ind));
ind = abs(x)/contraction_width * pi < pi & y > 0;
maxy = max(y(ind));

Hval = maxy-miny;
figure(1); clf;
plot(x,y,'linewidth',2)
axis equal
title(['H = ' num2str(Hval)])
