clear; clc;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


addpath ../../src/


fileName = 'parabolic_shape8_nv2000_auglag.bin';
[vesx, vesy, time, NN, nv, xinit, yinit] = loadShanVesFile(fileName);

% oc = curve;
ids = [1 1000 1500 2500 3500];
rpos = [[0 -1.25 2 1.5];...
    [19 2 2 1.5];...
    [47 1.5 2 1.5];...
    [47 -4.5 2 1.5];...
    [64.5 -4.5 2 1.5]];
for k = 1 : numel(ids)
figure(2*(k-1)+1); clf;
id = ids(k);
xvec = [interpft(vesx(:,:,id),32);vesx(1,:,id)];
yvec = [interpft(vesy(:,:,id),32);vesy(1,:,id)];
plot(xvec, yvec, 'Color',[202,0,32]/255,'linewidth',1)
hold on
plot(-0.99, -4.99, 'k.','markersize',0.001)
plot(79.99, 4.99, 'k.','markersize',0.001)

% hFill = fill(xvec, yvec, [202,0,32]/255);
% set(hFill,'EdgeColor', [202,0,32]/255);
r = rectangle('Position',rpos(k,:)');
r.LineWidth = 1;
r.LineStyle = '-';
r.EdgeColor = [0 0 0];
axis equal
ylim([-5 5])
xlim([-1 80])


set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')

ax = gca;
exportgraphics(ax,['~/Desktop/parabolAll_k' num2str(k) '.png'],'Resolution',300)

figure(2*k); clf;
xvec = [interpft(vesx(:,:,id),128);vesx(1,:,id)];
yvec = [interpft(vesy(:,:,id),128);vesy(1,:,id)];
plot(xvec, yvec, 'Color',[202,0,32]/255,'linewidth',1)
hold on
hFill = fill(xvec, yvec, [202,0,32]/255);
set(hFill,'EdgeColor', [202,0,32]/255);
plot(rpos(k,1), rpos(k,2), 'k.','markersize',0.001)
plot(rpos(k,1)+rpos(k,3), rpos(k,2)+rpos(k,4), 'k.','markersize',0.001)
axis equal
xlim([rpos(k,1) rpos(k,1)+rpos(k,3)])
ylim([rpos(k,2) rpos(k,2)+rpos(k,4)])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')

ax = gca;
exportgraphics(ax,['~/Desktop/parabolZoom_k' num2str(k) '.png'],'Resolution',300)

end
% 




