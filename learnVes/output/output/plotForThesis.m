clear all;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')


% NV = 5 -- DNN

load nv5DNNvelocity.mat

idxTsteps = 5;

radius = [1.1 2.1];
Nradii = 96;
Ntheta = 255;
thetapos = [0:Ntheta-1]'/Ntheta*2*pi;
thetapos = [thetapos;2*pi];
rpos = linspace(radius(1),radius(2),Nradii)';
[ttheta,rrpos] = meshgrid(thetapos,rpos);
      
xwalls = [interpft(Xwalls(1:end/2,:),256);Xwalls(1,:)];
ywalls = [interpft(Xwalls(end/2+1:end,:),256);Xwalls(end/2+1,:)];

speed = 100;
xlimits = [min(min(Xwalls(1:end/2,:))) max(max(Xwalls(1:end/2,:)))];
ylimits = [min(min(Xwalls(1+end/2:end,:))) max(max(Xwalls(1+end/2:end,:)))];
titY = 2.4; titX = -0.5;

vxCouette = 100*(ytracers/3.84).*(1-4.84./(xtracers.^2+ytracers.^2));
vyCouette = 100*(-xtracers/3.84).*(1-4.84./(xtracers.^2+ytracers.^2));

% put velocities together
velxTraAll{1} = velxTra;
velyTraAll{1} = velyTra;
XstoreAll{1} = Xstore;

% NV = 10 -- DNN
load nv10DNNvelocity
velxTraAll{2} = velxTra;
velyTraAll{2} = velyTra;
XstoreAll{2} = Xstore;

% NV = 5 -- TRUTH
load nv5Truevelocity
velxTraAll{3} = velxTra;
velyTraAll{3} = velyTra;
XstoreAll{3} = Xstore;

% NV = 10 -- TRUTH
load nv10N64Truevelocity
velxTraAll{4} = velxTra;
velyTraAll{4} = velyTra;
XstoreAll{4} = Xstore;


xc = [-2.4; 2.4];
yc = 0;

colors = [ 1 0 0; 0.5 0.5 0.5];
figure(1);clf;
for icase = 1 : 2
  velxTra = velxTraAll{icase};
  velyTra = velyTraAll{icase};
  Xstore = XstoreAll{icase};
  
  velxTraTrue = velxTraAll{icase+2};
  velyTraTrue = velyTraAll{icase+2};
  XstoreTrue = XstoreAll{icase+2};
  
  for it = 1 : numel(idxTsteps)
    idx = idxTsteps(it);
    time = (timeSteps(idx)-1)*1e-4;
    figure(1);hold on;
    plot(xwalls+xc(icase),ywalls,'Color',[0 0 0],'linewidth',1)
%     plot(cos(speed*time)+xc(icase),sin(speed*time)+yc,'o','Color',[0 0 0],...
%       'markersize',6,'markerfacecolor',[0 0 0])

    x = interpft(XstoreTrue(1:end/2,:,idx),256);
    y = interpft(XstoreTrue(end/2+1:end,:,idx),256);
    vecx = [x;x(1,:)]+xc(icase); vecy = [y;y(1,:)]+yc;

    plot(vecx,vecy,'Color',[0.5 0.5 0.5],'linewidth',2)
    hVes = fill(vecx,vecy,[0.5 0.5 0.5]);
    set(hVes,'edgecolor',[0.5 0.5 0.5]);
    
    x = interpft(Xstore(1:end/2,:,idx),256);
    y = interpft(Xstore(end/2+1:end,:,idx),256);
    vecx = [x;x(1,:)]+xc(icase); vecy = [y;y(1,:)]+yc;

    plot(vecx,vecy,'r','linewidth',2)
    hVes = fill(vecx,vecy,'r');
    set(hVes,'edgecolor','r');
    
    
    % find the mean radial position
    cx = mean(x); cy = mean(y); cr = sqrt(cx.^2+cy.^2); meanCr = mean(cr);
    x = interpft(XstoreTrue(1:end/2,:,idx),256);
    y = interpft(XstoreTrue(end/2+1:end,:,idx),256);
    cx = mean(x); cy = mean(y); cr = sqrt(cx.^2+cy.^2); meanCrTrue = mean(cr);
    
    thet = [0:255]'/256*2*pi;
    plot(meanCrTrue*cos(thet)+xc(icase),meanCrTrue*sin(thet)+yc,'--','Color',[.5 .5 .5]);
    plot(meanCr*cos(thet)+xc(icase),meanCr*sin(thet)+yc,'r--');

  end
end
% cb = colorbar('TickLabelInterpreter','latex');
% cb.Ticks = [4;8;12;16;20;24];
% cb.Label.Interpreter = 'latex';
axis equal
xlim([-4.8 4.8])
ylim([-2.4 2.4])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box off
set(gca,'visible','off')

    
  