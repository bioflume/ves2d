clear all;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')


% NV = 5 -- DNN

load nv5DNNvelocity.mat

idxTsteps = [1;3;5];

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


xc = [-7.2; -2.4; 2.4; 7.2; 12];
yc = [2.4; -2.4];

colors = [ 1 0 0; 0 0 0; 1 0 0; 0 0 0];
figure(1);clf;
for icase = 1 : 2
  velxTra = velxTraAll{icase};
  velyTra = velyTraAll{icase};
  Xstore = XstoreAll{icase};
  
  velxTraTrue = velxTraAll{icase+2};
  velyTraTrue = velyTraAll{icase+2};
  XstoreTrue = XstoreAll{icase+2};
  
  for it = 1 : 3
    idx = idxTsteps(it);
    time = (timeSteps(idx)-1)*1e-4;
    figure(1);hold on;
    plot(xwalls+xc(it),ywalls+yc(icase),'Color',[.5 .5 .5],'linewidth',1)
    plot(cos(speed*time)+xc(it),sin(speed*time)+yc(icase),'o','Color',[.5 .5 .5],...
      'markersize',6,'markerfacecolor',[.5 .5 .5])

    x = interpft(XstoreTrue(1:end/2,:,idx),256);
    y = interpft(XstoreTrue(end/2+1:end,:,idx),256);
    vecx = [x;x(1,:)]+xc(it); vecy = [y;y(1,:)]+yc(icase);

    plot(vecx,vecy,'k','linewidth',2)
    hVes = fill(vecx,vecy,'k');
    set(hVes,'edgecolor','k');
    
    x = interpft(Xstore(1:end/2,:,idx),256);
    y = interpft(Xstore(end/2+1:end,:,idx),256);
    vecx = [x;x(1,:)]+xc(it); vecy = [y;y(1,:)]+yc(icase);

    plot(vecx,vecy,'r','linewidth',2)
    hVes = fill(vecx,vecy,'r');
    set(hVes,'edgecolor','r');
    
    if it == 3    
    % find the mean radial position
    cx = mean(x); cy = mean(y); cr = sqrt(cx.^2+cy.^2); meanCr = mean(cr);
    x = interpft(XstoreTrue(1:end/2,:,idx),256);
    y = interpft(XstoreTrue(end/2+1:end,:,idx),256);
    cx = mean(x); cy = mean(y); cr = sqrt(cx.^2+cy.^2); meanCrTrue = mean(cr);
    
    thet = [0:255]'/256*2*pi;
    plot(meanCrTrue*cos(thet)+xc(it),meanCrTrue*sin(thet)+yc(icase),'k--');
    plot(meanCr*cos(thet)+xc(it),meanCr*sin(thet)+yc(icase),'r--');
    
    % First plot velocity for NNA   
    plot(xwalls+xc(it+1),ywalls+yc(icase),'Color',[.5 .5 .5],'linewidth',1)
    plot(cos(speed*time)+xc(it+1),sin(speed*time)+yc(icase),'o','Color',[.5 .5 .5],...
      'markersize',6,'markerfacecolor',[.5 .5 .5])
  
    x = interpft(Xstore(1:end/2,:,idx),256);
    y = interpft(Xstore(end/2+1:end,:,idx),256);
    vecx = [x;x(1,:)]+xc(it+1); vecy = [y;y(1,:)]+yc(icase);
    
    velx = velxTra(:,:,idx)-vxCouette; 
    vely = velyTra(:,:,idx)-vyCouette;
    
    [C,h] = contour(xtracers+xc(it+1),ytracers+yc(icase),sqrt(velx.^2+vely.^2));
    colormap(cool)
    h.LineWidth = 1.5;
    caxis([0.01 24])
    h.LevelStep = 2;
    
    plot(vecx,vecy,'r','linewidth',2)
    hVes = fill(vecx,vecy,'r');
    set(hVes,'edgecolor','r');
    
    % Then for True
    plot(xwalls+xc(it+2),ywalls+yc(icase),'Color',[.5 .5 .5],'linewidth',1)
    plot(cos(speed*time)+xc(it+2),sin(speed*time)+yc(icase),'o','Color',[.5 .5 .5],...
      'markersize',6,'markerfacecolor',[.5 .5 .5])
  
    x = interpft(XstoreTrue(1:end/2,:,idx),256);
    y = interpft(XstoreTrue(end/2+1:end,:,idx),256);
    vecx = [x;x(1,:)]+xc(it+2); vecy = [y;y(1,:)]+yc(icase);
    
    velx = velxTraTrue(:,:,idx)-vxCouette; 
    vely = velyTraTrue(:,:,idx)-vyCouette;
    
    [C,h] = contour(xtracers+xc(it+2),ytracers+yc(icase),sqrt(velx.^2+vely.^2));
    colormap(cool)
    h.LineWidth = 1.5;
    caxis([0.01 24])
    h.LevelStep = 2;
    
    plot(vecx,vecy,'k','linewidth',2)
    hVes = fill(vecx,vecy,'k');
    set(hVes,'edgecolor','k');
    end

  end
end
% cb = colorbar('TickLabelInterpreter','latex');
% cb.Ticks = [4;8;12;16;20;24];
% cb.Label.Interpreter = 'latex';
axis equal
xlim([-9.6 14.4])
ylim([-4.8 4.8])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box off
set(gca,'visible','off')

    
  