clear all;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

load nv5DNNvelocity.mat

time = (timeSteps-1)*1e-4;

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

if 1
vorticity = zeros(Nradii-1,Ntheta+1,numel(timeSteps));
for k = 1 : numel(timeSteps)
  vradial = (velxTra(:,:,k)-vxCouette).*cos(ttheta)+(velyTra(:,:,k)-vyCouette).*sin(ttheta);
  vangul = (velyTra(:,:,k)-vyCouette).*cos(ttheta)-(velxTra(:,:,k)-vxCouette).*sin(ttheta);
  
  radDer = zeros(Nradii-1,Ntheta);
  for ith = 1 : Ntheta
    radDer(:,ith) = (rrpos(2:end,ith).*vangul(2:end,ith)-...
        rrpos(1:end-1,ith).*vangul(1:end-1,ith))./...
        (rrpos(2:end,ith)-rrpos(1:end-1,ith));
    radDer(:,ith) = radDer(:,ith)./rrpos(2:end,ith);
  end
  
  angDer = zeros(Nradii-1,Ntheta);
  for ith = 1 : Nradii-1
    angDer(ith,:) = 1./rrpos(ith,2:end).*(vradial(ith,2:end)-vradial(ith,1:end-1))./...
        (ttheta(ith,2:end)-ttheta(ith,1:end-1));
  end
  vorticity(:,1:end-1,k) = radDer-angDer;
  vorticity(:,end,k) = vorticity(:,1,k);
end
else
    
vorticity = zeros(Nradii-1,Ntheta+1,numel(timeSteps));
for k = 1 : numel(timeSteps)
  vradial = (velxTra(:,:,k).*cos(ttheta)+velyTra(:,:,k).*sin(ttheta));
  vangul = velyTra(:,:,k).*cos(ttheta)-velxTra(:,:,k).*sin(ttheta);
  radDer = zeros(Nradii-1,Ntheta);
  for ith = 1 : Ntheta
    val = rrpos(:,ith).*vangul(:,ith);
    for ir = 1 : Nradii-2
      radDer(ir,ith) = (-1/2*val(ir+2)+2*val(ir+1)-...
          3/2*val(ir))./(rrpos(ir+1,ith)-rrpos(ir,ith));
    end
    radDer(Nradii-1,ith) = (val(Nradii)-val(Nradii-1))./...
        (rrpos(Nradii,ith)-rrpos(Nradii-1,ith));
  end
  
  angDer = zeros(Nradii-1,Ntheta);
  for ith = 1 : Nradii-1
    for ir = 1 : Ntheta-1
      angDer(ith,ir) = 1./rrpos(ith,ir).*(-1/2*vradial(ith,ir+2)+...
          2*vradial(ith,ir+1)-3/2*vradial(ith,ir))./...
          (ttheta(ith,ir+1)-ttheta(ith,ir));
    end
    angDer(ith,Ntheta) = 1./rrpos(ith,Ntheta).*(vradial(ith,Ntheta+1)-...
        vradial(ith,Ntheta))./(ttheta(ith,Ntheta+1)-ttheta(ith,Ntheta));
  end
  vorticity(:,1:end-1,k) = radDer-angDer;
  vorticity(:,end,k) = vorticity(:,1,k);
end    
    
    
    
end

for k = 1 : numel(timeSteps)
  figure;hold on;  
  plot(xwalls,ywalls,'k','linewidth',1.5)
  plot(cos(speed*time(k)),sin(speed*time(k)),'ko',...
    'markersize',8,'markerfacecolor','k')

  x = interpft(Xstore(1:end/2,:,k),256);
  y = interpft(Xstore(end/2+1:end,:,k),256);
  vecx = [x;x(1,:)]; vecy = [y;y(1,:)];
  
  plot(vecx,vecy,'r','linewidth',2)
  hVes = fill(vecx,vecy,'r');
  set(hVes,'edgecolor','r')
  
  if 1
  velx = velxTra(:,:,k)-vxCouette; vely = velyTra(:,:,k)-vyCouette;
  [C,h] = contour(xtracers,ytracers,sqrt(velx.^2+vely.^2));
  h.LineWidth = 2;
  caxis([2 20])
  else
  [C,h] = contour(xtracers(2:end,1:end),ytracers(2:end,1:end),...
      vorticity(:,:,k),200);
  h.LineWidth = 2;
  h.LineColor = [0 0 1];
  h.LevelStep = 10;
  h.LevelStepMode = 'manual';
  end
%   quiver(xtracers,ytracers,velx,vely)
  
%   starty = zeros(10,1);
%   startx = linspace(1.1,2.1,10)';
%   h = streamline(xtracers,ytracers,velx,vely,startx,starty);

  plot(vecx,vecy,'r','linewidth',2)
  hVes = fill(vecx,vecy,'r');
  set(hVes,'edgecolor','r')
  
  axis equal
  xlim(xlimits)
  ylim(ylimits)
  title(time(k))
  
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'ztick',[]);

  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  set(gca,'zcolor','w');
  box off
  set(gca,'visible','off')
  titleStr = ['t = ' num2str(time(k),'%.2f')];
  text(titX,titY,titleStr,'FontSize',28,'FontName','Palatino')

end