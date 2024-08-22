clear;
%% Load the vesicle configurations
X = zeros(256,4);

load('True_speed6000_width0.6455_FinalIC.mat')
X(:,1) = [Xic(1:end/2)-mean(Xic(1:end/2)); Xic(end/2+1:end)-mean(Xic(end/2+1:end))];

load('True_speed3750_width0.32275_FinalIC.mat')
X(:,2) = [Xic(1:end/2)-mean(Xic(1:end/2)); Xic(end/2+1:end)-mean(Xic(end/2+1:end))];

load('True_speed3000_width0.32275_FinalIC.mat')
X(:,3) = [Xic(1:end/2)-mean(Xic(1:end/2)); Xic(end/2+1:end)-mean(Xic(end/2+1:end))];

load('finalShearXclose.mat')
X(:,4) = [Xf(1:end/2,1)-mean(Xf(1:end/2,1)); Xf(end/2+1:end,1)-mean(Xf(end/2+1:end,1))];

%%
addpath ../src/
oc = curve;
N = 128;
op = poten(N);

% calculate bending
vesicle = capsules(X,[],[],1,1,0);
bendF = vesicle.tracJump(X,zeros(N,4));

% form layers around and calculate velocity there
[~,tang] = oc.diffProp(X);
% get x and y components of normal vector at each point
nx = tang(N+1:2*N,:);
ny = -tang(1:N,:);
h = vesicle.length/vesicle.N;
dlayer = [0 sqrt(h)/2 sqrt(h)];

G = op.stokesSLmatrix(vesicle);

for k = 1 : 4
  % Points where velocity is calculated involve the points on vesicle
  tracersX = zeros(2*N, 3);
  for il = 1 : 3
    tracersX(:,il) = [X(1:end/2,k)+nx(:,k)*dlayer(il);X(end/2+1:end,k)+ny(:,k)*dlayer(il)];
  end
  
  % build tracer class
  tracers.N = N;
  tracers.nv = 2;
  tracers.X = tracersX(:,2:3);

  vesicle1 = capsules(X(:,k),[],[],1,1,0);
  
  % Get the near zone
  [~,NearV2T] = vesicle1.getZone(tracers,2);

  % Calculate velocity on the layers
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
  SLP = @(X) op.exactStokesSLdiag(vesicle1,G(:,:,k),X);


  velLayers = zeros(2*N,3);
  velLayers(:,1) = G(:,:,k)*bendF(:,k);
  velLayers(:,2:3) = op.nearSingInt(vesicle1,bendF(:,k),SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);

  % pick random points in the layers
  Xlow = [interpft(X(1:end/2,k),16);interpft(X(end/2+1:end,k),16)];
  ves2 = capsules(Xlow,[],[],1,1,0);

  [~,tangLow] = oc.diffProp(Xlow);
  % get x and y components of normal vector at each point
  nxL = tangLow(17:32,:);
  nyL = -tangLow(1:16,:);
  tracersInterp = zeros(32,2);

  for il = 1 : 2
    tracersInterp(:,il) = [Xlow(1:end/2)+nxL*dlayer(il+1)+(-0.025+0.05*rand);Xlow(end/2+1:end)+nyL*dlayer(il+1)+(-0.025+0.05*rand)];
  end
  
  trac2.N = 16;
  trac2.nv = 2;
  trac2.X = tracersInterp;

  % Get the near zone
  [~,NearV2T] = vesicle1.getZone(trac2,2);
  velTracInterp = op.nearSingInt(vesicle1,bendF(:,k),SLP,[],NearV2T,kernel,kernelDirect,trac2,false,false);


  fname = ['./NearFieldTests_VesID' num2str(k) '.mat'];
  save(fname,'tracersX','velLayers','tracersInterp','velTracInterp');
  

  figure(1); clf;
  plot(X(1:end/2,k),X(end/2+1:end,k),'r','linewidth',2)
  hold on
  plot(tracersX(1:end/2,:),tracersX(end/2+1:end,:),'ko','markersize',5)
  quiver(tracersX(1:end/2,:),tracersX(end/2+1:end,:),velLayers(1:end/2,:),velLayers(end/2+1:end,:),'b')
  plot(tracersInterp(1:end/2,:),tracersInterp(end/2+1:end,:),'gs','markersize',5,'markerfacecolor','g')
  quiver(tracersInterp(1:end/2,:),tracersInterp(end/2+1:end,:),velTracInterp(1:end/2,:),velTracInterp(end/2+1:end,:),'g')
  axis equal
  pause

end



  
  

  


