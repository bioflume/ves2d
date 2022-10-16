clear all; clc

n = 64;                                 
nv = 192; %%55
prams.T = 10;                           
prams.ts = 1/25;
prams.m = prams.T/prams.ts;                          
prams.kappa = 1e-1;
prams.rhoOut = 1;                       
prams.order = 1;                        
prams.ts = prams.T/prams.m;             
prams.flowType = 'confined';            
prams.M = 64*[6 4];                     
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette');
prams.bc = @(x) forcing(x,'couette');  

options.useGPU = 1;
prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams,'direct',options.useGPU);

options.usePlot = 1; 
options.progressBar = 1;
options.AxesHandle = gca;
options.saveData = 1;                  
options.dataStride = floor(prams.m/20);
options.fileName = '../results/clusteringFull.mat';

% %% Checking the field velocity
% % domain = fixedBound(prams.M,prams.bd,1);
% % [theta R] = meshgrid(linspace(0,pi/2,16),linspace(22,30,16));
% % [xx yy] = pol2cart(theta(:),R(:));
% % 
% % ur = prams.bc([xx yy]);
% % u = prams.vInf([xx;yy],[]); u = reshape(u,[],2);
% % e = u-ur;
% % disp(['local inf error: ' num2str(sqrt(max(dot(e,e,2)./dot(ur,ur,2))))]);
% %  
% % viewer(domain);
% % hold on; quiver(xx,yy,u(:,1),u(:,2));
% % axis([0 32 0 32])
% 
% % %% Generating the distribution of vesicles
% [t prams options] = monitor(prams,options);
domain = fixedBound(prams.M,prams.bd,1);
% X = boundary(n,'sizes',domain,'nv',nv);
%  
% save '../results/clusteringXFull'
%  load '../results/clusteringXFull'
%  prams.T = 10;
%  prams.ts = 1/100;
%  prams.m = prams.T/prams.ts;                           
% viewer(X,[],prams,options);
%  [XF status] = Ves2D(X,prams,options,@monitor);
%  save '../results/clusteringFullFinal'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading data files and extracting data

fileId = fopen(options.fileName,'r');
Result = fread(fileId,'double');
fclose(fileId);

Result = reshape(Result,5*nv*n+1,[]); 

if(nv==192), Result = Result(:,2:end); end %The first entry is not
                                           %relevant -- from previous saves
Time = Result(end,:);
Xv = Result(1:2*nv*n,:);

a = .9;
lineColor = a*rand(nv,3);
X0 = reshape(Xv(:,1),2*n,nv);
[ra0 a0 l0] = reducedVolume(X0);
errorA = []; errorL = [];

%% making the box
X0b = linspace(13.9,22.87,50)'; 
Xb = [[X0b 13.9*ones(size(X0b))];[22.87*ones(size(X0b)) X0b]];
X0b = flipud(X0b);
Xb = [Xb;[X0b 22.87*ones(size(X0b))];[13.9*ones(size(X0b)) X0b]];

close all;

for shot = 1:size(Xv,2)
  h1 = figure;
  set(gcf,'Units','centimeters');
  set(gcf,'position',[5 5 7.5 7.5]);

  plot(domain(1).X([1:end 1],1),domain(1).X([1:end 1],2),'k-','linewidth',2);
  hold on;
  plot(domain(2).X([1:end 1],1),domain(2).X([1:end 1],2),'k-','linewidth',2);

  h2 = get(h1,'CurrentAxes');
  h3 = axes('pos',[.36 .36 .31 .31]);

  XX = reshape(Xv(:,shot),2*n,nv);
  [ra a l] = reducedVolume(XX);
  
  errorA = [errorA max(abs(a-a0)./a0)];
  errorL = [errorL max(abs(l-l0)./l0)];
  
  set(h1,'CurrentAxes',h2)
  plot(domain(1).X([1:end 1],1),domain(1).X([1:end 1],2),'k-','linewidth',2);
  hold on;
  plot(domain(2).X([1:end 1],1),domain(2).X([1:end 1],2),'k-','linewidth',2);
  plot(Xb(:,1),Xb(:,2),'r--','linewidth',3);
  for ves =1:nv
    plot(XX([1:n 1],ves),XX([n+1:2*n n+1],ves),'k-','linewidth',1.5,'color',lineColor(ves,:))
  end
  hold off
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  axis([-33 33 -33 33]);
  axis equal; axis off; 
  
  set(h1,'CurrentAxes',h3)
  plot(domain(1).X([1:end 1],1),domain(1).X([1:end 1],2),'k-','linewidth',2);
  hold on;
  plot(domain(2).X([1:end 1],1),domain(2).X([1:end 1],2),'k-','linewidth',2);
  plot(Xb(:,1),Xb(:,2),'r--','linewidth',3);
  for ves =1:nv
    plot(XX([1:n 1],ves),XX([n+1:2*n n+1],ves),'k-','linewidth',1.5,'color',lineColor(ves,:))
  end
  axis([13.84 22.89 13.88 22.94]);
  axis off;
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  hold off
  
  drawnow
  %pause(.2)
  %saveas(gcf,['../results/concentrate' num2str(shot) '.jpg']);
  %pause(.2)
end

disp(errorA)
