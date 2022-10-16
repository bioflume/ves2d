
clear all; clc

nv = 2; 
n = 64;                                                   
% X = boundary(n,'nv',nv,'center',[0 -10;0 0.5],'reducedArea',.98,'scale',.6);

prams.nv = nv;
prams.T = 20;                           
prams.ts = 1e-4;
prams.m = prams.T/prams.ts;
prams.kappa = .1;                        
prams.order = 1;                        
prams.ts = prams.T/prams.m;             
prams.Incompressibility = 1;            
prams.vInf = @(X) farFieldVel(X,'shear',2);

options.usePlot = 0;                    
options.AxesHandle = [];               
options.axis = [];
options.axisOn = 1;
options.progressBar = 0;
options.saveData = 1;
options.datastride = floor(prams.m/80);               

options.verbose = 0;
[trash prams options] = monitor(prams,options);
% 
% % hold on; viewer(X); axis on;
% % axis([-14 4 -4 4])
% % grid on
% %nu = [.1 .2 .5 1 2 5 10 20]
% nu = 50
% for ii=1:length(nu)
%   clear functions global
%   prams.viscCont = nu(ii)*[1 1];
%   options.fileName = ['../results/diffusivityAgain2_' num2str(nu(ii)) '.mat'];
%   [Xfinal status] = Ves2D(X,prams,options);
% end


%% Reading Data
% figure; 
% for nu =[.1 .2 .5 1 2 5 10]
% %  hold on
%   disp(nu);
%   C1 = [];C2 = [];
%   options.fileName = ['../results/diffusivityAgain2S4_' num2str(nu) '.mat'];
% 
%   fileId = fopen(options.fileName,'r');
%   Result = fread(fileId,'double');
%   fclose(fileId);
%   
%   Result = reshape(Result,5*nv*n+1,[]); 
%   X = Result(1:2*n*nv,:);
%   X1 = X(1:2*n,:);  X2 = X(2*n+1:4*n,:);
%   Time = Result(end,:);
%   
%   [ra0 a0 l0] = reducedVolume(reshape(X(:,1),[],2));
%   for ii=1:length(Time)
%     X11 = reshape(X1(:,ii),[],2);C1(:,ii) = mean(X11)';
%     X22 = reshape(X2(:,ii),[],2);C2(:,ii) = mean(X22)';
%     
%      plot(X11(:,1),X11(:,2),X22(:,1),X22(:,2),C1(1,:),C1(2,:),C2(1,:),C2(2,:));
%      axis([-12 12 -6 6]); axis equal
%      grid on; grid minor
% %  
% %      [ra a l] = reducedVolume(reshape(X(:,ii),[],2));
% %      disp('-----------------------------')
% %      disp(abs(a./a0-1));
% %      disp(abs(l./l0-1));
% %      pause(.1);
%   end
% 
%   if(nu ==.1)C2 = C2(:,2:end); end
%   xx = -10:.25:10;
%   yy = spline(C2(1,:)',C2(2,:)',xx);
% 
%  
% %  pause
% %  plot(C2(1,:),C2(2,:),'o',xx,yy);
% %  pause
% end
% hold off


%% Plotting
% ww = .17; W = .18;
% hh = .25; H = .26;
% nu = [.1 1 10];
% JJ = [23 30 34 37 32];
% [xg yg] = meshgrid(-6:.1:4,-3:.1:3); Xg = [xg(:) yg(:)];
% cline = -8:.2:8;
% 
% figure('Units','inches','position',[3 3 6 3]);
% for ii = 1:size(nu,2)
%   disp(nu(ii));
%   options.fileName = ['../results/diffusivityAgain2S4_' num2str(nu(ii)) '.mat'];
% 
%   fileId = fopen(options.fileName,'r');
%   Result = fread(fileId,'double');
%   fclose(fileId);
%   
%   Result = reshape(Result,5*nv*n+1,[]); 
%   NN = size(Result,2);
%   
%   for jj = 1:size(JJ,2)
%     X = reshape(Result(1:2*n*nv,JJ(jj)),[],nv);
%     sig = reshape(Result(2*n*nv+1:3*n*nv,JJ(jj)),[],nv);
%     u = reshape(Result(3*n*nv+1:5*n*nv,JJ(jj)),[],nv);
%     Time = Result(end,JJ(jj));
%     
%     position = [.05+W*(jj-1) .1+H*(ii-1) ww hh];
%     h = subplot('position',position); hold on; 
%     if(jj==5)
%       [F FX] = InteractionForce(X,u,sig,prams,options,prams.viscCont,Xg');
%       FX = reshape(prams.vInf(Xg(:)),[],2)+FX';
%       
%       uu = reshape(FX(:,1),size(xg));   
%       vv = reshape(FX(:,2),size(xg));  
%       h = streamslice(xg,yg,uu,vv,2,'cubic');
%       set(h,'Color',[.6 .6 .6]);
%     else
%       plot(cline,0,'.k','markersize',3); 
%     end
%     X1 = interpft(X(1:n,:),256);
%     X2 = interpft(X(n+1:2*n,:),256);
%     X = [X1;X2];
%     mm = 256;
%     plot(X([1:mm 1],1),X([mm+1:2*mm mm+1],1),'k','linewidth',3.2);
%     plot(X([1:mm 1],2),X([mm+1:2*mm mm+1],2),'k','linewidth',3.2);
%     hold off
%     
%     if(jj==1)
%       axis([-7 2 -1 1]); 
%     else
%       axis([-5 4 -1 1]); 
%     end
%     axis equal;
%     box on;
%     
%     if(ii==3),  title(['t = ' num2str(Time)]);end
%     if(jj==1), ylabel(['\nu = ' num2str(nu(ii))]);end
%     
%     if(ii~=1)
%       set(gca,'xtick',[]);
%     else
%       set(gca,'xtick',[-6 -3 0 3 6]);
%     end
%     if(jj~=1), set(gca,'ytick',[]);end
%   end 
% end

%% Centerline
lineSpec = {'ok','*k','sk','hk','dk','xk','pk'};
lineSpecGhost = {'-ok','-*k','-sk','-hk','-dk','-xk','-pk'};
nu = [.1 .2 .5 1 2 5 10];
  
h1 = figure;
set(gcf,'Units','centimeters');
%set(gcf,'position',[5 5 7.5 5]);
hold on;
for ii = 1:size(nu,2)
  options.fileName = ['../results/diffusivityAgain2S4_' num2str(nu(ii)) '.mat'];

  fileId = fopen(options.fileName,'r');
  Result = fread(fileId,'double');
  fclose(fileId);
  
  Result = reshape(Result,5*nv*n+1,[]); 
  X = Result(1:2*n*nv,:);
  
  C1x(:,ii) = mean(X(1:n,:));
  C1y(:,ii) = mean(X(n+1:2*n,:));
  C2x(:,ii) = mean(X(2*n+1:3*n,:));
  C2y(:,ii) = mean(X(3*n+1:4*n,:));

  xx = -10:.1:10;
  yy(:,ii) = spline(C2x(:,ii),C2y(:,ii)-C1y(:,ii),xx);
  
  plot(xx,yy(:,ii)+10,lineSpecGhost{ii});
  nuCell{ii} = ['\nu = ' num2str(nu(ii))];
end
axis([-10 10 0.4 2.6]); box on
legend(nuCell{:},'location','NW');%[.2 .6 .1 .1]); 

for ii = 1:size(nu,2)
  plot(xx,yy(:,ii),'k',xx(1:10:end),yy(1:10:end,ii)',lineSpec{ii});
end


%set(gca,'xtick',[]);
%set(gca,'ytick',[]);
xlabel('$x$'); ylabel('$\delta$');

h2 = get(h1,'CurrentAxes');
h3 = axes('pos',[.62 .6 .26 .26]);

set(h1,'CurrentAxes',h3)

nuInt = 10.^(-1:.1:1);
yInt = spline(nu,yy(end,:),nuInt);

semilogx(nu,yy(end,:),'ok',nuInt,yInt,'k');
axis([.1 10 0.5 1]); box on

xlabel('\nu'); ylabel('\delta_\infty');
set(gca,'ytick',[.6 .8 1]);

hold off
