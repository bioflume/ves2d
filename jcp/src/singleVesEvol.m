addpath ../../src
clear all;clc

prams.T = 1;
prams.kappa = 1e0;
prams.order = 1;   
prams.rhoIn = 2;
prams.rhoOut = 1;
prams.Incompressibility = 1;
prams.Case = 'schurComp';
prams.g = [0 -40];

options.GmresTol = 10^-8;  
options.usePlot = 0;
options.verbose = 0;                   
options.AxesHandle = [];               
options.axisOn = 0;                    
options.track = 0;          
options.axis = [];
options.saveFig = 0;                   
options.velProfile = 0;               
options.progressBar = 0;
options.saveData = 1;
options.showError = 0;               

[trash prams options] = monitor(prams,options);

% for ra = .75
%   for shear= 0%[0 4]
%     prams.vInf = @(X) farFieldVel(X,'shear',shear);
%     for ii=1
%       nu=[.01 1 100];
%       n = 128; prams.m = 2000;
%       X0 = boundary(n,'angle',-shear*pi/16);
%       prams.viscCont = nu(ii);
% 
%       iterate = true;
%       while(iterate)
%         options.fileName = ['../results/singVesEvol_' num2str(100*ra) '_' ...
%                             num2str(shear) '_' num2str(ii) ];
%         if(exist(options.fileName,'file')), fid = fopen(options.fileName,'w');fclose(fid);end
% 
%         prams.m = 2*prams.m; 
%         options.dataStride = floor(prams.m/5);               
%         
%         clear functions global;
%         [Xfinal status] = Ves2D(X0,prams,options);
%         if(~isempty(lastwarn))
%           [wmsg wid] = lastwarn;
%           lastwarn('');
%           warning('off',wid);
%         end
%         if(status.flag), iterate = false;end
%         fprintf('%3.3f\t %2.2f\t %2.2f\t %3.4e\t %1.1f\n',ra,shear,ii,prams.m,status.flag);
%       end
%     end
%   end
% end

% % Reading the saved data from the file.
status.fileName = ['../results/singVesEvol_RA2_2']; n =128;nu = [.05 1 50];
fileId = fopen(status.fileName,'r');
Result = fread(fileId,'double');
fclose(fileId);

Result = reshape(Result,5*n+1,[]); 
 
map0 = colormap(jet); map = map0(1:128/n:64,:); map = [map;flipud(map)];
Time = Result(end,:);
 
Result = Result(1:2*n,:);
NN = size(Result,2);

ww = .125;
W = .13;
figure('Units','inches','position',[3 3 5 4]);
for ii=0:NN-1
  col = mod(ii,NN/3);
  row = floor(3*ii/NN)+1;
 
  position = [.05+W*col W*row ww ww];
  h = subplot('position',position); 
  XX = Result(:,ii+1);
  plot(XX([1:n 1]),XX([n+1:end n+1]),'k','LineWidth',2.12);
%   plot(XX(1:n),XX(n+1:end),'k','LineWidth',.5);
%   hold on;
%   for k=1:2:n
%     plot(XX(k),XX(n+k),'-','MarkerSize',10,'Color',map(k,:));
%   end
%   hold off;
  if(row==3),  title(['t = ' num2str(Time(ii+1))]);end
  axis([-3.6 3.6 -3.2 3.2]);%axis equal;%axis off;
  if(mod(ii,NN/3)==0), ylabel(['\nu = ' num2str(nu(row))]);end
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  %if(ii<9), set(gca,'xtick',[]);endend
end

%% FREE FALL
% ww = .145; W = .15; H = .14; hh = .135;
% nu=[.01 1 100]; n = 128;     
% map0 = colormap(jet); map = map0(1:128/n:64,:); map = [map;flipud(map)];
% %figure('Units','inches','position',[3 3 6 4]);
% for ra = .75
%   for shear= 0
%     for row=1:3
%       options.fileName = ['../results/singVesEvol_' num2str(100*ra) '_' ...
%                           num2str(shear) '_' num2str(row) ];
%       fileId = fopen(options.fileName,'r');
%       Result = fread(fileId,'double');
%       fclose(fileId);
%       
%       Result = reshape(Result,5*n+1,[]); 
%       NN = size(Result,2);
%           
%       for col=0:NN-1
%         position = [.05+W*col H*row ww hh];
%         subplot('position',position); hold on
%         
%         X = Result(1:2*n,col+1); Cx =mean(X(1:n)); Cy=mean(X(n+1:2*n));
%         sig = Result(2*n+1:3*n,col+1);
%         u = Result(3*n+1:5*n,col+1);
%         Time = Result(end,col+1);
%         
%         if(col==5)
%           [xg yg] = meshgrid(Cx-3:.1:Cx+3,Cy-3.2:.1:Cy+2.2); Xg = [xg(:) yg(:)];
%           [F FX] = InteractionForce(X,u,sig,prams,options,prams.viscCont,Xg');
%           FX = FX';
%           uu = reshape(FX(:,1),size(xg));   
%           vv = reshape(FX(:,2),size(xg));  
%           h = streamslice(xg,yg,uu,vv,1,'cubic');
%           set(h,'Color',[.6 .6 .6]);
%         else
%           if(abs(Cy)<1e-3), Cy = 0;end
%           text(Cx,Cy-2.9,['c_y =' num2str(Cy,2)],'HorizontalAlignment', ...
%                'center');
%         end
%         X1 = interpft(X(1:n),256);
%         X2 = interpft(X(n+1:2*n),256);
%         plot(X1,X2,'k','linewidth',2.1);
%         hold off
% 
%         if(row==3),  title(['t = ' num2str(Time)]);end
%         axis([Cx-3 Cx+3 Cy-3.2 Cy+2.2]);%axis equal;%axis off;
%         if(col==0), ylabel(['\nu = ' num2str(nu(row))]);end
%         set(gca,'xtick',[]);
%         
%         set(gca,'ytick',[]);
%         box on;
%         %if(col<9), set(gca,'xtick',[]);endend
%       end
%     end      
%   end
% end
