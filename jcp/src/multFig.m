clear all; clc

prams.kappa = 1e0;  
prams.order = 2;
prams.Case = 'schurComp';

options.usePlot = 1;                    
options.progressBar = 0;
options.verbose = 0;
options.saveData = 1;
options.fileName = ['../results/multVesInteract.mat'];
options.useGPU = 0;

s = 1;
prams.T = 6/s;
prams.vInf = @(X) farFieldVel(X,'shear',s);

vc = 1;
prams.viscCont = vc;

n = 64;nv = 2; 
prams.n = n;
X = boundary(prams.n,'nv',nv,'angle',-pi/8,'center',[[-1.1;0] [1.1;0]]);

m = 180; prams.m = m;
options.dataStride = floor(prams.m/5);               

%[Xfinal status] = Ves2D(X,prams,options); 

fileId = fopen(options.fileName,'r');
Result = fread(fileId,'double');
fclose(fileId);

Result = reshape(Result,5*nv*n+1,[]);
Result = Result(:,2:end); %discarding the first config(it is not physical
[xg yg] = meshgrid(-3.2:.1:3.2,-2.5:.1:2.5); Xg = [xg(:) yg(:)];

ww = .15; W = .18; hh = .6;
figure('Units','inches','position',[3 3 5 1.3]);
for jj = 1:size(Result,2)
  X = reshape(Result(1:2*n*nv,jj),[],nv);
  sig = reshape(Result(2*n*nv+1:3*n*nv,jj),[],nv);
  u = reshape(Result(3*n*nv+1:5*n*nv,jj),[],nv);
  Time = Result(end,jj);
  
  position = [.05+W*(jj-1) .18 ww hh];
  h = subplot('position',position); hold on; 

  [F FX] = InteractionForce(X,u,sig,prams,options,prams.viscCont,Xg');
  FX = reshape(prams.vInf(Xg(:)),[],2)+FX';
  
  uu = reshape(FX(:,1),size(xg));   
  vv = reshape(FX(:,2),size(xg));  
  h = streamslice(xg,yg,uu,vv,2,'cubic');
  set(h,'Color',[.6 .6 .6]);

  X1 = interpft(X(1:n,:),256);
  X2 = interpft(X(n+1:2*n,:),256);
  X = [X1;X2];
  mm = 256;
  hold on;
  plot(X([1:mm 1],1),X([mm+1:2*mm mm+1],1),'k','linewidth',3.2);
  plot(X([1:mm 1],2),X([mm+1:2*mm mm+1],2),'k','linewidth',3.2);
  hold off;
  axis([-3.1 3.1 -2.4 2.4]); 
  box on;
  title(['t = ' num2str(Time)]);
  %if(jj==1), ylabel(['\nu = ']);end
  set(gca,'xtick',[-3 0 3]);
  if(jj~=1)
    set(gca,'ytick',[]);
  else
    set(gca,'ytick',[-2 0 2]);
  end
end
