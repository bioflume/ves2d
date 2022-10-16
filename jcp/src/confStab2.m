clear all;clc

%% Setting up the prams and options for the confined flow test
n = 128; nv = 1;
prams.T = 10;                           
prams.ts = .005;
prams.m = prams.T/prams.ts;
prams.order = 1;
prams.kappa = 5e-1;                     
prams.Incompressibility = 1;            
prams.method = [];
prams.M = 400;
prams.bd = @(ind,m) sampleBd(ind,m,1,'choke');      
prams.bc = @(x) forcing(x,'choke');
prams.flowType = 'confined';           
prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams);

bound = @(n) fixedBound(n,prams.bd,1);
domain = bound(prams.M);

options.usePlot = 1;                    
options.AxesHandle = [];
options.saveFig = 0;    
options.axisOn = 0;                     
options.axis = [];                      
options.progressBar = 1;
options.showError = 1;
options.verbose = 1;
options.saveData = 1;
options.dataStride = prams.m/20;               

prams.viscCont = 5;
%% Running cases of simulations and saving data to file
for scale = [.6:.1:1.1 1.15]
  clear functions global;
  options.fileName = ['../results/confStab5_' num2str(scale) '.mat'];
  gg = (1:n)'*2*pi/n; X = [-6+scale*1.5*cos(gg);scale*sin(gg)];
  [Xfinal status] = Ves2D(X,prams,options,@monitor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% changing prams and options for the unbounded flow benchmark
prams.flowType = 'unbounded';           
prams.vInf = @(X) farFieldVel(X,'poiseuille',1);  

% [x y] = meshgrid(-10:10,-4:4);
% u = prams.vInf([x(:);y(:)]); u = reshape(u,[],2);
% quiver(x(:),y(:),u(:,1),u(:,2));

options.axisOn = 1;                     
options.saveData = 1;
options.dataStride = prams.m/20;               
options.fileName = ['../results/confStabUnbound5.mat'];

scale = .9; 
gg = (1:n)'*2*pi/n; X = [-6+scale*1.5*cos(gg);1+scale*sin(gg)];
clear functions global;
[Xfinal status] = Ves2D(X,prams,options,@monitor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading data files and extracting data
Xd = domain.X; C = Xd(101,2); Xd = [Xd;Xd(1,:)];
Xb = Xd(abs(Xd(:,1))>9,:);
Xb = Xb(1:8:end,:);
U = prams.bc(Xb);
Scale = [.6:.1:1.1 1.15];
lineSpec = {'-k','--k',':k','-.k'};
A = []; L = []; Sprint = [];

figure('Units','inches','position',[3 3 6 3.8]);
ww = .45; hh = .21;
W = .45;H = .22;
COL = [0 1 0 1 0 1 0];
ROW = [3 3 2 2 1 1 0];

for jj = 1:length(Scale)
  scale = Scale(jj);
  clear functions global;
  fileName = ['../results/confStab5_' num2str(scale) '.mat'];
  
  fileId = fopen(fileName,'r');
  Result = fread(fileId,'double');
  fclose(fileId);
  
  Result = reshape(Result,5*nv*n+1,[]); 
  
  Time = Result(end,:);
  X = Result(1:2*nv*n,:);

  col = COL(jj);row = ROW(jj);
  
  position = [.08+W*col .1+H*row ww hh];
  h = subplot('position',position); 
  
  hold on;
  plot(Xd(:,1),Xd(:,2),'-k','linewidth',2);
  quiver(Xb(:,1),Xb(:,2),U(:,1),U(:,2),0,'k');
  
  NN = size(X,2);
  
  if(NN>10)
    KK = [1:3:7 8:1:11 11:3:NN];
  else
    KK = 1:2:NN;
  end
  
  for kk=KK 
    plot(X([1:n 1],kk),X([n+1:2*n n+1],kk),'-k');
  end
  hold off;
  
  axis([-11 12 -3.2 3.2]);axis off;
  DD = 1.25*scale/C; % DD =sum of ellipses semiaxis/choke size;
    
  text(-6,1.5,['r = ' num2str(DD,3)],'VerticalAlignment', 'bottom', ...
       'HorizontalAlignment', 'center');
  
  [ra aa ll] = reducedVolume(X);
  print = [Sprint DD];
  A = [A abs(aa(end)-aa(1))/aa(1)];
  L = [L abs(ll(end)-ll(1))/ll(1)];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading the benchmark
symLine = linspace(-10,10,160)'; 
symLine = [symLine 0*symLine];

Xb = Xb(Xb(:,1)<0,:); Xb(:,1) = -10;
%U = prams.vInf(Xb(:)); U = reshape(U,[],2);

fileName = ['../results/confStabUnbound5.mat'];
fileId = fopen(fileName,'r');
Result = fread(fileId,'double');
fclose(fileId);

Result = reshape(Result,5*nv*n+1,[]); 

Time = Result(end,:);
X = Result(1:2*nv*n,:);

col = 1;row = 0;
position = [.08+W*col .1+H*row ww hh];
h = subplot('position',position); 

hold on;
%quiver(Xb(:,1),Xb(:,2),U(:,1),U(:,2),0,'k');

plot(symLine(:,1),symLine(:,2),'.k','markersize',2);


NN = size(X,2);
KK = [1:3:NN-1 NN];

for kk=KK 
  plot(X([1:n 1],kk),X([n+1:2*n n+1],kk),'-k');
end
hold off;

axis([-11 12 -3.2 3.2]);axis off;

DD = 1.25*.9/C;
text(-6,2,['r = ' num2str(DD,3)],'VerticalAlignment', 'bottom', ...
     'HorizontalAlignment', 'center');

[ra aa ll] = reducedVolume(X);
Sprint = [DD Sprint];
A = [abs(aa(end)-aa(1))/aa(1) A];
L = [abs(ll(end)-ll(1))/ll(1) L];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Printing the errors
PS  = ['  %-2.3e & '];
PS = [repmat(PS,1,length(A)) '\n'];

disp('----------------------------------------------------------')
fprintf(['r     ' PS],Sprint)
disp('----------------------------------------------------------')
fprintf(['DA/A0 ' PS],A)
disp('----------------------------------------------------------')
fprintf(['DL/L0 ' PS],L)
disp('----------------------------------------------------------')

% vc = 1
% ----------------------------------------------------------
% r       5.833e-01 &   6.806e-01 &   7.778e-01 &   8.750e-01 &   9.722e-01 &   1.069e+00 &   1.118e+00 & 
% ----------------------------------------------------------
% DA/A0   7.127e-03 &   2.187e-04 &   3.694e-04 &   4.166e-04 &   4.267e-04 &   4.296e-04 &   2.983e-03 & 
% ----------------------------------------------------------
% DL/L0   9.895e-05 &   5.254e-04 &   8.208e-04 &   9.329e-04 &   1.048e-03 &   1.188e-03 &   1.231e-03 & 
% ----------------------------------------------------------

% vc = .2
% ----------------------------------------------------------
% r       8.750e-01 &   5.833e-01 &   6.806e-01 &   7.778e-01 &   8.750e-01 &   9.722e-01 &   1.069e+00 &   1.118e+00 & 
% ----------------------------------------------------------
% DA/A0   4.736e-04 &   1.293e-04 &   1.300e-04 &   1.175e-04 &   8.026e-05 &   3.188e-05 &   3.573e-04 &   1.903e-04 & 
% ----------------------------------------------------------
% DL/L0   1.951e-04 &   2.583e-04 &   2.478e-04 &   2.413e-04 &   2.396e-04 &   2.424e-04 &   1.933e-04 &   1.972e-04 & 
% ----------------------------------------------------------

% vc = 5
% ----------------------------------------------------------
% r       8.750e-01 &   5.833e-01 &   6.806e-01 &   7.778e-01 &   8.750e-01 &   9.722e-01 &   1.069e+00 &   1.118e+00 & 
% ----------------------------------------------------------
% DA/A0   2.784e-04 &   2.859e-05 &   2.022e-06 &   6.450e-06 &   1.460e-05 &   4.351e-05 &   1.893e-04 &   2.669e-04 & 
% ----------------------------------------------------------
% DL/L0   1.033e-04 &   7.174e-04 &   2.980e-04 &   1.683e-04 &   1.058e-04 &   1.436e-04 &   1.690e-04 &   1.950e-04 & 
% ----------------------------------------------------------
