clear all; clc
%dbstop error;

% Get the initial marker point location
nv = 5;                                 % Number of vesicles
n = 64;                                 % Number of discretization points
angle = 0;%-pi/8;                          % Initial inclination angle of vesicles
X = boundary(n,'nv',nv,'angle',angle);%,'center',[-2 0 2;0 0 0]);

% Setting modeling parameters
prams.nv = nv;
prams.T = 8;                            % Simulation time
prams.m = 200;                           % Number of time steps
prams.kappa = 1;                        % Bending modulus
prams.viscCont = ones(1,nv);
prams.order = 1;                        % order of the slover
prams.ts = prams.T/prams.m;             % step size (time)
prams.Incompressibility = 1;            % Incompressibility
prams.vInf = @(X) farFieldVel(X,'shear',2);
                                        % The far field ve locity
prams.Case = 'schurComp';
% Setting Options
options.usePlot = 1;                    % To plot the intermidate states   
options.AxesHandle = gca;               % The axes for the  plot
options.axis = [-8 8 -5 5];             % Force the range of the axes
                                        % inset of the figure                                       
options.progressBar = 1;
options.useGPU = 0;

options.fileName = ['./trash.mat'];
options.saveData = 1;
options.dataStride = prams.m/20;  
[t prams options] = monitor(prams,options);
% Calling update function with all of the above parameters 
%[Xfinal status] = Ves2D(X,prams,options);

fileId = fopen(options.fileName,'r');
Result = fread(fileId,'double');
fclose(fileId);

Result = reshape(Result,5*n*nv+1,[]);
[xg yg] = meshgrid(-8:.05:8,-3:.05:3); Xg = [xg(:) yg(:)];

for ii = 11
  X = reshape(Result(1:2*n*nv,ii),[],nv);
  sig = reshape(Result(2*n*nv+1:3*n*nv,ii),[],nv);
  u = reshape(Result(3*n*nv+1:5*n*nv,ii),[],nv);

  [F FX] = InteractionForce(X,u,sig,prams,options,prams.viscCont,Xg');
  FX = reshape(prams.vInf(Xg(:)),[],2)+FX';

  uu = reshape(FX(:,1),size(xg));   
  vv = reshape(FX(:,2),size(xg));   

  viewer(X);
  streamslice(xg,yg,uu,vv,10,'cubic');
  pause
end