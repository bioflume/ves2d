clear all; 

n = 64;                                 
nv = 1;
prams.T = 10;                           
prams.m = 100;                          
prams.kappa = 1e-1;
prams.rhoOut = 1;                       
prams.order = 1;                        
prams.ts = prams.T/prams.m;             
prams.flowType = 'confined';            
prams.M = 400;                     
prams.bd = @(ind,m) sampleBd(ind,m,1,'boxed');
prams.bc = @(x) forcing(x,'shear');  
prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams,'direct');
        
% [gg XV] = prams.bd(1,prams.M); XX = .98*XV;
% [x y] = meshgrid(3:.2:5); IN = inpolygon(x(:),y(:),XX(:,1),XX(:,2));
% x = x(IN); y = y(IN);
% 
% prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams,'direct');
% u = prams.vInf([x(:) y(:)],[]); u = reshape(u,[],2);
% e = u(:,1)-y(:); disp(max(abs(e)./abs(y(:)+(y(:)==0))));

options.usePlot = 1; 
options.progressBar = 1;
options.AxesHandle = gca;
options.saveFig = 0;   
options.saveStride = 5;
options.saveData = 0;                  
options.dataStride = 5;
X = boundary(n); 

[XF status] = Ves2D(X,prams,options,@monitor);
