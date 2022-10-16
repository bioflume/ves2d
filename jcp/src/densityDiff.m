clear functions global;

n = 64; 
X = boundary(n,'angle',0,'reducedArea',.2,'scale',.4);

% Setting modeling parameters
prams.T = 1;
prams.m = 800;
prams.viscCont = .01;
prams.kappa = 1e0;  
prams.order = 1;
prams.g = [0 -100];
prams.rhoIn = 2;
prams.rhoOut = 1;
prams.Incompressibility = 1;
prams.vInf = @(X) farFieldVel(X,'shear',10);
                            
% Setting Options
options.GmresTol = 10^-8;  
options.usePlot = 1;                   
options.AxesHandle = [];               
options.axisOn = 1;                    
options.track = 0;     
options.showAngle = 0;

options.axis = [-2 2 -2 2];
options.saveFig = 0;                   
                                       
options.progressBar = 0;
options.saveData = 0;                  
options.showError = 1;               

% Calling update function with all of the above parameters 
prams.n = n; prams.nv = 1;
prams.Case = 'schurComp';
[Xfinal status] = Ves2D(X,prams,options);
