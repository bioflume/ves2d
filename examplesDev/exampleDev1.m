clear functions global;

%dbstop in initTemp at 47
n = 64; 
X = boundary(n);

% Setting modeling parameters
prams.T = 1;
prams.m = 50;
%prams.m = 20000;
prams.m = 20000;
prams.viscCont = 1;
prams.kappa = 1e0;  
prams.order = 2;    
prams.Incompressibility = 1;
prams.vInf = @(X) farFieldVel(X,'shear',1);
                            
% Setting Options
options.GmresTol = 10^-12;  
options.usePlot = 1;                   
options.AxesHandle = [];               
options.axisOn = 0;                    
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
%prams.method = 'explicit';
[Xfinal status] = Ves2D(X,prams,options);
