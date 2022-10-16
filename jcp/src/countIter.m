clear functions global;

prams.kappa = 1e0;  
prams.order = 1;    
prams.Incompressibility = 1;
prams.vInf = @(X) farFieldVel(X,'shear',1);
prams.Case = 'schurComp';

% Setting Options
options.GmresTol = 10^-8;  
options.usePlot = 1;                   

clear functions global;
n = 64; 

X = boundary(n);

prams.n = n; 
prams.T = 6;
prams.m = 60;
prams.viscCont = 1e0;

[Xfinal status] = Ves2D(X,prams,options);

