ppclear all; clc

prams.T = 6;                           
prams.order = 1;
prams.kappa = 1e-1;                     
prams.Incompressibility = 1;            
prams.method = [];%'explicit';
                  %prams.Case = 'fullSystem';

options.usePlot = 1;                    
options.AxesHandle = [];
options.saveFig = 0;    
options.axisOn = 0;                     
options.axis = [];                      
options.progressBar = 1;
options.showError = 1;
options.verbose = 0;

prams.viscCont = 1;
n = 64;
prams.m = 400;
prams.ts = prams.T/prams.m;
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette');
prams.bc = @(x) forcing(x,'couette');   
prams.M = [128 196];

c = sqrt(30*16);

for s = 1:.2:3.6
  X = boundary(n,'couette','nv',1,'scale',s,'angle',pi/2,'center',[0;c]);

  clear functions global;
  prams.flowType = 'confined';           
  prams.vInf = @(X,bc) farFieldVel(X,bc,prams.flowType,prams);
  [XCouette status1] = Ves2D(X,prams,options,@monitor);

  clear functions global;
  X(n+1:end) = X(n+1:end)-c;
  prams.flowType = 'unbounded';            
  prams.vInf = @(X) farFieldVel(X,'shear',1);
  [XShear status2]= Ves2D(X,prams,options,@monitor);

  if(status1.flag && status2.flag)
    [r0 a0 l0] = reducedVolume(X);
    [rc ac lc] = reducedVolume(XCouette);
    [rs as ls] = reducedVolume(XShear);
    
    fprintf('%2.1f\t %4.3e\t %4.3e\t %4.3e\t %4.3e\n',s,abs(as-a0)/a0, ...
            abs(ac-a0)/a0,abs(ls-l0)/a0,abs(lc-l0)/l0); 
  else
    disp('diverged!');
  end
end


% --------------------------------------------------------------------------
%  scale    Area - shear   Area - couette  Length - shear  Length - couette
% --------------------------------------------------------------------------
%   1.0      8.367e-03       2.048e-03       1.050e-02       6.197e-03
%   1.3      8.388e-03       1.797e-03       7.882e-03       6.004e-03
%   1.6      8.401e-03       1.486e-03       6.336e-03       5.912e-03
%   1.9      8.408e-03       1.135e-03       5.309e-03       5.887e-03
%   2.2      8.413e-03       7.719e-04       4.574e-03       5.962e-03
%   2.5      8.416e-03       4.071e-03       4.020e-03       6.221e-03
% --------------------------------------------------------------------------
%   2.0      4.919e-03       1.347e-03       4.239e-03       6.369e-03
%   2.4      4.933e-03       1.338e-03       3.530e-03       6.223e-03
%   2.8      4.941e-03       1.264e-03       3.027e-03       6.073e-03
%   3.2      4.945e-03       1.135e-03       2.650e-03       5.941e-03
%   3.6      4.948e-03       9.625e-04       2.357e-03       5.835e-03
%   4.0      4.951e-03       7.524e-04       2.122e-03       5.761e-03
% --------------------------------------------------------------------------


% New simulation -- very similar to the unbounded flow:
% --------------------------------------------------------------------------
% scale    Area - shear   Area - couette  Length - shear  Length - couette
% --------------------------------------------------------------------------
% 1.0      1.584e-02       1.517e-02       1.703e-02       1.401e-02
% 1.2      1.584e-02       1.481e-02       1.412e-02       1.374e-02
% 1.4      1.584e-02       1.438e-02       1.208e-02       1.344e-02
% 1.6      1.584e-02       1.389e-02       1.056e-02       1.310e-02
% 1.8      1.584e-02       1.335e-02       9.377e-03       1.273e-02
% 2.0      1.584e-02       1.280e-02       8.437e-03       1.236e-02
% 2.2      1.584e-02       1.226e-02       7.670e-03       1.202e-02
% --------------------------------------------------------------------------
% 1.0      1.245e-02       1.405e-02       1.620e-02       1.554e-02
% 1.2      1.248e-02       1.382e-02       1.348e-02       1.530e-02
% 1.4      1.250e-02       1.344e-02       1.156e-02       1.503e-02
% 1.6      1.252e-02       1.302e-02       1.012e-02       1.473e-02
% 1.8      1.253e-02       1.266e-02       8.999e-03       1.444e-02
% 2.0      1.253e-02       1.237e-02       8.104e-03       1.420e-02
% 2.2      1.254e-02       1.250e-02       7.371e-03       1.439e-02
% 2.4      1.254e-02       1.520e-02       6.760e-03       1.715e-02
% 2.6      1.255e-02       1.567e-02       6.242e-03       1.813e-02
% 2.8      1.255e-02       8.445e-03       5.798e-03       1.943e-02
% 3.0      1.255e-02       1.971e-03       5.413e-03       1.551e-02
% 3.2      1.255e-02       1.691e-02       5.076e-03       1.863e-02
% 3.4      1.255e-02       2.491e-02       4.779e-03       3.444e-01
% 3.6      1.255e-02       4.538e-03       4.514e-03       1.239e-02
% --------------------------------------------------------------------------
%% Check FLAG
