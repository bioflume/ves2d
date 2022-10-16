% test1 test the convergence in case of shear flow by changing the
% order and the number of discretization points.

clear all; clc
prams.T = 6;
prams.viscCont = 1e0;
prams.kappa = 1e0;  
prams.vInf = @(X) farFieldVel(X,'shear',1);

options.usePlot = 0;                    
options.progressBar = 0;
options.verbose = 0;
  
prams.order = 1; prams.Case = 'schurComp';
aRef = 7; NRef = 2^aRef;
X = boundary(NRef);
[RV areaRef lengthRef] = reducedVolume(X);
[XRef status] = Ves2D(X,prams,options); 
XRef = reshape(XRef,[],2); 
XM = max(sqrt(dot(XRef,XRef,2)));

midLine = repmat('-',1,85); midLine = [midLine '\n'];
fprintf(' order\t method \t n,m\t FLAG\t |X-Xr|/|Xr|\t |A-Ar|/Ar\t |L-Lr|/Lr \n');
fprintf(midLine);

cases = {'schurComp','fullSystem'};
for a= 5:7
  N = 2^a;
  prams.m = N;
  prams.n = N; 
  
  X = boundary(prams.n);
  for O = 1:2
    prams.order = O;
    for c = 1:2
      prams.Case = cases{c};
      clear global functions;  
      [Xfinal status] = Ves2D(X,prams,options); 
      [RV Area Length] = reducedVolume(Xfinal);
      Xfinal = reshape(Xfinal,[],2);
  
      subInd = 1:2^(aRef-a):NRef;
      Xerr = Xfinal-XRef(subInd,:); Xerr = max(sqrt(dot(Xerr,Xerr,2)))/XM;
      printString = [' %-1.0f\t ' cases{c} '\t %-4.0f\t %-1.0f\t ' ...
      '%-5.4e\t %-5.4e\t %-5.4e\n'];
      fprintf(printString, prams.order,N, status.flag,Xerr,abs(Area - ...
              areaRef)/areaRef,abs(Length - lengthRef)/lengthRef); 
    end
  end
end
