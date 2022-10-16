
% case 1: one vesicle, shear flow, no contrast
% case 2: one vesicle, extensional flow, no contrast
% case 3: one vesicle, shear flow, contrast = 10 (enough to make it tumble)
% case 4: two vesicle, shear flow, no contrast

CASE = 4;
 
switch CASE
 case 1
  clear functions global;

  testArea = 50*pi;
  n = 64;                                 
  X = boundary(n);
  
  prams.T = 6;                            
  prams.m = 60;                           
  prams.viscCont = 1;                     
  prams.order = 3;                        
  prams.vInf = @(X) farFieldVel(X,'shear',1);

  options.usePlot = 0;                    
  options.axis = [];
  options.track = 1;                      
  options.progressBar = 0;

  [Xfinal status avgS] = Ves2D(X,prams,options);

  avgSinf = [0 1;1 0]/2;
  avgSinf = repmat(avgSinf(:),1,prams.m);
  avgS = avgSinf+avgS(:,2:end)/testArea;

 case 2
  clear functions global;

  testArea = 50*pi;
  n = 64;                                 
  X = boundary(n);
  
  prams.T = 3;                            
  prams.m = 120;                           
  prams.viscCont = 1;                     
  prams.order = 4;                        
  prams.vInf = @(X) farFieldVel(X,'extensional',1);

  options.usePlot = 0;    
  options.progressBar = 0;
  
  [Xfinal status avgS] = Ves2D(X,prams,options);

  avgSinf = [1 0;0 -1]/2;
  avgSinf = repmat(avgSinf(:),1,prams.m);
  avgS = avgSinf+avgS(:,2:end)/testArea;
 
 case 3
  clear functions global;

  testArea = 50*pi;
  n = 64;                                 
  X = boundary(n);

  prams.T = 12;                            
  prams.m = 240;                           
  prams.viscCont = 10;                     
  prams.order = 1;                        
  prams.vInf = @(X) farFieldVel(X,'shear',1);

  options.usePlot = 0;    
  options.progressBar = 0;
  
  [Xfinal status avgS] = Ves2D(X,prams,options);

  avgSinf = [0 1;1 0]/2;
  avgSinf = repmat(avgSinf(:),1,prams.m);
  avgS = avgSinf+avgS(:,2:end)/testArea;
    
 case 4
  clear functions global;

  testArea = 50*pi;
  n = 64;   
  nv = 2;
  X = boundary(n,'nv',nv);

  prams.T = 6;                            
  prams.m = 120;                           
  prams.viscCont = 1;                     
  prams.order = 1;                        
  prams.vInf = @(X) farFieldVel(X,'shear',1);

  options.usePlot = 0;
  options.progressBar = 0;

  [Xfinal status avgS] = Ves2D(X,prams,options);

  avgSinf = [0 1;1 0]/2;
  avgSinf = repmat(avgSinf(:),1,prams.m);
  avgS = avgSinf+avgS(:,2:end)/testArea;
end

plot(avgS');
