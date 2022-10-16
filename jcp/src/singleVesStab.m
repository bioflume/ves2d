clear all; clc

prams.T = .08;
prams.kappa = 1;  
prams.Case = 'schurComp';
prams.vInf = @(X) farFieldVel(X,'shear',250);

options.usePlot = 0;                    
options.progressBar = 0;
options.verbose = 0;
options.GmresTol = 10^-12;    

% Reference solve setup
aRef = 9; NRef = 2^aRef;
X0 = boundary(NRef,'nv',4);
[RV areaRef lengthRef] = reducedVolume(X0);

% Tests
fileName = '../multVesStabTemp.txt';
fileId = fopen(fileName,'a');
midLine = repmat('-',1,85); midLine = [midLine '\n'];

fprintf(fileId,[' lambda\t n,m\t FLAG\t |A-Ar|/Ar\t ' '|L-Lr|/Lr \n']);
fprintf(fileId,midLine);

for o = 2:3
  prams.order = o;
  for vc= [.04 1 25]
    aEr0 = 1; lEr0 = 1;
    for a = 5:8
      N = 2^a;
      prams.n = N; 
      prams.m = 100*N;
      X = boundary(prams.n,'nv',4);
       
      clear functions global;  
      [Xfinal status] = Ves2D(X,prams,options); 
      disp(['test solve ' num2str(status.flag) ' ' num2str(vc) ' ' num2str(N)]);
      
      [RV Area Length] = reducedVolume(Xfinal);
      %Xfinal = reshape(Xfinal,[],2);
      %subInd = 1:NRef/N:NRef;
      %Xerr = Xfinal-XRef(subInd,:); Xerr = max(sqrt(dot(Xerr,Xerr,2)))/XM;
      
      printString = ['%-1.4f\t %-3.0f\t %-1.0f\t %-5.4e(%-2.4f)\t %-5.4e(%-2.4f)\n'];
      aEr = max(abs(Area - areaRef)./areaRef);
      lEr = max(abs(Length - lengthRef)./lengthRef); 
      
      fprintf(fileId,printString,vc,N,status.flag,aEr,aEr0/aEr,lEr,lEr0/lEr); 
      aEr0 = aEr; lEr0 = lEr;
    end
  end
end

fclose(fileId);