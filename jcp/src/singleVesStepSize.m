addpath ../../src
clear all; clc

prams.kappa = 1;  
prams.order = 1; 
prams.Case = 'schurComp';
prams.method = 'explicit';

options.GmresTol = 10^-6;  
options.usePlot = 0;                    
options.progressBar = 0;
options.verbose = 0;
  
fileName = '../results/singleVesExp.txt';
fileId = fopen(fileName,'a');
if(fileId== -1);return;end
midLine = repmat('-',1,75); midLine = [midLine '\n'];

fprintf(fileId,[' shear\t\t lambda\t\t n\t\t m\t\t ts\n']);
fprintf(fileId,midLine);

for s = [10]
  prams.T = .5/s;
  prams.vInf = @(X) farFieldVel(X,'shear',s);
  for vc = [1]
    prams.viscCont = vc;
    for a = 5:8
      prams.n = 2^a;
      X = boundary(prams.n);
      
      iterate = true;
      m = 32;upperM=[];lowerM=[];
      while(iterate)
        clear functions global;
        m = ceil(m); prams.m = m;
        [Xfinal status] = Ves2D(X,prams,options); 
        
        if(~isempty(lastwarn))
          [wmsg wid] = lastwarn;
          lastwarn('');
          warning('off',wid);
        end
        
        if(~status.flag)
          m = 1.7*m;
        else
          iterate = false;
          prams.m = ceil(2.7*m/3.4);
        end
        
        if(m>1e10)
          iterate = false;
        end
        disp(prams.m);

%         if(status.flag)
%           if(isempty(upperM))
%             upperM = m;
%             lowerM = m/2;
%           else
%             upperM = m;
%           end
%           m = (upperM+lowerM)/2;
%         else
%           if(isempty(upperM))
%             m = 1.7*m;
%           else
%             lowerM = m;
%             m= (upperM+lowerM)/2;
%           end
%         end
%   
%         if(m>1e10 || ~isempty(upperM) && (upperM-lowerM)<max(.01*upperM,2))
%           iterate = false;
%         end
%         disp([lowerM prams.m upperM]);
      end
      disp('=============================================');
      printString = [' %-3.0f\t\t %0.2e\t %-3.0f\t\t %-4.0f\t\t %-5.4e\n'];
      fprintf(fileId,printString,s,vc,prams.n, prams.m,prams.T/prams.m); 
   end
 end
end

fclose(fileId);
exit