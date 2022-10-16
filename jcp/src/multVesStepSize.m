clear all; clc

%dbstop in multVesStepSize at 36
prams.kappa = 1e0;  
prams.order = 2;
prams.Case = 'schurComp';

options.usePlot = 0;                    
options.progressBar = 0;
options.verbose = 0;

fileName = '/home/abtin/works/Ves2D/jcp/results/multVesStepSizeO2.txt';
fileId = fopen(fileName,'a');
if(fileId== -1);return;end
midLine = repmat('-',1,75); midLine = [midLine '\n'];

fprintf(fileId,[' shear\t\t lambda\t\t n\t\t m\t\t ts\n']);
fprintf(fileId,midLine);

for s = [1 10]
  prams.T = 6/s;
  prams.vInf = @(X) farFieldVel(X,'shear',s);
  for vc = [.01 .1 1 10 100]
    prams.viscCont = vc;
    for a = 5:8
      prams.n = 2^a;
      X = boundary(prams.n,'nv',2,'angle',-pi/8);
      
      iterate = true;
      m = 16;upperM=[];lowerM=[];
      while(iterate)
        clear functions global;
        m = ceil(m); prams.m = m;
        [Xfinal status] = Ves2D(X,prams,options); 
        
        if(status.flag)
          if(isempty(upperM))
            upperM = m;
            lowerM = m/2;
          else
            upperM = m;
          end
          m = (upperM+lowerM)/2;
        else
          if(isempty(upperM))
            m = 1.4*m;
          else
            lowerM = m;
            m= (upperM+lowerM)/2;
          end
        end
        
        if(~isempty(upperM) && (upperM-lowerM)<max(ceil(.1*upperM),2))
          iterate = false;
        elseif(isempty(upperM) && m>1e7)
          iterate = false;
        end
        disp([lowerM m upperM]);
      end
      printString = [' %-3.0f\t\t %0.2e\t %-3.0f\t\t %-4.0f\t\t %-5.4e\n'];
      fprintf(fileId,printString,s,vc,prams.n, prams.m,prams.T/prams.m); 
    end
  end
end

fclose(fileId);
