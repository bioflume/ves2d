clear; clc;

% FIXED PARAMETERS IN ALL RUNS
N = 64;
Th = 0.75;
dt = 0.01;
speed = 1;
VC = 1;
kappa = 1e-2;
% vol frac for each set
volFrac = [0.05;0.10;0.15;0.20;0.25];
% number of simulations corresponding to the above volume fractions
numRuns = [200; 200; 200; 200; 200; 200];

iset = 5;

if iset == 1 % all VF = 25% runs
  
  for i = 1 : numRuns(iset)
    runName = ['freeCouetteVF25Run' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(5))
  end
  
elseif iset == 2 % all VF = 25% runs
   
  for i = numRuns(1) + 1 : sum(numRuns(1:2))
    runName = ['freeCouetteVF25Run' num2str(i)];
   freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(5))
  end
   
elseif iset == 3 % 25% runs
  istart = sum(numRuns(1:2))+1;
  iend = sum(numRuns(1:3));
  for i = istart : iend
    runName = ['freeCouetteVF25Run' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(5))
  end
    
elseif iset == 4 % 25% runs
  istart = sum(numRuns(1:3)) +  1;
  iend = sum(numRuns(1:4));
  for i = istart : iend
    runName = ['freeCouetteVF25Run' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(5))
  end
  
elseif iset == 5 % 25% runs
  istart = sum(numRuns(1:4)) +  1;
  iend = sum(numRuns(1:5));
  for i = istart : iend
    runName = ['freeCouetteVF25Run' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(5))
  end    

elseif iset == 6 % 25% runs
  istart = sum(numRuns(1:5)) +  1;
  iend = sum(numRuns(1:6));
  for i = istart : iend
    runName = ['freeCouetteVF25Run' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(5))
  end    
     
end
   














if 0
numRuns = [2000; 2000; 1000; 1000; 750];
if iset == 1 % all VF = 5% runs
  
  for i = 1 : numRuns(iset)
    runName = ['freeCouetteRun' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(1))
  end
  
elseif iset == 2 % all VF = 10% runs
   
  for i = numRuns(1) + 1 : sum(numRuns(1:2))
    runName = ['freeCouetteRun' num2str(i)];
   freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(2))
  end
   
elseif iset == 3 % half of VF = 15% runs
  istart = sum(numRuns(1:2))+1;
  iend = sum(numRuns(1:2)) + numRuns(3)/2;
  for i = istart : iend
    runName = ['freeCouetteRun' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(3))
  end
    
elseif iset == 4 % second half of VF = 15% runs
  istart = sum(numRuns(1:2)) + numRuns(3)/2 + 1;
  iend = sum(numRuns(1:3));
  for i = istart : iend
    runName = ['freeCouetteRun' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(3))
  end
  
elseif iset == 5 % half of VF = 20% runs
  istart = sum(numRuns(1:3)) +  1;
  iend = sum(numRuns(1:3)) + numRuns(4)/2;
  for i = istart : iend
    runName = ['freeCouetteRun' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(4))
  end    

elseif iset == 6 % second half of VF = 20% runs
  istart = sum(numRuns(1:3)) + numRuns(4)/2 + 1;
  iend = sum(numRuns(1:4));
  for i = istart : iend
    runName = ['freeCouetteRun' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(4))
  end    
     
elseif iset == 7 % half of VF = 25% runs
  istart = sum(numRuns(1:4)) +  1;
  iend = sum(numRuns(1:4)) + numRuns(5)/2;
  for i = istart : iend
    runName = ['freeCouetteRun' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(5))
  end    

elseif iset == 8 % second half of VF = 20% runs
  istart = sum(numRuns(1:4)) + numRuns(5)/2 + 1;
  iend = sum(numRuns(1:5));
  for i = istart : iend
    runName = ['freeCouetteRun' num2str(i)];
    freeCouetteRuns(runName,N,Th,dt,speed,VC,kappa,volFrac(5))
  end   
end
end % end if
       
    
    
    
    




% -------------------------------------------------------------------------
% Th = 10 RUNS
% -------------------------------------------------------------------------
% iset = 2;
% 
% % SET 1
% if iset == 1
% runNames = {['freeCouetteRun' num2str(1)];...
%     ['freeCouetteRun' num2str(2)];...
%     ['freeCouetteRun' num2str(3)];...
%     ['freeCouetteRun' num2str(4)];...
%     ['freeCouetteRun' num2str(5)];};
% for i = 4 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% elseif iset == 2
% runNames = {['freeCouetteRun' num2str(6)];...
%     ['freeCouetteRun' num2str(7)];...
%     ['freeCouetteRun' num2str(8)];...
%     ['freeCouetteRun' num2str(9)];...
%     ['freeCouetteRun' num2str(10)];};
% for i = 4 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% 
% elseif iset == 3
% runNames = {['freeCouetteRun' num2str(11)];...
%     ['freeCouetteRun' num2str(12)];...
%     ['freeCouetteRun' num2str(13)];...
%     ['freeCouetteRun' num2str(14)];...
%     ['freeCouetteRun' num2str(15)];};
% for i = 1 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% 
% elseif iset == 4
% runNames = {['freeCouetteRun' num2str(16)];...
%     ['freeCouetteRun' num2str(17)];...
%     ['freeCouetteRun' num2str(18)];...
%     ['freeCouetteRun' num2str(19)];...
%     ['freeCouetteRun' num2str(20)];};
% for i = 1 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% 
% elseif iset == 5
% runNames = {['freeCouetteRun' num2str(21)];...
%     ['freeCouetteRun' num2str(22)];...
%     ['freeCouetteRun' num2str(23)];...
%     ['freeCouetteRun' num2str(24)];...
%     ['freeCouetteRun' num2str(25)];};
% for i = 1 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% elseif iset == 6
% runNames = {['freeCouetteRun' num2str(26)];...
%     ['freeCouetteRun' num2str(27)];...
%     ['freeCouetteRun' num2str(28)];...
%     ['freeCouetteRun' num2str(29)];...
%     ['freeCouetteRun' num2str(30)];};
% for i = 1 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% elseif iset == 7
% runNames = {['freeCouetteRun' num2str(31)];...
%     ['freeCouetteRun' num2str(32)];...
%     ['freeCouetteRun' num2str(33)];...
%     ['freeCouetteRun' num2str(34)];...
%     ['freeCouetteRun' num2str(35)];};
% for i = 1 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% elseif iset == 8
% runNames = {['freeCouetteRun' num2str(36)];...
%     ['freeCouetteRun' num2str(37)];...
%     ['freeCouetteRun' num2str(38)];...
%     ['freeCouetteRun' num2str(39)];...
%     ['freeCouetteRun' num2str(40)];};
% for i = 1 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% elseif iset == 9
% runNames = {['freeCouetteRun' num2str(41)];...
%     ['freeCouetteRun' num2str(42)];...
%     ['freeCouetteRun' num2str(43)];...
%     ['freeCouetteRun' num2str(44)];...
%     ['freeCouetteRun' num2str(45)];};
% for i = 1 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% elseif iset == 10
% runNames = {['freeCouetteRun' num2str(46)];...
%     ['freeCouetteRun' num2str(47)];...
%     ['freeCouetteRun' num2str(48)];...
%     ['freeCouetteRun' num2str(49)];...
%     ['freeCouetteRun' num2str(50)];};
% for i = 1 : 5
%   freeCouetteRuns(runNames{i},N,Th,dt,speed,VC,kappa,volFrac(i))
% end
% 
% end
