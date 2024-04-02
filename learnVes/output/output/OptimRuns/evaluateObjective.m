function evaluateObjective(iter)

load(['lowFid_Iter' num2str(iter) '.mat'])

deltaCys = finalCys - initCys;
DelYscaling = 4*(Dpostx+Dx)*atan(1/6)/2+Dposty;
nvesTotal = numel(rem_ids)+nZigZags;
Yend = XhistStore{it}(end/2+1:end,1);
Feval = -nZigZags/nvesTotal+mean(deltaCys)/DelYscaling-abs(mean(Yend)+0.5*Dy)/(0.5*Dy);

disp(['nZZ/nTot = ' num2str(nZigZags/nvesTotal)])
disp(['delY/scaling = ' num2str(mean(deltaCys)/DelYscaling)])
disp(['meanY/Dy = ' num2str(abs(mean(Yend)+0.5*Dy)/(0.5*Dy))])
disp(['new Feval = ' num2str(Feval)])

