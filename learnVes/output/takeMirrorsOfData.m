load n128Dt1e-05RelaxMirrdDataSet.mat

XstandDup = zeros(2*N,nInstances);
XnewStandDup = zeros(2*N,nInstances);

for k = 1 : nInstances
  XstandDup(1:end/2,k) = -XstandStore(1:end/2,k);
  XstandDup(end/2+1:end,k) = XstandStore(end/2+1:end,k);
  XnewStandDup(1:end/2,k) = -XnewStandStore(1:end/2,k);
  XnewStandDup(end/2+1:end,k) = XnewStandStore(end/2+1:end,k);
  
  % figure(1); clf;
  % plot(XstandDup(1:end/2,k),XstandDup(end/2+1:end,k),'linewidth',2)
  % hold on
  % plot(XstandStore(1:end/2,k),XstandStore(end/2+1:end,k),'linewidth',2)
  % axis equal
  % 
  % figure(2); clf;
  % plot(XnewStandDup(1:end/2,k),XnewStandDup(end/2+1:end,k),'linewidth',2)
  % hold on
  % plot(XnewStandStore(1:end/2,k),XnewStandStore(end/2+1:end,k),'linewidth',2)
  % axis equal
  % 
  % pause
  
end

XstandStore = [XstandStore XstandDup];
XnewStandStore = [XnewStandStore XnewStandDup];
nInstances = 2*nInstances;

save n128Dt1e-05RelaxMirrdDataSet.mat XstandStore XnewStandStore nInstances