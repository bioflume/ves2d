clear;
kappa = 1;
idx_dt = [1:1:11]';
npar = [32*ones(6,1); 128; 256; 512; 1024; 2048];
dts = 1E-6*2.^(0:15)';

for jj = 1 : numel(idx_dt)
  dtEffect = dts(jj);   
  Xstore = [];
  XnewStore = [];
  nsampInSets = [];
  disp(['DtEffect id: ' num2str(jj)])
  for iset = 1 : npar(jj)
    disp(['Set id: ' num2str(iset)])
    fileName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/workspace/n256Dt' num2str(dtEffect) 'Relax100K_part' num2str(iset) '.mat'];
    load(fileName)
    nsampInSets(iset) = nInstances;
    Xstore = [Xstore XstandStore];
    XnewStore = [XnewStore XnewStandStore];
  end
  fileName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/workspace/completeData/n256Dt' num2sstr(dtEffect) 'RelaxAllDataSet.mat'];
  nInstances = sum(nsampInSets);
  XstandStore = Xstore;
  XnewStandStore = XnewStore;
  dt = dtEffect;
  save(fileName, 'nInstances','XstandStore','XnewStandStore','N','dt','kappa','-v7.3')
end

