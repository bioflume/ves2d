mul = 2.^(0:10)';
dts = 1E-5*mul;

for ip = 1 : 8
  fileAdd = ['n256Dt' num2str(dts(ip)) 'RelaxAddDataSet.mat'];
  load(fileAdd)

  XstandAdd = XstandStore;
  XnewStandAdd = XnewStandStore;

  fileAll = ['n256Dt' num2str(dts(ip)) 'RelaxAllDataSet.mat'];
  load(fileAll)

  XstandStore = [XstandStore XstandAdd];
  XnewStandStore = [XnewStandStore XnewStandAdd];

  nInstances = size(XnewStandStore,2);

  fileNew = ['n128Dt' num2str(dts(ip)) 'RelaxDataSet.mat'];
  save(fileNew,'nInstances','XstandStore','XnewStandStore','N','dt','kappa','-v7.3')
end