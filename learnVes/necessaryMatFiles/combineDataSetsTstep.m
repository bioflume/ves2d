clear;
load n256Dt2p5E3Kb1Relax100K_part1.mat
kappa = 1;

nsampInSets(1) = nInstances;
Xstore1 = XstandStore;
XnewStore1 = XnewStandStore;

load n256Dt2p5E3Kb1Relax100K_part2.mat
nsampInSets(2) = nInstances;
Xstore2 = XstandStore;
XnewStore2 = XnewStandStore;

load n256Dt2p5E3Kb1Relax100K_part3.mat
nsampInSets(3) = nInstances;
Xstore3 = XstandStore;
XnewStore3 = XnewStandStore;

load n256Dt2p5E3Kb1Relax100K_part4.mat
nsampInSets(4) = nInstances;
Xstore4 = XstandStore;
XnewStore4 = XnewStandStore;


nInstances = sum(nsampInSets);
XstandStore = zeros(2*N,nInstances);
XnewStandStore = XstandStore;

part1s = 1; part1e = nsampInSets(1);
part2s = nsampInSets(1)+1; part2e = sum(nsampInSets(1:2));
part3s = part2e+1; part3e = sum(nsampInSets(1:3));
part4s = part3e+1; part4e = nInstances;

XstandStore(:,part1s:part1e) = Xstore1;
XstandStore(:,part2s:part2e) = Xstore2;
XstandStore(:,part3s:part3e) = Xstore3;
XstandStore(:,part4s:part4e) = Xstore4;

XnewStandStore(:,part1s:part1e) = XnewStore1;
XnewStandStore(:,part2s:part2e) = XnewStore2;
XnewStandStore(:,part3s:part3e) = XnewStore3;
XnewStandStore(:,part4s:part4e) = XnewStore4;

save('n256Dt2p5E3Kb1RelaxAllDataSet.mat','nInstances','XstandStore',...
  'XnewStandStore','N','dt','kappa','-v7.3')
