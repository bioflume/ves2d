clear;

load ./selfTensionData/selfTensionData_1.mat
nsampInSets(1) = nsampInSet;
Xinput1 = XstandStore;
Toutput1 = selfTenStore;

load ./selfTensionData/selfTensionData_2.mat
nsampInSets(2) = nsampInSet;
Xinput2 = XstandStore;
Toutput2 = selfTenStore;

load ./selfTensionData/selfTensionData_3.mat
nsampInSets(3) = nsampInSet;
Xinput3 = XstandStore;
Toutput3 = selfTenStore;

load ./selfTensionData/selfTensionData_4.mat
nsampInSets(4) = nsampInSet;
Xinput4 = XstandStore;
Toutput4 = selfTenStore;

part1s = 1; part1e = nsampInSets(1);
part2s = nsampInSets(1)+1; part2e = sum(nsampInSets(1:2));
part3s = part2e+1; part3e = sum(nsampInSets(1:3));
part4s = part3e+1; part4e = nInstances;

Xinput = zeros(256,nInstances);
Toutput = zeros(128,nInstances);

Xinput(:,part1s:part1e) = Xinput1;
Xinput(:,part2s:part2e) = Xinput2;
Xinput(:,part3s:part3e) = Xinput3;
Xinput(:,part4s:part4e) = Xinput4;

Toutput(:,part1s:part1e) = Toutput1;
Toutput(:,part2s:part2e) = Toutput2;
Toutput(:,part3s:part3e) = Toutput3;
Toutput(:,part4s:part4e) = Toutput4;


save('selfTensionDataSet.mat','nInstances','Xinput','Toutput','-v7.3')
