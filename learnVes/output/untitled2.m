clear; clc;
dt = 1E-5;
nparts = 6;
kappa = 1;

nsampInSets = [];
XstoreIps = [];
XnewStoreIps = [];

for ip = 1 : nparts
disp(ip)
fileName = ['./relaxNewData/n128Dt' num2str(dt) 'RelaxMirrd_part' num2str(ip) '.mat'];
load(fileName)
nsampInSets = [nsampInSets;nInstances];
XstoreIps{ip} = XstandStore;
XnewStoreIps{ip} = XnewStandStore;
end

partIdStart = zeros(nparts,1);
partIdEnd = zeros(nparts,1);

nInstances = sum(nsampInSets);
XstandStore = zeros(2*N,nInstances);
XnewStandStore = XstandStore;

for ip = 1 : nparts
if ip == 1
partIdStart(ip) = 1;
else
partIdStart(ip) = partIdEnd(ip-1)+1;
end
partIdEnd(ip) = sum(nsampInSets(1:ip));
end

for ip = 1 : nparts
XstandStore(:,partIdStart(ip):partIdEnd(ip)) = XstoreIps{ip};
XnewStandStore(:,partIdStart(ip):partIdEnd(ip)) = XnewStoreIps{ip};
end

saveFileName = ['./relaxCombinedData/n128Dt' num2str(dt) 'RelaxMirrdDataSet.mat'];
save(saveFileName,'nInstances','XstandStore','XnewStandStore','N','dt','kappa','-v7.3')
