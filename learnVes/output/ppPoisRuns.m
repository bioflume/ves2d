clear; clc;
runs{1} = 'poisDNNnewSingVes_speed8000_newNet_exactAdv.bin';
runs{2} = 'poisDNNnewSingVes_speed12000_newNet_exactAdv.bin';
runs{3} = 'poisDNNnewSingVes_speed16000_newNet_exactAdv.bin';
runs{4} = 'poisDNNnewSingVes_speed24000_newNet_exactAdv.bin';
runs{5} = 'poisDNNnewSingVes_speed30000_newNet_exactAdv.bin';


for irun = 1 : 5
  [vesxN, vesyN, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runs{irun});
  times{irun} = time;
  for k = 1 : numel(time)
    centxs{irun}(k,1) = mean(vesxN(:,k));
    centys{irun}(k,1) = mean(vesyN(:,k));
  end
end
disp('Centers calculated')
    

% Plot the equilibrium shapes

equilXs = zeros(128,2,3);
equilTimes = [0.20; 0.12; 0.12; 0.049; 0.049];
for irun = 1 : 5
[vesxN, vesyN, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runs{irun});
id = find(time>=equilTimes(irun));
id = id(1);
equilXs(:,1,irun) = vesxN(:,id);
equilXs(:,2,irun) = vesyN(:,id);
end
%%
figure(1); clf;
cxs = [0; 0.5; 1; 1.5; 2];
for irun = 1 : 5
vecx = equilXs(:,1,irun)-mean(equilXs(:,1,irun))+cxs(irun);
vecy = equilXs(:,2,irun)-mean(equilXs(:,2,irun));
plot([vecx;vecx(1)], [vecy;vecy(1)],'linewidth',2)
hold on
end

axis equal