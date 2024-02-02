clear; clc;

filename = 'Newnv10N32n10.bin';
[vesxN, vesyN, tenN, etaxN, etayN, timeN, N, nv, xinit, yinit, xwalls, ywalls, ncountNN, ncountExact] = loadManyVesFile(filename);

filename = 'Oldnv10N32n10.bin';
[vesxO, vesyO, tenO, etaxO, etayO, timeO, N, nv, xinit, yinit, xwalls, ywalls, ncountNN, ncountExact] = loadManyVesFile(filename);

ntime = min(numel(timeN),numel(timeO));

for k = 1 :10: ntime

figure(1);clf;
plot(xwalls,ywalls,'k','linewidth',2)
hold on
plot(vesxN(:,:,k),vesyN(:,:,k),'r','linewidth',2);
plot(vesxO(:,:,k),vesyO(:,:,k),'b','linewidth',2);

axis equal
title(['Time = ' num2str(timeO(k))])
pause(0.1)

end