clear; clc;

theta = [atan(1/6);atan(1/8);atan(1/10)];
VCs = [1 10];
repulsion = true;
fmm = true;

runNames{1} = 'RA95_per6';
runNames{2} = 'RA95_per8';
runNames{3} = 'RA95_per10';

for i = 1 : 3
RA_fullDLDPutBack(runNames{i},0.95,theta(i),VCs,repulsion,fmm)
end
