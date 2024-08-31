npar = 4;
p = parcluster('local');
p.NumWorkers = npar;

for i = 1 : npar
  prepareVelTrainFFTData(i)
end

parfor i = 1 : npar
prepareSelfTensionData(i,npar);
end

for i = 1 : npar
prepareAdvectTensionData(i,npar)
end

for i = 1 : npar
prepareNearFieldVelocityData(i,npar)
end

