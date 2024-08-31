function submitDataPrepNearField
npar = 6;
p = parcluster('local');
p.NumWorkers = npar;

job = [];
for irun = 1 : npar
  disp(['submitting ' num2str(irun)])
  job{irun} = batch(p,@prepareNearFieldStokesletData,0,{irun,npar});
end
for j = 1 : npar
  wait(job{j},'finished')
  disp(['finished ' num2str(j)])
end
delete(findJob(p))
end
