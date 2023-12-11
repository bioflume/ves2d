function submitDataGenRuns(nv,idStart,runNameInit)

folderName = './output/';
speed = 10000;

npar = 5;
p = parcluster('local');
p.NumWorkers = npar;

job = [];
for irun = 1 : npar
  id = idStart - 1 + irun;
  runName = [runNameInit 'Id' num2str(id)];
  disp(['Submitting ' num2str(id)])
  job{irun} = batch(p,@generateVesShapesFlow, 0, {runName, folderName, nv, speed});
end

for j = 1 : npar
  wait(job{j},'finished')
end
delete(findJob(p))
end
