function submitDataPrep(idt,npar)
%mul = 2.^(0:10)';
%dts = 1E-5*mul;
dts = [1E-6;1E-5;1E-4];
p = parcluster('local');
p.NumWorkers = npar;
%parpool(npar);
%parfor irun = 1 : npar
%  prepareRelaxFlowMapData(dts(idt),irun,npar)
%end
job = [];
for irun = 1 : npar
  disp(['submitting ' num2str(irun)])
  job{irun} = batch(p,@prepareRelaxFlowMapData,0,{dts(idt),irun,npar});
end
for j = 1 : npar
  wait(job{j},'finished')
  disp(['finished ' num2str(j)])
end
delete(findJob(p))
end
