function prams = initMLARM(prams)
% Set a path pointing to src directory and set options and
% prams to default values if they are left blank

P = path; ii = find(pwd == filesep); ii = ii(end);
subPath = pwd; subPath = [subPath(1:ii) 'src'];
if isempty(strfind(P, subPath)),addpath(subPath);end
subPath = pwd; subPath = [subPath(1:ii) 'examples'];
if isempty(strfind(P, subPath)),addpath(subPath);end

PramList = {'N','nv','Nbd','nvbd','kappa','NbdInt','NbdExt','nvbdInt',...
    'nvbdExt','dt','Th','fmm','fmmDLP','bgFlow','speed','totnv','nrow',...
    'ncol','Dpostx','Dposty','Dx','Dy','epsilon','gPer','nCircular'};
defaultPram.N = 32;
defaultPram.nv = 1;
defaultPram.Nbd = 0;
defaultPram.NbdInt = 0;
defaultPram.NbdExt = 0;
defaultPram.nvbd = 0;
defaultPram.nvbdInt = 0;
defaultPram.nvbdExt = 0;
defaultPram.Th = 1;
defaultPram.dt = 1e-4;
defaultPram.kappa = 1;
defaultPram.fmm = false;
defaultPram.fmmDLP = false;
defaultPram.bgFlow = 'couette';
defaultPram.speed = 100;
defaultPram.totnv = 1;
defaultPram.nrow = 0;
defaultPram.ncol = 0;
defaultPram.Dpostx = 0;
defaultPram.Dposty = 0;
defaultPram.Dx = 0;
defaultPram.Dy = 0;
defaultPram.epsilon = 0;
defaultPram.gPer = 0;
defaultPram.nCircular = 0;

for k = 1:length(PramList)
  if ~isfield(prams,PramList{k})
    eval(['prams.' PramList{k} '=defaultPram.' PramList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end
end




