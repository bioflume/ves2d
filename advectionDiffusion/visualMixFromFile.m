function visualMixFromFile(fileName)
set(0,'DefaultAxesFontSize',22)
set(0,'DefaultAxesFontName', 'Computer Modern')

% Read File
fid = fopen(fileName,'r');
val = fread(fid,'double');
fclose(fid);
ntime  = val(1);
Th     = val(2);
H1Norm = val(3:end);
if numel(H1Norm) ~= ntime
    disp('PROBLEM with SIZE of H1Norm')
end

% Compute normalized mixing
mixing = H1Norm/H1Norm(1);

% Generate Time Domain
deltaT = Th/(ntime-1);
currTime = (numel(H1Norm)-1)*deltaT;
time = linspace(0,currTime,numel(H1Norm))';

plot(time,mixing,'linewidth',4)
xlim([0 currTime])
ylim([0 1])

xlabel('Time')
ylabel('Mixing Ratio')

end
