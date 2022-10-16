function [xImp,yImp,xExp,yExp,time,n,nv] = ...
    loadFileExplicit(runName,imovie,iExplicit)
% [xImp,yImp,xExp,yExp,time,n,nv] = loadFileExplicit(runName) extracts the
% data of the run with runName. 
% xImp = [number of points per vesicle, number of vesicles, number of time
% steps] is the format of the matrices. 
%
% xImp, yImp are the x and y coordinates of the vesicles from the
% simulation where we use the implicit time stepping.
%
% xExp, yExp are the counterparts of xImp and yImp, here we use the
% explicit time stepping with the implicit solution from the previous time
% step.
%
% n is the number of points per vesicle
% nv is the number of vesicles
% time is the vector which keeps the times at which we have a solution


file = [runName '_Data.bin'];

fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);
n = val(1);
nv = val(2);
nbd = val(3);
nvbd = val(4);
walls = val(5:5+2*nbd*nvbd-1);
val = val(5+2*nbd*nvbd:end);

ntime = numel(val)/(3*n*nv+3);
if ntime ~= ceil(ntime);
  disp('PROBLEM WITH VALUES FOR n AND nv');
end

wallx = zeros(nbd,nvbd);
wally = zeros(nbd,nvbd);
for k = 1:nvbd
  istart = (k-1)*nbd+1;
  iend = k*nbd;
  wallx(:,k) = walls(istart:iend);
  wally(:,k) = walls((istart:iend)+nbd*nvbd);
end
%wallx = walls(1:nbd*nvbd);
%wally = walls(nbd*nvbd+1:2*nbd*nvbd);
xImp = zeros(n,nv,ntime);
yImp = zeros(n,nv,ntime);
ten = zeros(n,nv,ntime);
time = zeros(ntime,1);
ea = zeros(ntime,1);
el = zeros(ntime,1);

istart = 1;
for m = 1:ntime
  for k=1:nv
    iend = istart + n - 1;
    xImp(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load x positions

  for k=1:nv
    iend = istart + n - 1;
    yImp(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load y positions

  for k=1:nv
    iend = istart + n - 1;
    ten(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load tensions

  ea(m) = val(istart);
  el(m) = val(istart+1);
  time(m) = val(istart+2);
  istart = istart + 3;
end

% Load the explicit solution
if iExplicit
fileExp = [runName '_Explicit_Data.bin'];
fid = fopen(fileExp,'r');
val = fread(fid,'double');
N = val(1);
nv = val(2);
ntime = val(3);
X = reshape(val(4:end),2*N,nv,ntime);
xExp = zeros(N,nv,ntime+1);
yExp = zeros(N,nv,ntime+1);
xExp(:,:,1) = xImp(:,:,1);
yExp(:,:,1) = yImp(:,:,1);
xExp(:,:,2:end) = X(1:end/2,:,:);
yExp(:,:,2:end) = X(end/2+1:end,:,:);
else
xExp = [];
yExp = [];
end
    
if imovie
  figure(1); clf; 
%   axs = [0 2*pi 0 2*pi] + [-0.5 0.5 -0.5 0.5];
  for k = 1 :10: ntime
    plot([xImp(:,:,k);xImp(1,:,k)],[yImp(:,:,k);yImp(1,:,k)],'r','linewidth',2);
    hold on
    if iExplicit
    plot(xExp(:,:,k),yExp(:,:,k),'b','linewidth',2);
    end
%     axis(axs)
    axis equal
    title(['t = ' num2str(time(k))])
    hold off
    pause(0.1)
  end
    
end

end


