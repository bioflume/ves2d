function [posx,posy,wallx,wally,extWallx,extWally,intWallx,intWally,...
    time,N,nv] = loadFile(file)

fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);
N = val(1);
nv = val(2);
Nbd = val(3);
nvbd = val(4);
NbdExt = val(5);
NbdInt = val(6);
nvbdExt = val(7);
nvbdInt = val(8);

% do we have two different sets of walls with different discretization?
diffDiscWalls = false;
if NbdExt ~= 0
  diffDiscWalls = true;
end

if ~diffDiscWalls
  walls = val(9:9+2*Nbd*nvbd-1);
  val = val(9+2*Nbd*nvbd:end);
else
  wallsExt = val(9:9+2*NbdExt*nvbdExt-1);
  wallsInt = val(9+2*NbdExt*nvbdExt:9+2*NbdExt*nvbdExt+2*NbdInt*nvbdInt-1);
  val = val(9+2*NbdExt*nvbdExt+2*NbdInt*nvbdInt:end);
end

% Every time step we save X(size:2*N*nv), sigma(size:N*nv), time(size:1)
ntime = numel(val)/(3*N*nv+1);
if ntime ~= ceil(ntime)
  disp('PROBLEM WITH VALUES FOR n AND nv');
end

if~diffDiscWalls
  wallx = zeros(Nbd,nvbd);
  wally = zeros(Nbd,nvbd);
  for k = 1:nvbd
    istart = (k-1)*Nbd+1;
    iend = k*Nbd;
    wallx(:,k) = walls(istart:iend);
    wally(:,k) = walls((istart:iend)+Nbd*nvbd);
  end
  extWallx = []; extWally = [];
  intWallx = []; intWally = [];
else
  extWallx = zeros(NbdExt,nvbdExt); extWally = zeros(NbdExt,nvbdExt);
  for k = 1 : nvbdExt
    istart = (k-1)*NbdExt+1;
    iend = k*NbdExt;
    extWallx(:,k) = wallsExt(istart:iend);
    extWally(:,k) = wallsExt((istart:iend)+NbdExt*nvbdExt);
  end
  intWallx = zeros(NbdInt,nvbdInt); intWally = zeros(NbdInt,nvbdInt);
  for k = 1 : nvbdInt
    istart = (k-1)*NbdInt+1;
    iend = k*NbdInt;
    intWallx(:,k) = wallsInt(istart:iend);
    intWally(:,k) = wallsInt((istart:iend)+NbdInt*nvbdInt);
  end
  wallx = []; wally = [];
end

posx = zeros(N,nv,ntime);
posy = zeros(N,nv,ntime);
ten = zeros(N,nv,ntime);
time = zeros(ntime,1);

istart = 1;
for m = 1:ntime
  for k=1:nv
    iend = istart + N - 1;
    posx(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load x positions

  for k=1:nv
    iend = istart + N - 1;
    posy(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load y positions

  for k=1:nv
    iend = istart + N - 1;
    ten(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load tensions

  time(m) = val(istart);
  istart = istart + 1;
end



