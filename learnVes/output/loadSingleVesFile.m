function [vesx, vesy, ten, time, N, nv, xinit, yinit, ncountNN, ncountExact] = loadSingleVesFile(filename)

fid = fopen(filename,'r');
val = fread(fid,'double');
fclose(fid);
N = val(1);
nv = val(2);
Xinit = val(3:3+2*N*nv-1);
xinit = zeros(N,nv); yinit = zeros(N,nv);
istart = 1;

for iv = 1 : nv
  iend = istart + N - 1;
  xinit(:,iv) = Xinit(istart:iend);
  istart = iend + 1;
end
for iv = 1 : nv
  iend = istart + N - 1;
  yinit(:,iv) = Xinit(istart:iend);
  istart = iend + 1;
end
val = val(3+2*N*nv:end);

ntime = floor(numel(val)/(3*N*nv+3));

vesx = zeros(N,nv,ntime);
vesy = zeros(N,nv,ntime);
ten = zeros(N,nv,ntime);
time = zeros(ntime,1);
ncountNN = zeros(ntime,1);
ncountExact = zeros(ntime,1);

istart = 1;
for it = 1 : ntime
  time(it) = val(istart);
  ncountNN(it) = val(istart+1);
  ncountExact(it) = val(istart+2);
  istart = istart+3;
  for iv = 1 : nv
    iend = istart + N - 1; 
    vesx(:,iv,it) = val(istart:iend);
    istart = iend + 1;
  end
  
  for iv = 1 : nv
    iend = istart + N - 1; 
    vesy(:,iv,it) = val(istart:iend);
    istart = iend + 1;
  end
  
  for iv = 1 : nv
    iend = istart + N - 1; 
    ten(:,iv,it) = val(istart:iend);
    istart = iend + 1;
  end  
end


end