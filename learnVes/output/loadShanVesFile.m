function [vesx, vesy, time, N, nv, xinit, yinit] = loadShanVesFile(filename)

fid = fopen(filename,'r');
val = fread(fid,'double');
fclose(fid);
N = val(1);
nv = val(2);
Xinit = val(3:3+2*(N*nv)-1);
xinit = zeros(N,nv); yinit = zeros(N,nv);

istart = 1;
for iv = 1 : nv
  iend = istart + N - 1;
  xinit(:,iv) = Xinit(istart:iend);
  istart = iend + 1;
  iend = istart + N - 1;
  yinit(:,iv) = Xinit(istart:iend);
  istart = iend + 1;
end


val = val(3+2*(N*nv):end);

ntime = floor(numel(val)/(2*N*nv+1));

vesx = zeros(N,nv,ntime);
vesy = zeros(N,nv,ntime);
time = zeros(ntime,1);

istart = 1;
for it = 1 : ntime
  time(it) = val(istart);
  istart = istart+1;
  for iv = 1 : nv
    iend = istart + N - 1; 
    vesx(:,iv,it) = val(istart:iend);
    istart = iend + 1;
    iend = istart + N - 1; 
    vesy(:,iv,it) = val(istart:iend);
    istart = iend + 1;
  end
  
end
end
