function [vesx, vesy, ten, etax, etay, time, N, nv, xinit, yinit, xwalls, ywalls, ncountNN, ncountExact] = loadManyVesFile(filename)

fid = fopen(filename,'r');
val = fread(fid,'double');
fclose(fid);
N = val(1);
nv = val(2);
Nbd = val(3);
nvbd = val(4);
Xinit = val(5:5+2*(N*nv+Nbd*nvbd)-1);
xinit = zeros(N,nv); yinit = zeros(N,nv);
xwalls = zeros(Nbd,nvbd); ywalls = zeros(Nbd, nvbd);

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
for iv = 1 : nvbd
  iend = istart + Nbd - 1;
  xwalls(:,iv) = Xinit(istart:iend);
  istart = iend + 1;
end
for iv = 1 : nvbd
  iend = istart + Nbd - 1;
  ywalls(:,iv) = Xinit(istart:iend);
  istart = iend + 1;
end

val = val(5+2*(N*nv+Nbd*nvbd):end);

ntime = floor(numel(val)/(3*N*nv+2*Nbd*nvbd+3*nvbd+3));

vesx = zeros(N,nv,ntime);
vesy = zeros(N,nv,ntime);
ten = zeros(N,nv,ntime);
etax = zeros(Nbd,nvbd,ntime);
etay = zeros(Nbd,nvbd,ntime);
RS = zeros(3,nvbd,ntime);

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
  
  for iv = 1 : nvbd 
    iend = istart + Nbd - 1;
    etax(:,iv,it) = val(istart:iend);
    istart = iend + 1;
  end
  
  for iv = 1 : nvbd 
    iend = istart + Nbd - 1;
    etay(:,iv,it) = val(istart:iend);
    istart = iend + 1;
  end
  
  for iv = 1 : nvbd 
    iend = istart + 3 - 1;
    RS(:,iv,it) = val(istart:iend);
    istart = iend + 1;
  end
  
end
end
