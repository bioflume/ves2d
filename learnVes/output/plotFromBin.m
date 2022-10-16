filename = './repulsion/repStr2000_repScale1.5VF30.bin';

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
%%


wallx = [xwalls;xwalls(1,:)];
wally = [ywalls;ywalls(1,:)];
count = 1;
frameFile = ['frames/image'];
for k = 1 : 50 : ntime
  figure(1); clf;
  plot(wallx, wally, 'k', 'linewidth', 2);
  hold on
  vecx = [vesx(:,:,k);vesx(1,:,k)];
  vecy = [vesy(:,:,k);vesy(1,:,k)];
  plot(vecx, vecy, 'r', 'linewidth', 2)
  axis equal
  title(time(k))
  
  figname = [frameFile, sprintf('%04d',count),'.png'];
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
  ylim([-2.5 2.5])
  xlim([-2.5 2.5])
  pause(0.1)
  ax = gca;
  exportgraphics(ax, figname, 'Resolution', 300)
  count = count + 1;
  
    
end

%%
if 1
fileGT = 'nv70N64VF30GT_bgFlowcouette_speed100';
load(fileGT)
XTrue = XhistTrue;
tTrue = timeTrue;

wallx = [Xwalls(1:end/2,:);Xwalls(1,:)];
wally = [Xwalls(end/2+1:end,:);Xwalls(end/2+1,:)];
count = 1;
frameFile = ['frames_true/image'];
for k = 1 : 50 : ntime
  figure(1); clf;
  plot(wallx, wally, 'k', 'linewidth', 2);
  hold on
  vecx = [XTrue(1:end/2,:,k);XTrue(1,:,k)];
  vecy = [XTrue(end/2+1:end,:,k);XTrue(end/2+1,:,k)];
  plot(vecx, vecy, 'r', 'linewidth', 2)
  axis equal
  title(time(k))
  
  figname = [frameFile, sprintf('%04d',count),'.png'];
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
  ylim([-2.5 2.5])
  xlim([-2.5 2.5])
  pause(0.1)
  ax = gca;
  exportgraphics(ax, figname, 'Resolution', 300)
  count = count + 1;
  
    
end
    
    
    
    
end