function [posx,posy,ten,velx,vely,fx,fy,den1,den2,const,...
    time,ea,el] = loadfile_confined(file,n,nv,M,nholes)
%M is the number of points per hole
%nholes is the number of holes

fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);
data_per_timestep = 7*n*nv + 3;
%For everything but the density functions on the boundary
data_per_timestep = data_per_timestep + 2*M + 3*nholes;
%For density function (x and y coordinate) and the 3 constants per hole
ntime = numel(val)/data_per_timestep;
if (ntime ~= ceil(ntime))
  disp('PROBLEM WITH VALUES FOR n AND nv')
end

posx = zeros(n,nv,ntime);
posy = zeros(n,nv,ntime);
ten = zeros(n,nv,ntime);
velx = zeros(n,nv,ntime);
vely = zeros(n,nv,ntime);
fx = zeros(n,nv,ntime);
fy = zeros(n,nv,ntime);
time = zeros(ntime,1);
ea = zeros(ntime,1);
el = zeros(ntime,1);
den1 = zeros(M,ntime);
den2 = zeros(M,ntime);
const = zeros(3,nholes,ntime);

istart = 0;
for kk=1:ntime
    istart = istart + 1;
    for j=1:nv
        iend = istart + n - 1;
        posx(:,j,kk) = val(istart:iend);
        istart = iend+1;
        iend = istart + n - 1;
        posy(:,j,kk) = val(istart:iend);
        istart = iend+1;
    end
   
    for j=1:nv
        iend = istart + n - 1;
        ten(:,j,kk) = val(istart:iend);
        istart = iend + 1;
    end

    for j=1:nv
        iend = istart + n - 1;
        velx(:,j,kk) = val(istart:iend);
        istart = iend+1;
        iend = istart + n - 1;
        vely(:,j,kk) = val(istart:iend);
        istart = iend+1;
    end

    iend = istart + M-1;
    den1(:,kk) = val(istart:iend);
    istart = iend+1;
    iend = istart + M-1;
    den2(:,kk) = val(istart:iend);
    istart = iend+1;

    for j=1:nholes
      iend = istart + 3 - 1;
      const(:,j,kk) = val(istart:iend);
      istart = iend+1;
    end

    for j=1:nv
        iend = istart + n -1;
        fx(:,j,kk) = val(istart:iend);
        istart = iend+1;
        iend = istart+n-1;
        fy(:,j,kk) = val(istart:iend);
        istart = iend+1;
    end
    
    time(kk) = val(istart);
    ea(kk) = val(istart+1);
    el(kk) = val(istart+2);  
    istart = istart + 2;
end
