function [posx,posy,ten,velx,vely,fx,fy,time,ea,el] = loadfile_withdensity(file,n,nv)

fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);

ntime = numel(val)/(7*n*nv+3);
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

for kk=1:ntime
    istart = (7*n*nv+3)*(kk-1)+1;
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
    
end
