directory = '/scratch/quaife/tracerSimulations/';
file = 'example6_data.bin';
%file = 'couetteFlow.bin';


if (strcmp(file,'example1_data.bin'))
    options.n = 128;
    options.nv = 1;
    options.confined = false;
    load tracers_example1.mat;
elseif (strcmp(file,'example2_data.bin'))
    options.n = 64;
    options.nv = 1;
    options.confined = false;
    load tracers_example2.mat;
elseif (strcmp(file,'example3_data.bin'))
    options.n = 32;
    options.nv = 4;
    options.confined = false;
    load tracers_example3.mat;
elseif (strcmp(file,'example4_data.bin'))
    options.n = 64;
    options.nv = 2;
    options.confined = false;
    load tracers_example4.mat;
elseif (strcmp(file,'example6_data.bin'))
    options.n = 64;
    options.nv = 9;
    options.confined = false;
    load tracers_example6.mat;
elseif (strcmp(file,'couetteFlow.bin'))
    options.n = 64;
    options.nv = 120;
    options.M = [256 128];
    options.nholes = 1;
    %Total number of points on two solid walls
    options.confined = true;
    load tracers_couetteF.mat;
end
%need to specify the problem size


if (options.confined)
  [posx,posy,ten,velx,vely,fx,fy,den1,den2,const,time,ea,el] = ...
    loadfile_confined([directory file],options.n,options.nv,...
    sum(options.M),options.nholes);
else
  [posx,posy,ten,velx,vely,fx,fy,time,ea,el] = ...
    loadfile_withdensity([directory file],options.n,options.nv);
end



[ntime ntrac]=size(y);
ntrac = ntrac/2;
xx = y(:,1:ntrac);
yy = y(:,ntrac+1:2*ntrac);

p = 2;

I = find(time <= t(end) + 1e-12);
%find index of times in vesicle simulation that are less
%than the final time of the tracer simulation
%Add on a buffer as time tends to have non-zero entries in the 1e-16 decimals

X = zeros(numel(I),ntrac);
Y = zeros(numel(I),ntrac);
for i = 1:numel(I) 
  %disp(numel(I)-i)
  tau = time(i);
  %tau is the time step of the vesicles
  [m,j] = min(abs(tau-t));
  %find the closest time step of the tracers to the current vesicle time

  if (m < 1e-12)
    X(i,:) = xx(j,:);
    Y(i,:) = yy(j,:);
    %If the two time meshes agree at points, don't need to do interpolation
  else
    j = max(p+1,j);
    %make sure we have enough data points below tau to do the interpolation
    if (j+p > numel(t) && m > 1e-12)
      j = j - p;
    end
    tint = t(j-p:j+p);
    for k=1:ntrac
      xint = xx(j-p:j+p,k);
      yint = yy(j-p:j+p,k);
      px = polyfit(tint,xint,2*p-2);
      py = polyfit(tint,yint,2*p-2);
      X(i,k) = polyval(px,tau);
      Y(i,k) = polyval(py,tau);
    end
  end
end


