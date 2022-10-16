clear all
%format long
%fid = fopen('time.log','w');
%fclose(fid);

P = path; ii = find(pwd == filesep); ii = ii(end);
subPath = pwd; subPath = [subPath(1:ii) 'src'];
if isempty(strfind(P,subPath))
    addpath(subPath); 
end;
subPath = pwd; subPath = [subPath(1:ii) 'src/neareval'];
if isempty(strfind(P, subPath))
    addpath(subPath);
end
%Load paths for src directory and neareval directory

directory = '/scratch/quaife/tracerSimulations/';
%Directory where we keep the vesicle simulations
%options.file = 'example6_data.bin';
options.file = 'couetteFlow.bin';
%which file you want to load

if (strcmp(options.file,'example1_data.bin'))
    options.n = 128;
    options.nv = 1;
    options.confined = false;
elseif (strcmp(options.file,'example2_data.bin'))
    options.n = 64;
    options.nv = 1;
    options.confined = false;
elseif (strcmp(options.file,'example3_data.bin'))
    options.n = 32;
    options.nv = 4;
    options.confined = false;
elseif (strcmp(options.file,'example4_data.bin'))
    options.n = 64;
    options.nv = 2;
    options.confined = false;
elseif (strcmp(options.file,'example6_data.bin'))
    options.n = 64;
    options.nv = 9;
    options.confined = false;
elseif (strcmp(options.file,'couetteFlow.bin'))
    options.n = 64;
    options.nv = 120;
    options.M = [256 128];
    options.nholes = 1;
    %Total number of points on two solid walls
    options.confined = true;
end
%need to specify the problem size

if ~options.confined
  [posx,posy,ten,velx,vely,fx,fy,time,ea,el] = loadfile_withdensity...
      ([directory options.file],options.n,options.nv);
  den1 = [];
  den2 = [];
  const = [];
else
  [posx,posy,ten,velx,vely,fx,fy,den1,den2,const,time,ea,el] = ...
      loadfile_confined([directory options.file],options.n,options.nv,...
      sum(options.M),options.nholes);
end
%extract all the data
%   posx is x-location of vesicle
%   posy is y-location of vesicle
%   ten is tension of vesicle
%   velx is the x-component of the velocity on the vesicle
%   vely is the y-component of the velocity on the vesicle
%   fx is the first component of the density function
%   fy is the second component of the density function
%   den1 is the x component of the density function defined on the solid walls
%   den2 is the y component of the density function defined on the solid walls
%   const stores the Rotlet and Stokeslet terms
%   time is the time steps
%   ea is the error in area
%   el is the error in length

dt = time(2)-time(1);
%Get the time step size

[n,nv,ntime] = size(posx);
if (n ~= options.n | nv ~= options.nv)
    disp('PROBLEM')
end
%Make sure that n and nv agree with those specified by the user

if (strcmp(options.file,'example1_data.bin'))
    [ptx pty] = meshgrid(-2:.1:3,-2:.1:2);
    ptx = ptx(:)';
    pty = pty(:)';
elseif (strcmp(options.file,'example2_data.bin'))
    [ptx pty] = meshgrid(-2:.5:2,-3:.5:3);
    ptx = ptx(:)';
    pty = pty(:)';
elseif (strcmp(options.file,'example3_data.bin'))
%    [ptx pty] = meshgrid(-1:0.5:12,-2:0.2:2);
     [ptx pty] = meshgrid(2:.1:3,0:.1:1);
    ptx = ptx(:)';
    pty = pty(:)';
elseif (strcmp(options.file,'example4_data.bin'))
  [r,theta] = meshgrid(.09:.1:.99,0:2*pi/10:2*pi-2*pi/10);
  ptx = (r.*cos(theta)-1.2);
  pty = 2/3*sin(theta);
  %[ptx pty] = meshgrid(-3.05:.1:3.05,-1:0.02:1);
  %[ptx pty] = meshgrid(-3.05:.4:3.05,-1:0.1:1);
  %ptx = [-2.05 -2.05];
  %pty = [-0.2 -0.1];
  %ptx = ptx(1);
  %pty = pty(1);
  ptx = ptx(95);
  pty = pty(95);
  ptx = ptx(:)';
  pty = pty(:)';
elseif (strcmp(options.file,'example6_data.bin'))
  [ptx pty] = meshgrid(0:pi/10:pi,0:pi/10:pi);
  ptx = ptx(:)';
  pty = pty(:)';
elseif (strcmp(options.file,'couetteFlow.bin')) 

  den1(:,1) = den1(:,2);
  den2(:,1) = den2(:,2);
  const(:,1) = const(:,2);
  %Set density function and stokeslet/rotlet terms at time step 1 to the
  %values at time step 2.  At time step 1, they are always 0.

  Ri = 10;
  Ro = 20;
  opts.bd = @(ind,m) sampleBd(ind,m,1,'couette','Ri',Ri,'Ro',Ro);
  bound = @(n) fixedBound(n,opts.bd,1);
  evalFieldDirect(options.M,[],[],bound,[],[],[],[]);
  %Need to call to pass the persistent variables


  ptx = [];
  pty = [];
  r = linspace(10.1,19.9,15);
  Ntheta = 4;
  ra = pi/2;
  ptx = [ptx;r(1)*cos((0:Ntheta-1)*2*pi/Ntheta/20+ra)'];
  pty = [pty;r(1)*sin((0:Ntheta-1)*2*pi/Ntheta/20+ra)'];
  h = 2*pi*r(1)/Ntheta;
  for k=2:numel(r)
    ra = pi/2;
    ptx = [ptx;r(k)*cos((0:Ntheta-1)*2*pi/Ntheta/20+ra)'];
    pty = [pty;r(k)*sin((0:Ntheta-1)*2*pi/Ntheta/20+ra)'];
  end
  %[r,theta] = meshgrid(11.5:0.5:19.5,0:2*pi/20:2*pi-2*pi/20);
  %ptx = r.*cos(theta);
  %pty = r.*sin(theta);
  ptx = ptx(:)';
  pty = pty(:)';
end
ntrac = numel(ptx);
%set of initial conditions for tracers


pt = [ptx pty]';
vel = @(t,x) (velocity1(t,x,time,options,ntrac,nv,posx,posy,fx,fy,...
    den1,den2,const,n));
%Build a function handle that computes the velocity at a point x 
%and a time t


%time = time(1:1001);
euler = 1;
profile on
tic
if (~euler)
  options45 = odeset('RelTol',1e-8,'AbsTol',1e-10);
  [t y] = ode45(vel,time(1:1001),pt,options45);
  %[t y] = ode45(vel,time,pt);
end
if (euler)
  time = time(1:21);


  y = zeros(numel(time),numel(pt));
  y(1,:) = [ptx pty]';
  t = time(1);
  for k=2:numel(time)
    y(k,:) = y(k-1,:) + dt * vel(time(k-1),y(k-1,:)')';
    t = [t;time(k)];
  end


  %y(2,:) = y(1,:) + dt * vel(time(1),y(1,:)')';
  %t = [t;time(2)];
  %First time step is Explicit Euler
  %for k=3:numel(time)
  %  y(k,:) = y(k-1,:) + dt/2 * ...
  %    (vel(time(k-1),y(k-1,:)')' + vel(time(k-2),y(k-2,:)')');
  %    disp(vel(time(k-1),y(k-1,:)'))
  %  t = [t;time(k)];
  %end
  %Use Crank-Nicholson for all other time steps
end


toc
profile off
profsave(profile('info'),'tracersProfile');
save(['tracers_' options.file(1:8)],'t','y');


