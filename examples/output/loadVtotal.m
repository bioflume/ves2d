function [VselfDecoup,VnearDecoup,VfarDecoup,Vinf] = loadVtotal(runName)

% Load the simulation data
[posx,posy,~,~,~,~,~,~,time,N,nv] = loadFile([runName 'Data.bin']);
ntime = numel(time);

% Load the velocity data
fid = fopen([runName '_vTotalDecoup.bin'],'r'); % DECOUPLED TENSION
% fid = fopen([runName '_vTotalPreComp.bin'],'r'); % FIXED TENSION

valDecoup = fread(fid,'double');
ntimeVel = numel(valDecoup)/(4*nv*2*N); % this should be ntime-1

% fid = fopen([runName '_vTotal.bin'],'r');
% valDecoup = fread(fid,'double');

% Vself = zeros(2*N,nv,ntimeVel);
% Vnear = Vself;
% Vfar = Vself;
% Vinf = Vself;

VselfDecoup = zeros(2*N,nv,ntimeVel);
VnearDecoup = zeros(2*N,nv,ntimeVel);
VfarDecoup = zeros(2*N,nv,ntimeVel);
Vinf = zeros(2*N,nv,ntimeVel);

ives = 5;
idebug = false;

for t = 1 : ntimeVel
%   vel = val((t-1)*nv*2*N*4+1:t*nv*2*N*4);
%   Vself(:,:,t) = reshape(vel(1:2*N*nv),2*N,nv); 
%   Vnear(:,:,t) = reshape(vel(2*N*nv+1:2*(2*N*nv)),2*N,nv);
%   Vfar(:,:,t) = reshape(vel(2*(2*N*nv)+1:3*(2*N*nv)),2*N,nv);
%   Vinf(:,:,t) = reshape(vel(3*(2*N*nv)+1:4*(2*N*nv)),2*N,nv);
  
  vel = valDecoup((t-1)*nv*2*N*4+1:t*nv*2*N*4);
  VselfDecoup(:,:,t) = reshape(vel(1:2*N*nv),2*N,nv); 
  VnearDecoup(:,:,t) = reshape(vel(2*N*nv+1:2*(2*N*nv)),2*N,nv);
  VfarDecoup(:,:,t) = reshape(vel(2*(2*N*nv)+1:3*(2*N*nv)),2*N,nv);
  Vinf(:,:,t) = reshape(vel(3*(2*N*nv)+1:4*(2*N*nv)),2*N,nv);
  
  if idebug
  figure(1);clf;hold on;
  plot(posx(:,:,t),posy(:,:,t),'r')
  axis([-21 21 -21 21])
  axis equal

%   Vtotal = Vself+Vnear+Vfar+Vinf;
  VtotalDecoup = VselfDecoup+VnearDecoup+VfarDecoup+Vinf;
  figure(2); clf; hold on;
  plot(posx(:,ives,t),posy(:,ives,t),'r')
%   quiver(posx(:,ives,t),posy(:,ives,t),Vtotal(1:end/2,ives,t),Vtotal(end/2+1:end,ives,t))
  quiver(posx(:,ives,t),posy(:,ives,t),VtotalDecoup(1:end/2,ives,t),VtotalDecoup(end/2+1:end,ives,t))
  axis equal
  
  disp('SLP')
%   norm(Vself(:,ives,t)+Vnear(:,ives,t)+Vfar(:,ives,t))
  disp('Decoupled SLP')
  norm(VselfDecoup(:,ives,t)+VnearDecoup(:,ives,t)+VfarDecoup(:,ives,t))
  disp('Vinf')
  norm(Vinf(:,ives,t))
  
  pause(0.1)
  end
end




