function [V,mmean,timeAveRMS,spatAv,spatRMS] = WeakStatEq(velFile,Nr,Nt)

set(0,'DefaultAxesFontSize',20)
set(0,'DefaultAxesFontName', 'Helvetica')

% V1 = @(x,y,t)  sin(2*pi*10/200*(t)) * 1/3 * y .* (1-400./(x.^2 + y.^2)); % x-velocity
% V2 = @(x,y,t) -sin(2*pi*10/200*(t)) * 1/3 * x .* (1-400./(x.^2 + y.^2)); % y-velocity

%  V1 = @(x,y,t)  y/3 .* (1-400./(x.^2+y.^2)); % x-velocity
%  V2 = @(x,y,t) -x/3 .* (1-400./(x.^2+y.^2)); % y-velocity


fid = fopen(velFile,'r');
val = fread(fid,2,'double');
Nrow  = val(1);
ntime = val(2);
% Nrow = 128*513*2;
% ntime = 13001;
% Generate grid and compute spatial average
[x,y,r1D,theta,~,~,~] = generateGrid(Nt,Nr,[10 20]);
time = linspace(0,(ntime-1)*0.01,ntime);
spatAv = zeros(ntime,1);
spatRMS = zeros(ntime,1);
rMat = repmat(r1D(1:end),1,Nt);

for tt = 1 : ntime
    vel = fread(fid,Nrow,'double');

%     VxAnt = V1(x,y,time(tt));
%     VyAnt = V2(x,y,time(tt));
    
    % Compute magnitude of velocity
    vx = vel(1:end/2);
    vy = vel(end/2+1:end);
    
%     Vx = reshape(vx,Nr+1,Nt+1);
%     Vy = reshape(vy,Nr+1,Nt+1);

%     Vx0 = Vx - VxAnt;
%     Vy0 = Vy - VyAnt;
    
%     V = (Vx0.^2 + Vy0.^2).^0.5;
    
    V  = (vx.^2 + vy.^2).^0.5;  
    V = reshape(V,Nr+1,Nt+1);
    V = V(:,1:end-1);
 
%     plot(x(:,1),Vy0(:,1),'linewidth',4)
%     pause

    spatAv(tt) = 1/(pi*(r1D(end).^2-r1D(1).^2))*trapz(theta,trapz(r1D,rMat.*V(:,:),1));
    spatRMS(tt) = sqrt(1/(pi*(r1D(end).^2-r1D(1).^2))*trapz(theta,trapz(r1D,rMat.*(V(:,:)).^2,1)));
    disp(['Computing mean and RMS, step: ' num2str(tt) ' of ' num2str(ntime)])
end

% Compute mean
mmean = zeros(ntime-1001,1);
timeAveRMS = zeros(ntime-1001,1);
for i = 1 : ntime-1001
    mmean(i) = 1/(time(i+1001)-time(i)) * trapz(time(i:i+1001),spatAv(i:i+1001));
    timeAveRMS(i) = 1/(time(i+1001)-time(i)) * trapz(time(i:i+1001),spatRMS(i:i+1001));
end

% subplot(1,2,1)
figure; hold on;
plot(time(1:ntime-1001),mmean,'linewidth',4)
xlabel('t')
ylabel('Time Average of ||v||')
axis square

% subplot(1,2,2)
figure; hold on;
plot(time(1:ntime-1001),timeAveRMS,'linewidth',4)
xlabel('t')
ylabel('Time Average of ||v||_2')
axis square
% errorbar(time(2:30:end),mmean(1:30:end),sqrt(timeAveRMS(1:30:end)))
% semilogy(time(2:end),abs(mmean-mmean(end))/mmean(end),'linewidth',4)
