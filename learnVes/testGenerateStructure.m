%function [XwallsInt,XwallsExt,area0,len0] = initializeWalls(prams,oc)
clear;
addpath ../src/
oc = curve;
prams.outWallRad = 2;
prams.inWallScale = 0.45;
prams.NbdExt = 256;
prams.NbdInt = 64;
prams.nvbdInt = 10;
prams.nvbdExt = 1;

% outer boundary
thet = (0:prams.NbdExt-1)'*2*pi/prams.NbdExt;
XwallsExt = [prams.outWallRad * cos(thet); prams.outWallRad*sin(thet)];
angle = pi * rand(prams.nvbdInt,1);
% inner boundary 
thet = (0:prams.NbdInt-1)'*2*pi/prams.NbdInt;
% sample shapes
% 1: circle, 2: curly, 3: cShape, 4: star3, 5: star4,
% 6: wigglyStar

ishapes = randi(9,[1,prams.nvbdInt]);
XwallsInt = zeros(prams.NbdInt*2,prams.nvbdInt);
Xps = zeros(prams.NbdInt*2,prams.nvbdInt);
sizes = zeros(prams.nvbdInt,1);
for ik = 1 : prams.nvbdInt
  if ishapes(ik) == 1 % circle
    Xp = [cos(-thet); sin(-thet)];
  elseif ishapes(ik) == 2 % optim-shape1
    points = [[6.6 3.9]; [15.1 6.3]; [5.6 6.2]; [-10.6 5.6]; [-2.0 4.1]; ...
        [-2.8 -6.5]; [6.0 -13.7]; [1.8 -10.8]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);
  elseif ishapes(ik) == 3 % cShape
    Xp = [-(1.5+sin(-thet)).*(cos(0.85*pi*cos(-thet))); (1.5+sin(-thet)).*(sin(0.85*pi*cos(-thet)))];
  elseif ishapes(ik) == 4 % star3
    folds = 3;
    radius = 1 + 0.3*cos(-folds*thet);
    Xp = [radius.*cos(-thet);radius.*sin(-thet)];
  elseif ishapes(ik) == 5 % star4
    folds = 4;
    radius = 1 + 0.3*cos(-folds*thet);
    Xp = [radius.*cos(-thet);radius.*sin(-thet)];  
  elseif ishapes(ik) == 6 % optim-shape2
    points = [[6.3 2.7]; [11.6 6.0]; [4.3 7.7]; [-12.8 5.1]; [-6.3 4.9]; ...
        [-4.7 -6.1]; [5.3 -13.6]; [1.0 -8.6]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);
  elseif ishapes(ik) == 7 % optim-shape3
    points = [[0.3 4.2]; [13.8 11.6]; [-2.0 10.6]; [-12.0 10.4]; [-6.9 -4.3]; ...
        [-6.3 -0.9]; [2.5 -9.2]; [1.0 -6.1]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);
  elseif ishapes(ik) == 8 % optim-shape4
    points = [[0.6 4.4]; [14.5 10.1]; [-1.8 10.2]; [-12.1 7.9]; [-6.9 1.4]; ...
        [-10.3 -6.9]; [3.6 -7.5]; [2.3 -9.9]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);
  elseif ishapes(ik) == 9 % optim-shape5
    points = [[10.0 9.1]; [0.0 9.1]; [-10.0 9.1]; [-9.1 7.5]; [-0.9 -7.5]; ...
        [0 -9.1]; [0.9 -7.5]; [9.1 7.5]];
    [Xp,Dx,Dy] = oc.samplePost(prams.NbdInt, points);  
  end
  % Scale down the arc-length
  Xp = [Xp(1:end/2)-mean(Xp(1:end/2)); Xp(end/2+1:end)-mean(Xp(end/2+1:end))];
  Lx = max(Xp(1:end/2))-min(Xp(1:end/2));
  Ly = max(Xp(end/2+1:end))-min(Xp(end/2+1:end));
  sizes(ik) = max(Lx,Ly);
  Xp = prams.inWallScale * Xp/sizes(ik);
  Xps(1:end/2,ik) = Xp(1:end/2)*cos(angle(ik)) - Xp(end/2+1:end)*sin(angle(ik));
  Xps(end/2+1:end,ik) = Xp(1:end/2)*sin(angle(ik)) + Xp(end/2+1:end)*cos(angle(ik));
  Lx = max(Xps(1:end/2,ik))-min(Xps(1:end/2,ik));
  Ly = max(Xps(end/2+1:end,ik))-min(Xps(end/2+1:end,ik));
  sizes(ik) = max(Lx,Ly);
end

% Place the interior walls
xcs = zeros(prams.nvbdInt,1); ycs = zeros(prams.nvbdInt,1); 
for ik = 1 : prams.nvbdInt
  K = [(1:ik-1) (ik+1:prams.nvbdInt)];
  iplace = false; iter = 0;
  while ~iplace
    xc = -(0.8*prams.outWallRad) + 1.6*prams.outWallRad*rand();
    yc = -(0.8*prams.outWallRad) + 1.6*prams.outWallRad*rand();
    dx = xc-xcs(K); dy = yc-ycs(K);
    dists = sqrt(dx.^2 + dy.^2)-(sizes(ik)/2+sizes(K)/2);

    % Check with the wall
    if sqrt(xc^2 + yc^2) + sizes(ik)/2 >= 0.9*prams.outWallRad
        iplace = false; 
    else
        iplace = true;
    end
    if iplace == true
        if any(dists < 0.2*prams.inWallScale)
            iplace = false; 
        end
    end
    iter = iter + 1;
    disp([num2str(ik) 'th post, iteration: ' num2str(iter)])
  end
  if iplace
    
    disp('Placed.') 
    xcs(ik) = xc; ycs(ik) = yc;
    XwallsInt(:,ik) = [Xps(1:end/2,ik)+xc; Xps(end/2+1:end,ik)+yc];
  end
end

clf
plot(XwallsExt(1:end/2), XwallsExt(end/2+1:end),'linewidth',2)
hold on
axis equal
plot(XwallsInt(1:end/2,:), XwallsInt(end/2+1:end,:),'linewidth',2)

