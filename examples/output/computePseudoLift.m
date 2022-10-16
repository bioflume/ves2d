clear; 
addpath ../../src/
oc = curve;

% Mat-file storing stress and other data
fileName = 'N6K3VC3_Stress';
load(fileName)

% number of target points
N = size(posx,1);

% if DLD, then get the center of the 1st pillar
if iDiffDiscWalls
  pillarCx = mean(intWallx(:,2));
end

% Compute the total stress
if iDiffDiscWalls
  stress1 = stress1Store+stress1ExtStore+stress1SelfStore;
  stress2 = stress2Store+stress2ExtStore+stress2SelfStore;
else
  stress1 = stress1Store+stress1SelfStore;
  stress2 = stress2Store+stress2SelfStore;  
end

% ADD STRESSES
stress11 = stress1(1:N,:);
stress12 = stress1(N+1:end,:);
stress21 = stress2(1:N,:);
stress22 = stress2(N+1:end,:);

% compute the velocity on the vesicle boundary
velx = zeros(N,numel(time)); vely = velx;
centx = zeros(numel(time),1); centy = centx;
for tt = 2:numel(time)
  velx(:,tt) = (posx(:,1,tt)-posx(:,1,tt-1))/(time(tt)-time(tt-1)); 
  vely(:,tt) = (posy(:,1,tt)-posy(:,1,tt-1))/(time(tt)-time(tt-1)); 
  
  centx(tt,1) = mean(interpft(posx(:,1,tt),256));
  centy(tt,1) = mean(interpft(posy(:,1,tt),256));
end

centx(1) = mean(interpft(posx(:,1,1),256));
centy(1) = mean(interpft(posy(:,1,1),256));

% also compute the velocity of the center
centVelX = [0; (centx(2:end)-centx(1:end-1))./(time(2:end)-time(1:end-1))];
centVelY = [0; (centy(2:end)-centy(1:end-1))./(time(2:end)-time(1:end-1))];

% GET UNIT VECTORS IN VELOCITY DIRECTION AND PERP. TO VEL. DIRECTION
% IN VEL. DIR.--> DRAG, IN PERP. TO VEL. DIR.--> LIFT
% unit vector in velocity direction
normvelx = velx(:,2:end)./(sqrt(velx(:,2:end).^2+vely(:,2:end).^2));
normvely = vely(:,2:end)./(sqrt(velx(:,2:end).^2+vely(:,2:end).^2));
% unit vector perpendicular to the velocity
perpnormvelx = -normvely; perpnormvely = normvelx;
% since velocity at time t = 0 is zero.
perpnormvelx = [zeros(N,1) perpnormvelx];
perpnormvely = [zeros(N,1) perpnormvely];
normvelx = [zeros(N,1) normvelx];
normvely = [zeros(N,1) normvely];

% COMPUTE LIFT AT EACH POINT, THEN INTEGRATE OVER SURFACES
timeSkipped = time(1:skip:end); 
LSurfAve = zeros(numel(timeSkipped),1);

for tt = 1 : numel(timeSkipped)
  % Get the directions at time step tt
  pnormx = perpnormvelx(:,(tt-1)*skip+1);
  pnormy = perpnormvely(:,(tt-1)*skip+1);
  vnormx = normvelx(:,(tt-1)*skip+1);
  vnormy = normvely(:,(tt-1)*skip+1);
 
  % Get the configuration at time step tt
  X = [posx(:,1,(tt-1)*skip+1);posy(:,1,(tt-1)*skip+1)];
  % jacobian, tangent
  [jac,tang,~] = oc.diffProp(X); 
  tanx = tang(1:end/2 ); tany = tang(end/2+1:end);
  % normal
  nx = -tany; ny = tanx;
    
  LiftAtPnt = zeros(N,1); 
  Tn1AtPnt = [];
  Tn2AtPnt = [];
  for idx = 1 : N 
    % Total stress tensor at time tt and point idx multiplied by normal
    Tn = [stress11(idx,tt) stress12(idx,tt);stress21(idx,tt) stress22(idx,tt)]*...
        [nx(idx);ny(idx)]; 
    % Lift at a point (Tn dotted with unit vector perpendicular to velocity
    % vector
    LiftAtPnt(idx,1) = sum(Tn.*[pnormx(idx);pnormy(idx)]);
  end

  % Take the surface average of the lift
  LSurfAve(tt,1) = 2*pi/N*sum(LiftAtPnt.*jac);
end

% Decide when to start averaging 
centxSkipped = centx(1:skip:end);
centySkipped = centy(1:skip:end);
if iDiffDiscWalls % if DLD
    
  % Start averaging from the region between the first two gaps
  whenStart = find(centxSkipped>=pillarCx+1.25); 
  
  % -------------------------------------------------------------------------
  % Start averaging above the first pillar
  % whenStart = find(centxSkipped>=pillarCx); 
  % -------------------------------------------------------------------------
  
  whenStart = whenStart(1);
  
  % Find the centers of pillars
  % Get the centers of each pillar
  pillCentx = zeros(size(intWallx,2),1); pillCenty = pillCentx;
  for k = 1 : size(intWallx,2)
    pillCentx(k) = mean(interpft(intWallx(:,k),256));
    pillCenty(k) = mean(interpft(intWally(:,k),256));  
  end
  
  % Find distance to the pillars
  closeDistx = zeros(numel(timeSkipped,1)); closeDisty = closeDistx;
  for tt = 1 : numel(timeSkipped)
    dist2pills = sqrt((centxSkipped(tt)-pillCentx).^2+(centySkipped(tt)-...
        pillCenty).^2);
    [closeDist,id] = min(dist2pills);
    closeDistx(tt) = centxSkipped(tt)-pillCentx(id);
    closeDisty(tt) = centySkipped(tt)-pillCenty(id);
  end
  % Average until
  whenEnd = find(closeDisty<-0.1); % when zig-zagged
  
  if isempty(whenEnd) % then displaced so average to the end
    whenEnd = numel(centxSkipped);
  else
    disp('It zig-zags')  
    whenEnd = whenEnd(1);  
  end

% -------------------------------------------------------------------------
%   whenEndZZ = find(closeDisty<-0.1); % when zig-zagged
%   whenEnd = find(centxSkipped>=pillarCx+6.25); % until the gap between 3rd and 4th pillars
%   whenEnd = find(centxSkipped>=pillarCx+5); % until above the third pillar
  
%   if ~isempty(whenEndZZ)
%     whenEnd = min(whenEnd(1),whenEndZZ(1)); % if zzs earlier, then average until that.
%   else
%     whenEnd = whenEnd(1);
%   end
% -------------------------------------------------------------------------
else
  % if a channel flow
  whenStart = find(centySkipped>=centySkipped(1)*0.8); % start from begining
  if ~isempty(whenStart)
    whenStart = whenStart(1);
  else
    whenStart = 10;
  end
  whenEnd = find(centySkipped<=-0.1); % when it reaches center stop averaging
  if ~isempty(whenEnd)
    whenEnd = whenEnd(end);
  else
    whenEnd = numel(centySkipped);  
  end
end

% Time averaged pseudo-lift
LtimeAve = 1/kappa*1/(timeSkipped(whenEnd)-timeSkipped(whenStart))*...
    trapz(timeSkipped(whenStart:whenEnd),LSurfAve(whenStart:whenEnd));
% LtimeAve = 1/kappa*1/numel(LSurfAve(whenStart:whenEnd))*...
%     sum(LSurfAve(whenStart:whenEnd));
fprintf('Time averaged lift is %2.2e.\n',LtimeAve);

% Time averaged migration velocity
migVel = centVelY(1:skip:end);
VyTimeAve = 1/(timeSkipped(whenEnd)-timeSkipped(whenStart))*...
    trapz(timeSkipped(whenStart:whenEnd),migVel(whenStart:whenEnd));
% VyTimeAve = 1/numel(migVel(whenStart:whenEnd))*...
%     sum(migVel(whenStart:whenEnd));
fprintf('Time averaged migration velocity is %2.2e.\n',VyTimeAve);

figure(1);clf;hold on;
if iDiffDiscWalls
  plot(extWallx,extWally,'k','linewidth',3);
  plot(intWallx,intWally,'k','linewidth',2);
else
  plot(wallx,wally,'k','linewidth',3)
end
plot(centx,centy,'b','linewidth',2)
plot(centxSkipped(whenStart:whenEnd),centySkipped(whenStart:whenEnd),'r')
axis equal;
figure(2);clf;
plot(centxSkipped,LSurfAve,'b','linewidth',2)
hold on;
plot(centxSkipped(whenStart:whenEnd),LSurfAve(whenStart:whenEnd),'r','linewidth',2)
xlabel('time')
ylabel('Pseudo-lift')
legend('Entire simulation','Duration for averaging')

figure(3);clf;
plot(centxSkipped,migVel,'b','linewidth',2)
hold on;
plot(centxSkipped(whenStart:whenEnd),migVel(whenStart:whenEnd),'r','linewidth',2)
xlabel('time')
ylabel('Migration velocity')
legend('Entire simulation','Duration for averaging')

figure(2);
