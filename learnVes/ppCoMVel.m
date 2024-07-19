clear; clc;
addpath ../src
addpath ./output/


R = 0.1291;
chanWidths = R ./ [0.2; 0.4; 0.6; 0.75];
speeds = [1500 3000 4500 6000 7500;
    750 1500 2250 3000 3750;
    500 1000 1500 2000 2500;
    400 800 1200 1600 2000];
Cks = zeros(numel(chanWidths)*3,1);
Cns = Cks;
count = 1;
for iw = 1 : numel(chanWidths)
    for is = 1 : 5
      w = chanWidths(iw);
      vmax = speeds(iw,is);
      Cks(count,1) = 2*vmax*R^3/w;
      Cns(count,1) = R/w;
      count = count + 1;
    end
end

N = 128;
oc = curve;
op = poten(N);

for iw = 1 : numel(chanWidths)
    for is = 1 : 5
    speed = speeds(iw,is);
    chanWidth = chanWidths(iw);

    runNew = ['./output/truePoisRuns/poisTrueRuns_speed' num2str(speed) '_width' num2str(chanWidth) '.bin'];
    [vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runNew);
    vesxN = vesxN(:,end-150:end-50);
    vesyN = vesyN(:,end-150:end-50);
    timeN = timeN(end-150:end-50);
    
    Xs = [vesxN;vesyN];
    Xs1 = Xs;
    Xs = filterShape(Xs);
    % for k = 1 : nv
    %   figure(4);clf;
    %   plot(Xs1(1:end/2,k),Xs1(end/2+1:end,k),'k','linewidth',2)
    %   hold on
    %   plot(Xs(1:end/2,k),Xs(end/2+1:end,k),'r','linewidth',2)
    %   axis equal
    %   pause
    % end
    nv = numel(Xs(1,:));

    vesicle = capsules(Xs,[],[],1,1,0);
    [Ben, Ten, Div] = vesicle.computeDerivs;
    G = op.stokesSLmatrix(vesicle);

    tension = zeros(N,nv);
    for k = 1 : nv
      LHS = (Div(:,:,k)*G(:,:,k)*Ten(:,:,k));
      RHS = -Div(:,:,k)*(G(:,:,k)*(Ben(:,:,k)*Xs(:,k)));

      tension(:,k) = LHS\RHS;
    end

    tracJump = vesicle.tracJump(Xs,tension);
    vinf = @(X) [speed*(1-(X(end/2+1:end,:)/chanWidth).^2);...
      zeros(size(X(1:end/2,:)))];
    selfVel = zeros(2*N,nv);
    magVel = zeros(N,nv);
    for k = 1 : nv
      selfVel(:,k) = G(:,:,k)*tracJump(:,k);
      vback = vinf(Xs(:,k));
      selfVel = selfVel + vback;
      magVel(:,k) = sqrt(selfVel(1:end/2,k).^2 + selfVel(end/2+1:end,k).^2);
    end

    ppDataFile = ['./VelocityData/ppTrueData_Speed' num2str(speed) '_width' num2str(chanWidth) '.mat'];
    save(ppDataFile, 'Xs', 'selfVel')

    cxV = abs(movmean(mean(selfVel(1:end/2,:),1),1));
    cyV = abs(movmean(mean(selfVel(end/2+1:end,:),1),1));
    % mV = movmean(mean(magVel,1),500);
    mV = sqrt(cxV.^2 + cyV.^2);
    
    % figure(1);clf;
    % yyaxis left
    % plot(timeN,cxV./mV,'linewidth',2)
    % ylabel('Vx/V')
    % ylim([0.95 1])
    % axis square
    % yyaxis right
    % plot(timeN,cyV./mV,'linewidth',2)
    % ylim([0  0.05])
    % ylabel('Vy/V')
    % axis square
    % xlabel('Time')
    % grid
    % 
    % ax = gca;
    % figName = ['AveVel_speed' num2str(speed) '_width' num2str(chanWidth) '.png'];
    % exportgraphics(ax,figName,'Resolution',300)

    % figure(2); clf;
    % plot([Xs(1:end/2,end);Xs(1,end)]-mean(Xs(1:end/2,end)),[Xs(end/2+1:end,end);Xs(end/2+1,end)]-mean(Xs(end/2+1:end)),'r','linewidth',2)
    % xT = [Xs(1:end/2,end);Xs(1,end)]-mean(Xs(1:end/2,end));
    % yT = [Xs(end/2+1:end,end);Xs(end/2+1,end)]-mean(Xs(end/2+1:end));
    % hFill = fill(xT, yT,'r');
    % set(hFill,'EdgeColor','r')
    % axis equal
    % box on

    % pause

    
    end
end


function X = filterShape(X)
N = size(X,1)/2;
nv = size(X,2);


Nup = 1024;
modes = [(0:Nup/2-1) (-Nup/2:-1)];

for k = 1 : nv
    x = interpft(X(1:end/2,k),Nup);
    y = interpft(X(end/2+1:end,k),Nup);
    z = x + 1i*y;
    z = fft(z);
    z(abs(modes) > 16) = 0;
    z = ifft(z);
    X(:,k) = [interpft(real(z),N);interpft(imag(z),N)];
end

end