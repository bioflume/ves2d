clear;
% 
xcenters = [-8 0];
ycenters = [0.25 0.25];
shearStrength = 1;

% xcenters = [-4 4];
% ycenters = [0 0];
% force = [-1;0;0];
dt = 0.4;
Ns = [16;32;64;128;256];

for in = 1 : numel(Ns)
[trajectories,velocities,ave_iter,max_iter] = twoSpheresShear_TestSymmAlpert(xcenters, ycenters, shearStrength, dt, Ns(in));
% [trajectories,velocities,ave_iter,max_iter] = twoSpheresShear_Test1stKind(xcenters, ycenters, shearStrength, dt, Ns(in));
% [trajectories,velocities,ave_iter,max_iter] = twoSpheresForce_TestSymmAlpert(xcenters, ycenters, force, dt, Ns(in));
% [trajectories,velocities,ave_iter,max_iter] = twoSpheresForce_Test1stKind(xcenters, ycenters, force, dt, Ns(in));
ave_iters(in,1) = ave_iter;
max_iters(in,1) = max_iter;
trajs{in} = trajectories;
vels{in} = velocities;
end


%%
theta = [0:63]'/64 * 2 * pi;
vecx = cos(theta)*0.5;
vecy = sin(theta)*0.5;
vecx = [vecx;vecx(1)];
vecy = [vecy;vecy(1)];
minDists = [];
if 1
for in = 1 : numel(Ns)
figure(1); clf;
plot(trajs{in}(1,:),trajs{in}(2,:),'r','linewidth',2)
hold on
plot(trajs{in}(3,:),trajs{in}(4,:),'b','linewidth',2)

dx = (trajs{in}(1,:)-trajs{in}(3,:));
dy = (trajs{in}(2,:)-trajs{in}(4,:));
mindist = sqrt(dx.^2 + dy.^2)-1;
[val,id] = min(mindist);
minDists = [minDists;val];
diskx1 = vecx + trajs{in}(1,id);
disky1 = vecy + trajs{in}(2,id);

diskx2 = vecx + trajs{in}(3,id);
disky2 = vecy + trajs{in}(4,id);
    
plot(diskx1,disky1,'r','linewidth',2)
plot(diskx2,disky2,'b','linewidth',2)

axis equal
xlim([-4 4])
ylim([-1.5 1.5])

title(['N = ' num2str(Ns(in))])
ax = gca;
exportgraphics(ax,['~/Desktop/N' num2str(Ns(in)) '.png'],'Resolution',300)
end
end

figure(1); clf;
for in = 1 : numel(Ns)
plot(trajs{in}(1,:),trajs{in}(2,:),'linewidth',2)
hold on
end
axis equal
xlim([-4 4])
ylim([0 1.5])
legend('N = 16', 'N = 32','N = 64','N = 128','N = 256')
legend boxoff