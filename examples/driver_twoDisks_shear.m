clear;

xcenters = [-8 0];
ycenters = [0.5 0];
shearStrength = 1;
dt = 0.4;
Ns = [16;32;64;128;256];

for in = 1 : numel(Ns)
[trajectories,velocities,ave_iter,max_iter] = twoSpheresShear_TestSymmAlpert(xcenters, ycenters, shearStrength, dt, Ns(in));
ave_iters(in,1) = ave_iter;
max_iters(in,1) = max_iter;
trajs{in} = trajectories;
vels{in} = velocities;
end


%%
for in = 1 : numel(Ns)
figure(1);clf;
plot(trajs{in}(1,:),trajs{in}(2,:),'r','linewidth',2)
hold on
plot(trajs{in}(3,:),trajs{in}(4,:),'b','linewidth',2)
axis equal
grid

xlim([-10 2])
ylim([-1 1])

title(['N = ' num2str(Ns(in))])
ax = gca;
exportgraphics(ax,['~/Desktop/N' num2str(Ns(in)) '.png'],'Resolution',300)

end


save LukasShearData
