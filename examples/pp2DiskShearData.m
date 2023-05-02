load LukasShearData1stKind.mat

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
