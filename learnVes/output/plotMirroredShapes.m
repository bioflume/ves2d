load n128Dt1e-05RelaxMirrdDataSet.mat
for k = 1 : nInstances/2
figure(1);clf;
id1 = 2*k-1;
id2 = 2*k;

vecx1 = XnewStandStore(1:end/2,id1);
vecy1 = XnewStandStore(end/2+1:end,id1);

vecx2 = XnewStandStore(1:end/2,id2);
vecy2 = -XnewStandStore(end/2+1:end,id2);

plot([vecx1;vecx1(1)], [vecy1;vecy1(1)],'linewidth',2)
hold on
plot([vecx2;vecx2(1)], [vecy2;vecy2(1)],'linewidth',2)
axis equal
% pause;
% ax = gca;
% exportgraphics(ax,['~/Desktop/ves' num2str(k) '.png'],'Resolution',300)


% figure(2);clf;
% 
% mirrVecx1 = -vecx1;
% mirrVecx2 = -vecx2;
% plot([mirrVecx1;mirrVecx1(1)],[vecy1;vecy1(1)],'linewidth',2)
% hold on
% plot([mirrVecx2;mirrVecx2(1)],[vecy2;vecy2(1)],'linewidth',2)
% axis equal
% ax = gca;
% exportgraphics(ax,['~/Desktop/xMirrves' num2str(k) '.png'],'Resolution',300)
pause


end


%% 
for k = 1 : nInstances
cx(k,1) = mean(XnewStandStore(1:end/2,k))-mean(XstandStore(1:end/2,k));
cy(k,1) = mean(XnewStandStore(end/2+1:end,k))-mean(XstandStore(end/2+1:end,k));
end