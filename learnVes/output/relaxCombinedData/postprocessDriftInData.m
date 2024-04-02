clear; clc;
mul = 2.^(0:10)';
dts = 1E-5*mul;

% load n128Dt1e-05RelaxDataSet.mat
Xnew1 = XnewStandStore;
Xold1 = XstandStore;
nInstances1 = nInstances;

load n128Dt1e-05RelaxMirrdDataSet.mat
% XnewStandStore = [Xnew1 XnewStandStore];
% XstandStore = [Xold1 XstandStore];
% nInstances = nInstances1 + nInstances;



if 0
dcxs = zeros(nInstances,1);
dcys = zeros(nInstances,1);

for ives = 1 : nInstances
  dcxs(ives) = mean(XnewStandStore(1:end/2,ives)) - mean(XstandStore(1:end/2,ives));
  dcys(ives) = mean(XnewStandStore(end/2+1:end,ives)) - mean(XstandStore(end/2+1:end,ives));
end

[fX,xiX] = ksdensity(dcxs,linspace(min(dcxs),max(dcxs),1000)');
[fY,xiY] = ksdensity(dcys,linspace(min(dcys),max(dcys),1000)');
fX = fX/trapz(xiX,fX);
fY = fY/trapz(xiY,fY);

disp(['Mean drift in X: ' num2str(mean(dcxs))])
disp(['Mean drift in Y: ' num2str(mean(dcys))])

subplot(1,2,1)
plot(xiX,fX,'linewidth',4)
axis square
xlabel('Drift in x-direction')
ylabel('PDF')
grid on

subplot(1,2,2)
plot(xiY,fY,'linewidth',4)
axis square
xlabel('Drift in y-direction')
ylabel('PDF')
grid on
end

