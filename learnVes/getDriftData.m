clear; clc;

load parabolic_n128Dt1e-05RelaxDataSet.mat

XinFFT = zeros(18,2,nInstances);
DriftXY = zeros(2,nInstances);
modes = [(1:N/2-1) (-N/2:-1)];

for k = 1 : nInstances
z = XstandStore(1:end/2,k) + 1i*XstandStore(end/2+1:end,k);
zh = fft(z)/N;

XinFFT(:,1,k) = real(zh(abs(modes)<10));
XinFFT(:,2,k) = imag(zh(abs(modes)<10));
DriftXY(1,k) = mean(XnewStandStore(1:end/2,k)) - mean(XstandStore(1:end/2,k));
DriftXY(2,k) = mean(XnewStandStore(end/2+1:end,k)) - mean(XstandStore(end/2+1:end,k));
end

save parabolic_Dt1e-05Relax_DriftData XinFFT DriftXY nInstances 
