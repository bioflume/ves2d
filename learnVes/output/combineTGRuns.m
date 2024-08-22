fileNameTR = 'taylorGreen_IC3_trueFiner_diff625kNetJune8_dt5e-06_speed500.bin';
[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);

fileNameTR = 'taylorGreen_IC3_trueFiner_resume_diff625kNetJune8_dt5e-06_speed500.bin';
[vesxT2, vesyT2, ten, timeT2, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);
tLast = timeT(end);
for it = 2 : numel(timeT2)

  vesxT(:,:,end+1) = vesxT2(:,:,it);
  vesyT(:,:,end+1) = vesyT2(:,:,it);
  timeT(end+1) = tLast + timeT(it);
end

save taylorGreen_IC3_trueFiner_long_dt5E6_speed500 vesxT vesyT timeT


%%
% fileNameNN = 'taylorGreen_IC3_nearNet_diff625kNetJune8_dt1e-05_speed500.bin';
% [vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameNN);

load taylorGreen_IC3_nearNet_diff625kNetJune8_long_dt1E5_speed500

fileNameNN = 'taylorGreen_IC3_nearNet_resume2_diff625kNetJune8_dt1e-05_speed500.bin';
[vesxN2, vesyN2, ten, timeN2, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameNN);

tLast = timeN(end);
for it = 2 : numel(timeN2)

  vesxN(:,:,end+1) = vesxN2(:,:,it);
  vesyN(:,:,end+1) = vesyN2(:,:,it);
  timeN(end+1) = tLast + timeN(it);
end

save taylorGreen_IC3_nearNetLonger_diff625kNetJune8_long_dt1E5_speed500 vesxN vesyN timeN

%% 
fileNameTR = 'taylorGreen_IC3_true_diff625kNetJune8_dt1e-05_speed500.bin'; % BIEM with Near but resolution is low
[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);

fileNameTR = 'taylorGreen_IC3_true_resume_diff625kNetJune8_dt1e-05_speed500.bin'; % BIEM with Near but resolution is low
[vesxT2, vesyT2, ten, timeT2, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);
tLast = timeT(end);
for it = 2 : numel(timeT2)

  vesxT(:,:,end+1) = vesxT2(:,:,it);
  vesyT(:,:,end+1) = vesyT2(:,:,it);
  timeT(end+1) = tLast + timeT(it);
end

save taylorGreen_IC3_true_long_dt5E6_speed500 vesxT vesyT timeT
