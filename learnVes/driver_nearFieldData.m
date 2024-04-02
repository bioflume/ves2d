
% This is for the self interaction part
%npar = 4;
%parfor i = 1 : npar
%prepareNearFieldVelocityData(i,npar);
%end


% Below is for the action of near-field operators on Vinf
npar = 4;

parfor i = 1 : npar
prepareNearFieldVinfFFTData(i,npar);
end