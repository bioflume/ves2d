clear;

load n256Dt1E4Kb1E1nModes64PCAData.mat
XstandAll = XoldCoeffs(:,1:64);
XnewStandAll = XnewCoeffs(:,1:64);
kappas = 1e-1*ones(nInstances,1);
Dts = 1e-4*ones(nInstances,1);
DtTimesKappa = 1e-5*ones(nInstances,1);
nSamples = nInstances;
nInstances

load n256Dt1E4Kb1E2nModes64PCAData.mat
XstandAll = [XstandAll;XoldCoeffs(:,1:64)];
XnewStandAll = [XnewStandAll;XnewCoeffs(:,1:64)];
kappas = [kappas;ones(nInstances,1)*1e-2];
Dts = [Dts;ones(nInstances,1)*1e-4];
DtTimesKappa = [DtTimesKappa;ones(nInstances,1)*1e-6];
nSamples = nSamples + nInstances;
nInstances

load n256Dt5E4Kb1E1nModes64PCAData.mat
XstandAll = [XstandAll;XoldCoeffs(:,1:64)];
XnewStandAll = [XnewStandAll;XnewCoeffs(:,1:64)];
kappas = [kappas;ones(nInstances,1)*1e-1];
Dts = [Dts;ones(nInstances,1)*5e-4];
DtTimesKappa = [DtTimesKappa;ones(nInstances,1)*5e-5];
nSamples = nSamples + nInstances;
nInstances

load n256Dt5E4nModes64PCAData.mat
XstandAll = [XstandAll;XoldCoeffs(:,1:64)];
XnewStandAll = [XnewStandAll;XnewCoeffs(:,1:64)];
kappas = [kappas;ones(nInstances,1)*1];
Dts = [Dts;ones(nInstances,1)*5e-4];
DtTimesKappa = [DtTimesKappa;ones(nInstances,1)*5e-4];
nSamples = nSamples + nInstances;
nInstances

%load n256Dt5E5Kb1E1nModes64PCAData.mat
%XstandAll = [XstandAll;XoldCoeffs(:,1:64)];
%XnewStandAll = [XnewStandAll;XnewCoeffs(:,1:64)];
%kappas = [kappas;ones(nInstances,1)*1e-1];
%Dts = [Dts;ones(nInstances,1)*5e-5];
%nSamples = nSamples + nInstances;
%nInstances

%load n256Dt5E5nModes64PCAData.mat
%XstandAll = [XstandAll;XoldCoeffs(:,1:64)];
%XnewStandAll = [XnewStandAll;XnewCoeffs(:,1:64)];
%kappas = [kappas;ones(nInstances,1)*1];
%Dts = [Dts;ones(nInstances,1)*5e-5];
%DtTimesKappa = [DtTimesKappa;ones(nInstances,1)*5e-5];
%nSamples = nSamples + nInstances;
%nInstances

load n256Dt5E4Kb1E2nModes64PCAData.mat
XstandAll = [XstandAll;XoldCoeffs(:,1:64)];
XnewStandAll = [XnewStandAll;XnewCoeffs(:,1:64)];
kappas = [kappas;ones(nInstances,1)*1e-2];
Dts = [Dts;ones(nInstances,1)*5e-4];
DtTimesKappa = [DtTimesKappa;ones(nInstances,1)*5e-6];
nSamples = nSamples + nInstances;
nInstances


load n256Dt5E5Kb1E2nModes64PCAData.mat
XstandAll = [XstandAll;XoldCoeffs(:,1:64)];
XnewStandAll = [XnewStandAll;XnewCoeffs(:,1:64)];
kappas = [kappas;ones(nInstances,1)*1e-2];
Dts = [Dts;ones(nInstances,1)*5e-5];
DtTimesKappa = [DtTimesKappa;ones(nInstances,1)*5e-7];
nSamples = nSamples + nInstances;
nInstances

load only1timeStepMoreModesPCAData.mat
XstandAll = [XstandAll;XoldCoeffs(:,1:64)];
XnewStandAll = [XnewStandAll;XnewCoeffs(:,1:64)];
kappas = [kappas;ones(nInstances,1)*1];
Dts = [Dts;ones(nInstances,1)*1e-4];
DtTimesKappa = [DtTimesKappa;ones(nInstances,1)*1e-4];
nSamples = nSamples + nInstances;
nInstances

nInstances = nSamples;
XoldCoeffs = XstandAll;
XnewCoeffs = XnewStandAll;
N = 256;

save('n256AllDtTimesKappasPCAData.mat','nInstances','XoldCoeffs',...
  'XnewCoeffs','N','Dts','kappas','DtTimesKappa','evects','colMeans','-v7.3')
