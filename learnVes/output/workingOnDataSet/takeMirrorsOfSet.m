clear; clc;
addpath ../../src/
oc = curve;
% load n128Dt0.0001RelaxAllDataSet
% load n128Dt1e-05RelaxAllDataSet
load n128Dt5e-05RelaxAllDataSet

XstandStoreMirrY = zeros(2*N,nInstances);
XstandStoreMirrYX = zeros(2*N,nInstances);
XstandStoreMirrX = zeros(2*N,nInstances);

XnewStandStoreMirrY = zeros(2*N,nInstances);
XnewStandStoreMirrYX = zeros(2*N,nInstances);
XnewStandStoreMirrX = zeros(2*N,nInstances);
N = 128;

for k = 1 : nInstances
  
  Xref = [-XstandStore(1:end/2,k);XstandStore(end/2+1:end,k)];
  [XstandStoreMirrY(:,k),sortIdx,translation] = sortPoints(Xref,oc);
  Xref = [-XnewStandStore(1:end/2,k);XnewStandStore(end/2+1:end,k)];
  Xref = [Xref(1:end/2)+translation(1);Xref(end/2+1:end)+translation(2)];
  XnewStandStoreMirrY(:,k) = [Xref(sortIdx);Xref(sortIdx+N)];

  Xref = [XstandStoreMirrY(1:end/2,k);-XstandStoreMirrY(end/2+1:end,k)];
  [XstandStoreMirrYX(:,k),sortIdx,translation] = sortPoints(Xref,oc);
  Xref = [XnewStandStoreMirrY(1:end/2,k);-XnewStandStoreMirrY(end/2+1:end,k)];
  Xref = [Xref(1:end/2)+translation(1);Xref(end/2+1:end)+translation(2)];
  XnewStandStoreMirrYX(:,k) = [Xref(sortIdx);Xref(sortIdx+N)];

  Xref = [-XstandStoreMirrYX(1:end/2,k);XstandStoreMirrYX(end/2+1:end,k)];
  [XstandStoreMirrX(:,k),sortIdx,translation] = sortPoints(Xref,oc);
  Xref = [-XnewStandStoreMirrYX(1:end/2,k);XnewStandStoreMirrYX(end/2+1:end,k)];
  Xref = [Xref(1:end/2)+translation(1);Xref(end/2+1:end)+translation(2)];
  XnewStandStoreMirrX(:,k) = [Xref(sortIdx);Xref(sortIdx+N)];

  % 
  % figure(1); clf;
  % plot(XstandStore(1:end/2,k),XstandStore(end/2+1:end,k),'k','linewidth',2)
  % hold on
  % plot(XnewStandStore(1:end/2,k),XnewStandStore(end/2+1:end,k),'r','linewidth',2)
  % axis equal
  % 
  % figure(2);clf;
  % plot(XstandStoreMirrY(1:end/2,k),XstandStoreMirrY(end/2+1:end,k),'k','linewidth',2)
  % hold on
  % plot(XnewStandStoreMirrY(1:end/2,k),XnewStandStoreMirrY(end/2+1:end,k),'r','linewidth',2)
  % axis equal
  % 
  % figure(3);clf;
  % plot(XstandStoreMirrYX(1:end/2,k),XstandStoreMirrYX(end/2+1:end,k),'k','linewidth',2)
  % hold on
  % plot(XnewStandStoreMirrYX(1:end/2,k),XnewStandStoreMirrYX(end/2+1:end,k),'r','linewidth',2)
  % axis equal
  % 
  % figure(4);clf;
  % plot(XstandStoreMirrX(1:end/2,k),XstandStoreMirrX(end/2+1:end,k),'k','linewidth',2)
  % hold on
  % plot(XnewStandStoreMirrX(1:end/2,k),XnewStandStoreMirrX(end/2+1:end,k),'r','linewidth',2)
  % axis equal
  % title(k)
  % pause
  
end

XstandStore = [XstandStore XstandStoreMirrY XstandStoreMirrYX XstandStoreMirrX];
XnewStandStore = [XnewStandStore XnewStandStoreMirrY XnewStandStoreMirrYX XnewStandStoreMirrX];

nInstances = nInstances*4;

save newStand625K_n128Dt5e-05RelaxDataSet_June8 XstandStore XnewStandStore nInstances dt kappa N


function [Xsort,sortIdx,translation] = sortPoints(Xref,oc)
N = 128;
center = oc.getPhysicalCenterShan(Xref);
translation = -center;
Xref = [Xref(1:end/2)-center(1);Xref(end/2+1:end)-center(2)];

firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];
Xsort = [Xref(sortIdx);Xref(sortIdx+N)];
end