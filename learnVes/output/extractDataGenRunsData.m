clear; clc;

iplot = 1;
fileName = 'dataGenRun.bin';
fid = fopen(fileName,'r');
val = fread(fid,'double');
fclose(fid);

% Some important numbers
dt = val(1);
nv = val(2);
N = val(3);
nvbdExt = val(4);
NbdExt = val(5);
nvbdInt = val(6);
NbdInt = val(7);

% Remove them and move on with rigid boundaries
val = val(8:end);
XwallsExt = val(1:2*NbdExt);
xxwallsExt = XwallsExt(1:end/2);
yywallsExt = XwallsExt(end/2+1:end);
val = val(2*NbdExt+1:end);
XwallsInt = val(1:2*nvbdInt*NbdInt);
XwallsInt = reshape(XwallsInt,[2*NbdInt nvbdInt]);
val = val(2*nvbdInt*NbdInt+1:end);

% Now read in time steps
X = []; sigma = []; time = [];
etaExt = []; etaInt = []; RS = [];
nEntry2Read = 1 + 3 * N * nv + 2 * NbdExt + 2 * nvbdInt * NbdInt + 3*nvbdInt;
ist = 1;
while(numel(val)>=nEntry2Read)
time(ist) = val(1);
Xst = reshape(val(2:2*N*nv+1),[2*N nv]);
X(:,:,ist) = Xst;

val = val(2*N*nv+2:end);
sigma(:,:,ist) = reshape(val(1:N*nv), [N nv]);

val = val(N*nv+1:end);
etaExt(:,ist) = val(1:2*NbdExt);
val = val(2*NbdExt+1:end);
etaInt(:,:,ist) = reshape(val(1:2*nvbdInt*NbdInt), [2*NbdInt nvbdInt]);
val = val(2*NbdInt*nvbdInt+1:end);
RS(:,:,ist) = reshape(val(1:3*nvbdInt+3), [3 nvbdInt+1]);
val = val(3*nvbdInt+4:end);

ist = ist + 1;
end

% plot data
if iplot
for it = 1 : numel(time)
figure(1); clf;
plot([XwallsExt(1:end/2); XwallsExt(1)], [XwallsExt(end/2+1:end); XwallsExt(end/2+1)],'k','linewidth',2)
hold on
axis equal
plot([XwallsInt(1:end/2,:); XwallsInt(1,:)], [XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k','linewidth',2)
hWalls = fill([XwallsInt(1:end/2,:); XwallsInt(1,:)],[XwallsInt(end/2+1:end,:); XwallsInt(end/2+1,:)],'k');
set(hWalls,'edgecolor','k')
Xst = X(:,:,it);
plot([Xst(1:end/2,:); Xst(1,:)], [Xst(end/2+1:end,:);Xst(end/2+1,:)], 'r', 'linewidth',2)
hVes = fill([Xst(1:end/2,:); Xst(1,:)],[Xst(end/2+1:end,:);Xst(end/2+1,:)],'r');
set(hVes,'edgecolor','r')
title(['Time = ' num2str(time(it))])
pause(0.1)
end


end