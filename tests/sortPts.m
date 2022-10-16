% test routine that will sort points according to being inside or
% outside of a vesicle

addpath ../src

N = 64; theta = (0:N-1)'*2*pi/N;

X = [[cos(theta);3*sin(theta)]...
     [3*cos(theta)+5;sin(theta)-2]...
     [3*cos(theta)+4;sin(theta)+2]];
% some arbitrary collection of ellipses
vesicle = capsules(X,[],[],1,1,0);


ntra = 200;
[xtar,ytar] = meshgrid(linspace(-2,10,ntra),linspace(-4,4,ntra));
xtar = xtar(:);ytar = ytar(:);
Xtar = [xtar;ytar];
% set of target points that we want to bin as being outside or inside

fmm = true;
fprintf('USING FMM\n')
tic
InOut1 = vesicle.sortPts(Xtar,fmm);
toc

fmm = false;
fprintf('USING DIRECT\n')
tic
InOut2 = vesicle.sortPts(Xtar,fmm);
toc

fprintf('Error between two results is %4.2e\n',norm(InOut1 - InOut2))

if 0
InOut = InOut1;
figure(1); clf; hold on
plot(X(1:end/2,:),X(end/2+1:end,:),'r')
plot(xtar,ytar,'b.')
plot(xtar(InOut==0),ytar(InOut==0),'g.')
end



