disp('============================================')
clear all

%dbstop in evalTraction at 71
% TRAC = @(x,y) [2*x 1-2*y;1-2*y 2*x];
% PR = @(x,y) -2*x;

prams.M = 4*[32 32]; np = sum(prams.M);
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette');
flow = 'mixed';
prams.bc = @(x) forcing(x,flow);
domain = fixedBound(prams.M,prams.bd,1);
[x y] = meshgrid(-15:5:15);
IN = inDomain(domain,x(:),y(:),.1);
X = [x(IN) y(IN)];
[U trash mu] = farFieldVel(X,[],'confined',prams);
U = reshape(U,[],2);
XX = []; n = [];
for ind = 1:length(domain)
    XX = [XX;domain(ind).X];
    n = [n;domain(ind).n];
end
% Uref = forcing(X,flow);
% viewer(domain);
% hold on; quiver(X(:,1),X(:,2),U(:,1),U(:,2));
% hold on; quiver(X(:,1),X(:,2),Uref(:,1),Uref(:,2),'r');

% err = U-Uref;
% err = myNorm(err)./max(myNorm(Uref));
% fprintf('Velocity error: %8.7e\n',max(err));

uBound = forcing(XX,flow);
st1 = surfaceTraction(prams.M,domain,uBound(:));st1 = reshape(st1,[],2);
% for ind = 1:size(XX,1)
%     st1(ind,:) = TRAC(XX(ind,1),XX(ind,2))*n(ind,:)';
%     p1(ind,1) = PR(XX(ind,1),XX(ind,2));
% end
p1 = 0;
%avgS1 = avgStress(prams.M,domain,uBound(:),[],[]);

[st2 p2] = evalTraction(prams.M,domain,mu);
%avgS2 = avgStress(prams.M,domain,mu,[],[]);

%figure; plot(myNorm(st1)); hold on; plot(myNorm(st2),'r')
%figure; %plot(st1);
%hold on; plot(st2,'r')
err = st1-st2; err = myNorm(err)/max(myNorm(st1));
fprintf('%8.7e \t %8.7e\n',max(abs(p1-p2)),max(err));
figure;
subplot(1,2,1); plot(1:np,st1(:,1),'k',1:np,st2(:,1),'b');
subplot(1,2,2); plot(1:np,st1(:,2),'k',1:np,st2(:,2),'b');

%%=============================================================
%
% Ro = 15; Ri = 5;
% c = Ro/Ri; c = c/(c^2-1);
% vTheta = @(r) c*(Ro*dot(r,r,2).^(-1/2)-sqrt(dot(r,r,2))/Ro);
%
% vt = vTheta(X);
% Uref = zeros(size(U));
% Uref(IN,:) = [vt vt].*[-X(:,2)./sqrt(dot(X,X,2)) X(:,1)./sqrt(dot(X,X,2))];
%
% viewer(domain);
% hold on; quiver(x(:),y(:),U(:,1),U(:,2));
% hold on; quiver(x(:),y(:),Uref(:,1),Uref(:,2),'r');

% err = U-Uref;
% err = sqrt(dot(err,err,2)./dot(Uref,Uref,2));
% figure;surf(x,y,reshape(err,size(x))); view(2); colorbar;

% [st p] = evalTraction(prams.M,domain,mu);
% stSize = sqrt(dot(st,st,2))
%hold on; quiver(XX(:,1),XX(:,2),st(:,1),st(:,2));
%figure; plot(p);
%avgS = avgStress(prams.M,domain,mu,[],[]);

% options.usePlot = 1;
% options.AxesHandle = gca;
%
% X = [];
% Results = Ves2D(X,prams,options,@monitor);
% syms r
% Ri = 5; Ro = 15; c = Ro/Ri;
% vt =c*(Ro/r-r/Ro)/(c^2-1);
% dv = simple(diff(vt,'r'));
% subs(dv,'r',15);
%
%
%