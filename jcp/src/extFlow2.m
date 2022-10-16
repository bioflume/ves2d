clear all;clc

nv = 1;   
n = 64;   
angle = 0;
b = 3; a = b/2; c=1;
t = (0:n-1)'*2*pi/n;
r = 0.5*sqrt( (a*cos(t)).^2 + (b*sin(t)).^2);
X0 = [c*r.*cos(t); r.*sin(t)];
% [RV0 A0 L0] = reducedVolume(X0)
% viewer(X0);
prams.nv = nv;
prams.kappa = 0;
prams.viscCont = 20;
prams.order = 1;
prams.vInf = @(X) farFieldVel(X,'extensional',-20);
prams.Case = 'fullSystem';

options.verbose = 0;
options.usePlot = 1;   
options.AxesHandle = gca;
options.axis = [-2 2 -2 2];      
options.progressBar = 1;
options.axisOn = 1;
options.saveData = 0;
 options.fileName = ['../results/extFlow2.mat'];

prams.T = 5;
[RV0 A0 L0] = reducedVolume(X0)
[X status] = Ves2D(X0,prams,options);
[RV A L] = reducedVolume(X)

M = [32 32 32 32  64  64  64  64  64  128];
N = [64 64 64 128 128 128 128 128 128 128];
for kk=1:10
  prams.T = .5;
  prams.m = M(kk);
  prams.n = N(kk);
  disp([kk prams.m prams.n]);
  clear functions global
  
  [X status] = Ves2D(X0,prams,options);
  
  X = interpft(reshape(X,[],2),prams.n); X = X(:);  
%   fileId = fopen(options.fileName,'a');
%   fwrite(fileId,X,'double');
%   fclose(fileId);
  X0 = X;
end
% 
% fileId = fopen(options.fileName,'r');
% Result = fread(fileId,'double');
% fclose(fileId);
% 
% head = 0;
% hold on;
% NN = 10;
% for ii =1:NN
%   tail = head+2*N(ii);
%   X = Result(head+1:tail); 
%   [RV A L] = reducedVolume(X);
%   disp([RV A L]);
%   n = size(X,1)/2;
%   plot(X([1:n 1]),X([n+1:2*n n+1]),'k','LineWidth',1.5);
%   head = tail;
% end
% hold off;
% axis equal



