% clear all;clc
% for n = 2.^(6:14)
%   disp(['-- n = ' num2str(n) ' ------------------------']);
%   x = 20*rand(1,n); y = 20*rand(1,n);
%   fx = rand(1,n); fy = rand(1,n);
% 
%   tic; F = exact_stokes([x;y],[fx;fy],0); t1=toc;
%   F = F*4*pi;
%   tic; [px py] = fastsummation(x,y,fx,fy); t2=toc;
%   P = [px;py];
%   e = F-P; disp(['error: ' num2str(sqrt(max(dot(e,e)./dot(F,F))))]);
%   disp(['Time: ' num2str(t1/t2)]);
% end

clear all;clc

ns = 256;
for nt = 2.^(6:14)
  %clear functions
  disp(['-- n =' num2str(nt) ' ------------------------']);
  
  bd = @(ind,m) sampleBd(ind,m,1,'couette');      
  bound = @(n) fixedBound(n,bd,1);
  
  [xx yy] = pol2cart(2*pi*rand(nt,1),26+4*(rand(nt,1)-.5));
  
  mu = 10*rand(4*ns+3,1);

  tic; you1 = evalField([ns ns],mu,bound,xx,yy,[]); t1=toc;
  tic; you2 = evalFieldDirect([ns ns],mu,bound,xx,yy,[]); t2=toc;
  
  you1 = reshape(you1,[],2);    you2 = reshape(you2,[],2);  
  e = you1-you2; disp(['error: ' num2str(sqrt(max(dot(e,e)./dot(you1,you1))))]);
  disp(['Speedup: ' num2str(t1/t2)]);
end



