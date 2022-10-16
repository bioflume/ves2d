clear; clc;
ishape = 'super';
Nbdi = 256;
Nbdo = 512;
nvbd = 2;
WsizeX = 6; WsizeY = 10;
ai = 10; bi = 4; n = 4;
ao = ai + WsizeX; bo = bi + WsizeY;

ti = (0:Nbdi-1)'/Nbdi * 2 * pi;
to = (0:Nbdo-1)'/Nbdo * 2 * pi;

if strcmp(ishape,'ellipse')
xi = ai*cos(ti); yi = bi*sin(ti);
xo = ao*cos(to); yo = bo*sin(to);
elseif strcmp(ishape,'super')
xi = abs(cos(ti)).^(2/n).*ai.*sign(cos(ti));
yi = abs(sin(ti)).^(2/n).*bi.*sign(sin(ti));
xo = abs(cos(to)).^(2/n).*ao.*sign(cos(to));
yo = abs(sin(to)).^(2/n).*bo.*sign(sin(to));  
end



figure(1); clf;
plot([xi;xi(1)], [yi;yi(1)],'k','linewidth',2)
hold on
plot([xo;xo(1)], [yo;yo(1)],'k','linewidth',2)
axis equal

NbdP = 32; Dpost = 1; Gx = 1; Gy = 1; 
Xpost = [Dpost/2*cos(-(0:NbdP-1)'/NbdP*2*pi);...
      Dpost/2*sin(-(0:NbdP-1)'/NbdP*2*pi)];  

nrow = floor((WsizeY-Gy)/(Dpost+Gy));
ncol = floor((ai+Gx)/(Dpost+Gx));
H = WsizeY;
L = ai; 
period = ncol; 
epsilon = 1/period;
delLat = (Gy+Dpost)*epsilon;

pnrow = ceil(3*nrow);
pHeight = (pnrow-1)*(Gy+Dpost);

% FIRST THE TOP PART
centy1stCol = linspace(-pHeight/2,pHeight/2,pnrow)' + bi + WsizeY/2;
centx = linspace(-ai/2,ai/2,ncol);

centx = centx(ones(pnrow,1),:);

delLatVect = [0:ncol-1]*delLat;
centy = delLatVect(ones(pnrow,1),:) + centy1stCol(:,ones(1,ncol));

Xobs = zeros(2*NbdP,pnrow*ncol);
for iwall = 1 : pnrow*ncol
  Xobs(:,iwall) = [Xpost(1:end/2)+centx(iwall); Xpost(end/2+1:end)+centy(iwall)];    
end

jdx = 1; XwallsInt = [];
for iwall = 1 : pnrow*ncol
    xwall = Xobs(1:end/2,iwall); ywall = Xobs(end/2+1:end,iwall);
    idx1 = isempty(find(ywall <= 1.02*bi,1));
    idx2 = isempty(find(ywall >= 0.98*bo,1));
    if idx1 && idx2
      XwallsInt(:,jdx) = [xwall;ywall];
      jdx = jdx + 1;
    end
end


% THEN THE BOTTOM PART
centy1stCol = linspace(-pHeight/2,pHeight/2,pnrow)' - bi - WsizeY/2;
% centx = linspace(-ai/2,ai/2,ncol);
centx = linspace(ai/2,-ai/2,ncol);

centx = centx(ones(pnrow,1),:);

delLatVect = [0:ncol-1]*delLat;
centy = delLatVect(ones(pnrow,1),:) + centy1stCol(:,ones(1,ncol));

Xobs = zeros(2*NbdP,pnrow*ncol);
for iwall = 1 : pnrow*ncol
  Xobs(:,iwall) = [Xpost(1:end/2)+centx(iwall); Xpost(end/2+1:end)+centy(iwall)];    
end


for iwall = 1 : pnrow*ncol
    xwall = Xobs(1:end/2,iwall); ywall = Xobs(end/2+1:end,iwall);
    idx1 = isempty(find(ywall >= -1.02*bi,1));
    idx2 = isempty(find(ywall <= -0.98*bo,1));
    if idx1 && idx2
      XwallsInt(:,jdx) = [xwall;ywall];
      jdx = jdx + 1;
    end
end

figure(1); 
plot([XwallsInt(1:end/2,:);XwallsInt(1,:)],[XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)],'k','linewidth',2)





