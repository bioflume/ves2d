clear all

prams.M = 4*[32 32];
prams.bd = @(ind,m) sampleBd(ind,m,1,'couette');
Flow = {'poiseulle','cubic','couette'};

Ro = 15; Ri = 5;c = Ro*Ri/(Ro^2-Ri^2);

T_POS = @(x,y) [2*x 1-2*y;1-2*y 2*x];
T_CUB = @(x,y) [-6*x*y 3*(x^2+y^2);3*(x^2+y^2) -6*x*y];
T_COU = @(x,y) -2*c*Ro/(x^2+y^2)*[0 1;1 0];

for dd = 3:size(Flow,2)
    flow = Flow{dd};
    prams.bc = @(x) forcing(x,flow);
    domain = fixedBound(prams.M,prams.bd,1);
    
    [U trash mu] = farFieldVel([],[],'confined',prams);
    avgS = avgStress(prams.M,domain,mu,[],[]);
    disp(avgS);
end
% cubic:    avgS = [0 3/400(Ro^4-Ri^4);3/400(Ro^4-Ri^4) 0]
% poisulle: avgS = [0 1;1 0]
% Symbolic derivation of couette flow average stress = 0
% Ro=15; Ri=5;c=Ri*Ro/(Ro^2-Ri^2);
% 
% syms x y
% r = sqrt(x^2+y^2);
% t = atan(y/x);
% ut = c*(Ro/r-r/Ro);
% u = simple(ut*[-sin(t);cos(t)]);
% 
% T = [diff(u,'x') diff(u,'y')];
% T = simple(T+T.');
% T1 = T(1);
% T2 = T(2);
% 
% syms t r
% T1 = simple(subs(subs(T1,'x','r*cos(t)'),'y','r*sin(t)'));
% T2 = simple(subs(subs(T2,'x','r*cos(t)'),'y','r*sin(t)'));
% 
% int(int(T1,'r',Ri,Ro),t,0,2*pi)
% int(int(T2,'r',Ri,Ro),t,0,2*pi)