clear all

addpath ../src
oc = curve;

dt = 1e0;

N = 192;
theta = (0:N-1)'*2*pi/N;
if 0
  X = [cos(theta);3*sin(theta)];
end

if 1
  X = oc.initConfig(N,'curly');
end

vesicle = capsules(X,[],[],1,1);

[Ben,Ten,Div] = vesicle.computeDerivs;
% Bending, tension and divergence operators
op = poten(N);
G = op.stokesSLmatrix(vesicle);
% Stokes SLP
LinSys = [eye(2*N) + dt*G*Ben -dt*G*Ten;Div zeros(N)];
% Time stepping matrix for vesicle suspension




[ra,a,l] = oc.geomProp(X);
% reduced area, area, and length of vesicle
perim = sqrt(a/pi); % circle has same area
%perim = l/2/pi; % circle has soame perimeter

%perim = 1;
Xcirc = [perim*cos(theta);perim*sin(theta)];
vesicleCirc = capsules(Xcirc,[],[],1,1);
% build circular vesicle

[Benc,Tenc,Divc] = vesicleCirc.computeDerivs;
% Bending, tension and divergence operators
Gc = op.stokesSLmatrix(vesicleCirc);
% Stokes SLP

LinSysCirc = [eye(2*N) + dt*Gc*Benc -dt*Gc*Tenc; ...
      Divc zeros(N)];
% scaling matrix

[V,D] = eig(LinSysCirc);
s = find(abs(diag(D)) < 1e-12);
D(s,s) = 1;
LinSysCirc = real(V*D*inv(V));
% remove rank-one null space of scaling matrix





fprintf('Condition # of origional linear system is            %8.2e\n', ...
    cond(LinSys));

fprintf('Condition # of left preconditioned linear system is  %8.2e\n', ...
    cond(inv(LinSysCirc)*LinSys));

fprintf('Condition # of right preconditioned linear system is %8.2e\n', ...
    cond(LinSys*inv(LinSysCirc)));

fprintf('Condition # of scaling matrix is                     %8.2e\n', ...
    cond(LinSysCirc));









