clear; 
clc;

addpath ../src/
oc = curve;
N = 32;
op = poten(N);

vinf = @(Y) [Y(end/2+1:end); zeros(size(Y(1:end/2)))];

load crashTGtestICsMore

id = 1;
X = Xics(:,id);

X1 = X; % original X
vesicle1 = capsules(X1,[],[],1,1,0);
[Ben1,Ten1,Div1] = vesicle1.computeDerivs;
G1 = op.stokesSLmatrix(vesicle1);
V1 = vinf(X1);
z = V1(1:end/2)+1i*V1(end/2+1:end);
zh = fft(z);
V11 = real(zh); V12 = imag(zh);

cent = oc.getPhysicalCenterShan(X); % centered to origin
X2 = [X(1:end/2)-cent(1); X(end/2+1:end)-cent(2)];
vesicle2 = capsules(X2,[],[],1,1,0);
[Ben2,Ten2,Div2] = vesicle2.computeDerivs;
G2 = op.stokesSLmatrix(vesicle2);
V2 = vinf(X2);
z = V2(1:end/2)+1i*V2(end/2+1:end);
zh = fft(z);
V21 = real(zh); V22 = imag(zh);

V = oc.getPrincAxesGivenCentroid(X, cent);
w = [0;1];
rotation = atan2(w(2)*V(1)-w(1)*V(2), w(1)*V(1)+w(2)*V(2));
X3 = rotationOperator(X, rotation, cent); % rotated 
vesicle3 = capsules(X3,[],[],1,1,0);
[Ben3,Ten3,Div3] = vesicle3.computeDerivs;
G3 = op.stokesSLmatrix(vesicle3);
V3 = vinf(X3);
z = V3(1:end/2)+1i*V3(end/2+1:end);
zh = fft(z);
V31 = real(zh); V32 = imag(zh);

[~,~,length] = oc.geomProp(X);
X4 = X/length * 2; % length doubled 
vesicle4 = capsules(X4,[],[],1,1,0);
[Ben4,Ten4,Div4] = vesicle4.computeDerivs;
G4 = op.stokesSLmatrix(vesicle4);
V4 = vinf(X4);
z = V4(1:end/2)+1i*V4(end/2+1:end);
zh = fft(z);
V41 = real(zh); V42 = imag(zh);

% Build velocity matrices for those vesicles
theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N*exp(1i*theta*ks');
activeModes = [(1:N/2)';(N/2+1:N)'];
B1 = real(basis(:,activeModes)); B2 = imag(basis(:,activeModes));

M1 = ((Div1*G1*Ten1)\eye(vesicle1.N))*Div1;
M1_11 = M1(:,1:end/2); M1_12 = M1(:,end/2+1:end);
Z1_11 = M1_11*B1+M1_12*B2; Z1_12 = M1_12*B1-M1_11*B2;
MV1 = Z1_11*V11 + Z1_12*V12;
ten1 = M1*G1*(-Ben1)*X1;


M2 = ((Div2*G2*Ten2)\eye(vesicle2.N))*Div2;
M2_11 = M2(:,1:end/2); M2_12 = M2(:,end/2+1:end);
Z2_11 = M2_11*B1+M2_12*B2; Z2_12 = M2_12*B1-M2_11*B2;
MV2 = Z2_11*V21 + Z2_12*V22;
ten2 = M2*G2*(-Ben2)*X2;

M3 = ((Div3*G3*Ten3)\eye(vesicle3.N))*Div3;
M3_11 = M3(:,1:end/2); M3_12 = M3(:,end/2+1:end);
Z3_11 = M3_11*B1+M3_12*B2; Z3_12 = M3_12*B1-M3_11*B2;
MV3 = Z3_11*V31 + Z3_12*V32;
ten3 = M3*G3*(-Ben3)*X3;

M4 = ((Div4*G4*Ten4)\eye(vesicle4.N))*Div4;
M4_11 = M4(:,1:end/2); M4_12 = M4(:,end/2+1:end);
Z4_11 = M4_11*B1+M4_12*B2; Z4_12 = M4_12*B1-M4_11*B2;
MV4 = Z4_11*V41 + Z4_12*V42;
ten4 = M4*G4*(-Ben4)*X4;


function Xrot = rotationOperator(X,theta, rotCent)
% Get x-y coordinates
Xrot = zeros(size(X));
x = X(1:end/2)-rotCent(1); y = X(end/2+1:end)-rotCent(2);

% Rotated shape
xrot = (x)*cos(theta) - (y)*sin(theta);
yrot = (x)*sin(theta) + (y)*cos(theta);

Xrot(1:end/2) = xrot+rotCent(1);
Xrot(end/2+1:end) = yrot+rotCent(2);
end % rotationOperator

