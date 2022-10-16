clear;

load triXpostSharper
Dx = 1;
Dy = 1;
runName = 'streamVCs10a1_triPost';
folderName = './output/streamVCs10a1_triPost/';
wbox = 2.5;
theta = atan(1/6);
Next = 3072;
Nint = 72;
Nves = 64;
Ufar = 1.2;
VCmaj = 1;
VCmin = 10;
seedRate = 1;
useFMM = true;
totnv = 100;
%Xpost = [];

DLDStreamRunsWpost(runName,folderName,Xpost,wbox,Dx,Dy,theta,...
    Next,Nint,Nves,Ufar,VCmaj,VCmin,seedRate,useFMM,totnv)
