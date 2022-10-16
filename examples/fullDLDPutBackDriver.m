clear; clc;

list_postFile = {'OptVCs5a10PostsNew.mat','OptVCs5a50Posts.mat',...
  'OptVCs4a10Posts.mat','OptVCs5a10p5Posts.mat'};
list_VC = {[5 10], [5 50], [4 10], [5 10]};

runNames{1} = {'fullVCs5a10_Circ','fullVCs5a10_P1','fullVCs5a10_P2','fullVCs5a10_P3'};
runNames{2} = {'fullVCs5a50_Circ','fullVCs5a50_P1','fullVCs5a50_P2','fullVCs5a50_P3'};
runNames{3} = {'fullVCs4a10_Circ','fullVCs4a10_P1','fullVCs4a10_P2','fullVCs4a10_P3'};
runNames{4} = {'fullVCs5a10p5_Circ','fullVCs5a10p5_P1Gy8.5'};

runName = 'testPB';
postFile = []; % empty: circular pillar
whichPost = []; % choose which of the best shapes
% if whichPost > numel(bestIds), then we run with jiggled shapes
if isempty(postFile)
  Dx = 1; 
  Dy = 1;
  Dpost = 1.5;
else
  Dx = []; Dy = []; Dpost = [];  
end

whichCell = []; % empty: both cells, otherwise 1 or 2
theta = atan(1/6);
nRepeatPeriod = 4;
VCs = [5 10];
rhoAL = 5E-4;
NbdOut = 3584;
NbdPill = 64;
Nves = 64;
repulsion = true;
fmm = true;

% TRIANGULAR POST
postFile = 'triXpostSharper.mat';
Dx = 0.75;
Dy = 0.75;

% fullDLDPutBack(runName,postFile,whichPost,Dx,Dy,Dpost,whichCell,...
%     theta,nRepeatPeriod,VCs,rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)

repulsion = true;
NbdPill = 64;
Nves = 64;
rhoAL = 5E-4;


%postFile = 'jiggledTriPostGreaterp10.mat';

fullDLDPutBack('fullVCs5.5a8_triXpost',postFile,[],Dx,Dy,[],[],theta,nRepeatPeriod,...
  [5.5 8],rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)

fullDLDPutBack('fullVCs5a8_triXpost',postFile,[],Dx,Dy,[],[],theta,nRepeatPeriod,...
  [5 8],rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)

%fullDLDPutBackJigg('fullVCs5a10_jiggTri6',postFile,2,[],[],[],[],theta,nRepeatPeriod,...
%  [5 10],rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm);

%fullDLDPutBackJigg('fullVCs5a10_jiggTri7',postFile,3,[],[],[],[],theta,nRepeatPeriod,...
%  [5 10],rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm);

if 0
% VCs = [5; 10] with p = 5
repulsion = true;
NbdPill = 72;
Nves = 64;
rhoAL = 2.5E-4;
theta = atan(1/5);
% CIRCULAR PILLAR
%fullDLDPutBack(runNames{4}{1},[],[],1,1,1.5,[],theta,nRepeatPeriod,list_VC{4},...
%  rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)
% BEST PILLAR
fullDLDPutBack(runNames{4}{2},list_postFile{4},1,0.7,0.85,[],[],...
  theta,nRepeatPeriod,list_VC{4},rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)
end

if 0
% VCs = [5; 10]
repulsion = true;
NbdPill = 72;
Nves = 64;
rhoAL = 5E-4;
%fullDLDPutBack(runNames{1}{1},[],[],1,1,1.5,[],theta,nRepeatPeriod,list_VC{1},...
%  rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)
%for k = 1 : 1
%  fullDLDPutBack(runNames{1}{k+1},list_postFile{1},k,[],[],[],[],...
%    theta,nRepeatPeriod,list_VC{1},rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)
%end
fullDLDPutBack('fullVCs5a10_P3',list_postFile{1},3,[],[],[],[],...
  theta,nRepeatPeriod,list_VC{1},rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm);
end

if 0
% VCs = [5; 50]
repulsion = true;
rhoAL = 5E-4;
fullDLDPutBack('fullVCs5a50_P1from5a10',list_postFile{1},1,[],[],[],[],...
  theta,nRepeatPeriod,list_VC{2},rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)
%for k = 1 : 3
%  fullDLDPutBack(runNames{2}{k+1},list_postFile{2},k,[],[],[],[],...
%    theta,nRepeatPeriod,list_VC{2},rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)
%end
end

if 0
repulsion = true;
rhoAL = 5E-4;
NbdPill = 64;
% VCs = [4; 10]
fullDLDPutBack('fullVCs4a50_P1from5a10',list_postFile{1},1,[],[],[],[],...
  theta,nRepeatPeriod,list_VC{3},rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)
%for k = 2 : 2
%  fullDLDPutBack(runNames{3}{k+1},list_postFile{3},k,[],[],[],[],...
%    theta,nRepeatPeriod,list_VC{3},rhoAL,NbdOut,NbdPill,Nves,repulsion,fmm)
%end
end
