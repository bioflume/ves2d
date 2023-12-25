% LOAD THE DATA one by one or all

load nearFieldData_1.mat % set 1
%load nearFieldData_2.mat % set 2
%load nearFieldData_3.mat % set 3
%load nearFieldData_4.mat % set 4

% Each .mat file includes
% - XstandStore, matrix of 512 x 25045 storing x and y coordinates of 25045
% vesicles discretized with 256 points, each column is [x;y] -- this is the
% input to the network

% - nearVelocity, matrix of 512 x 5 x 25045 storing velocity in x and y
% coordinates at 5 layers surrounding the vesicle for 25045 vesicles
% discretized with 256 points, each column is [vx; vy]. The points where
% the velocity is calculated are stored in tracersXstore which has the same
% size as nearVelocity. The layers are (1) on the vesicle, (2) at h/2, (3)
% at h, (4) at 1.5 * h, (5) at 2*h with h is the arc-length spacing equal
% to 1/512. 

% In total there are ~ 100K data

% you can visualize the velocity as 

nves = size(XstandStore,2);

for ives = 1 : nves
figure(1); clf; hold on;

plot(XstandStore(1:end/2,ives), XstandStore(end/2+1:end,ives),'k','linewidth',2)
plot(tracersXstore(1:end/2,:,ives),tracersXstore(end/2+1:end,:,ives),'ro','markersize',8)
quiver(tracersXstore(1:end/2,:,ives),tracersXstore(end/2+1:end,:,ives), nearVelocity(1:end/2,:,ives),nearVelocity(end/2+1:end,:,ives))

axis equal

pause
end


%% 
% If we want to reduce the input size to 256 from 512, we need to
% downsample everything to 128 points on vesicle, that can be done as
% follows
Xdown = zeros(256,nves);
velDown = zeros(256,5,nves);

for ives = 1 : nves
Xdown(:,ives) = [interpft(XstandStore(1:end/2,ives),128);interpft(XstandStore(end/2+1:end,ives),128)];
velDown(:,ives) = [interpft(nearVelocity(1:end/2,ives),128);interpft(nearVelocity(end/2+1:end,ives),128)];
end