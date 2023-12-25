% LOAD THE DATA one by one or all

load selfTensionData_1.mat % set 1
%load selfTensionData_2.mat % set 2
%load selfTensionData_3.mat % set 3
%load selfTensionData_4.mat % set 4

% Each .mat file includes
% - XstandStore, matrix of 256 x 26485 storing x and y coordinates of nsampInSet = 26485
% vesicles discretized with 128 points, each column is [x;y] -- this is the
% input to the network

% - selfTenStore is a matrix of 128 x 26485 
% nearVelocity, matrix of 512 x 5 x 25045 storing (scalar) tension of vesicle for 26485 vesicles
% discretized with 128 points

% nsampInSet can be different in the last set (4)

% In total there are ~ 106K data

%% 
% Input is a vesicle configuration stored in XstandStore 
% Output is tension stored in selfTenStore 


