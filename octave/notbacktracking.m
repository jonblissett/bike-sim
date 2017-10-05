clear all
close all

load '../data_export/TT_SpeedLimits/SpeedLimitsTest.mat'

[n y] = size(T); % number of corners

%same = (1:n)(!same);
same = (1:n)(change);

T = T(same,:);
E = E(same,:);

TE = T./E;

% find best time saving * energy saving??

%score each on timeSaving/laptime*energySaving/LapEnergy