clear all
close all
clc
load '../data_export/Python_Sims_FW_TT_125mph_test'

P = TT_Sim.P.Drive;
t = TT_Sim.t;

P(P < 0) = 0;

t2 = t(1):0.2:t(end);

P = interp1(t, P, t2);
t = t2;

plot(t,P)

figure

plot(diff(t))

Pprofile = [t' P'];

clear t P t2

save -7 'TT_Pprofile_125mph_test.mat' 