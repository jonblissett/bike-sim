close all
hold on
plot(TT_Sim.v*2.23, TT_Sim.Iq*TT_Sim.constants.Km.*TT_Sim.Rpm*pi/30/750,[0 TT_Sim.v'*2.23],[0 TT_Sim.v'*2.23*1.53], 'r')
xlim([0 170])
ylim([0 200])
xlabel('Speed (mph)')
ylabel('Power (hp)')
figure
hold on
plot(TT_Sim.v*2.23, TT_Sim.Iq, [0 max(TT_Sim.v)*2.23], [800 800], 'r')
xlim([0 170])
ylim([0 1200])
xlabel('Speed (mph)')
ylabel('Motor current (Iq)')
figure
hold on
plot(TT_Sim.t,TT_Sim.v*2.23, TT_Sim.t,TT_Sim.Iq*TT_Sim.constants.Km/TT_Sim.motor.T_max*100)
xlim([0 42])
xlabel('Time (s)')
ylabel('Speed (mph), Throttle position (%)')

