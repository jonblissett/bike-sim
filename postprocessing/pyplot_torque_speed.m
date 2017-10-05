close all
hold on
plot(PP_sim.v*2.23, PP_sim.Iq*PP_sim.constants.Km.*PP_sim.Rpm*pi/30/750,[0 PP_sim.v'*2.23],[0 PP_sim.v'*2.23*1.53], 'r')
xlim([0 170])
ylim([0 200])
xlabel('Speed (mph)')
ylabel('Power (hp)')
figure
hold on
plot(PP_sim.v*2.23, PP_sim.Iq, [0 max(PP_sim.v)*2.23], [800 800], 'r')
xlim([0 170])
ylim([0 1200])
xlabel('Speed (mph)')
ylabel('Motor current (Iq)')
figure
hold on
plot(PP_sim.t,PP_sim.v*2.23, PP_sim.t,PP_sim.Iq*PP_sim.constants.Km/PP_sim.motor.T_max*100)
xlim([0 42])
xlabel('Time (s)')
ylabel('Speed (mph), Throttle position (%)')

