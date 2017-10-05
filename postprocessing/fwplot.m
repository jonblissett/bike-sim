figure
plot(fwsim.rpm,[fwsim.Id fwsim.Iq fwsim.Is].*0.707)
grid on
title('Figure 3, Stator current vs Motor speed')
xlabel('Speed (rpm)')
ylabel('Id,Iq,Is (Arms)')
legend('Id', 'Iq', 'Is', 'location', 'southwest')
xlim([0 8400])
ylim([-220 410])

figure
plot(fwsim.rpm,[fwsim.torque fwsim.power/1e3])
grid on
title('Figure 4, Motor Torque vs Speed')
xlabel('Speed (rpm)')
ylabel('Torque (Nm), Power (kW)')
legend('Motor Torque', 'Motor Power', 'location', 'southeast')
ylim([0 250])
xlim([0 8400])