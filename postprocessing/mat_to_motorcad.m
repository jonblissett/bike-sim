clear all
close all

file_exp = '../data_export/splat.dat'

load '../data_export/Python_Sim_import_motor_power_varied_mph_Mr25.mat'

t = TT_Sim.t;
rpm = TT_Sim.Rpm;
torque = interp1(TT_Sim.motor.I_interp,TT_Sim.motor.T_interp,TT_Sim.Iq);

subplot(2,1,1)
plot(t,rpm)
xlabel('Time (s)')
ylabel('Motor RPM')
subplot(2,1,2)
plot(t,torque)
xlabel('Time (s)')
ylabel('Motor Torque')

fi = fopen (file_exp, 'w');

t_in = t(1):0.5:t(end);
rpm_in = interp1(t,rpm,t_in);
torque_in = interp1(t,torque,t_in);

for i = 1:length(t_in)
  fprintf(fi, '%.2f\t%.0f\t%.0f\n',t_in(i),torque_in(i),rpm_in(i));
endfor

fclose(fi);



%save -7 'TT_Pprofile_125mph_test.mat' 