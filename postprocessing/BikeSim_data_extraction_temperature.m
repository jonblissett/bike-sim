file = 'Python_TT_sim_UoN01_max';
load(file);

fnames = fieldnames(TT_sim_UoN01_max);

keep = [5,9,11,16,19,21];
remove = 1:length(fnames);


for i = 1:length(keep)
  remove = remove(remove~=keep(i));
end

TT_sim_UoN01_max = rmfield(TT_sim_UoN01_max,fnames(remove));

plot(TT_sim_UoN01_max.t, [TT_sim_UoN01_max.P.DriveLosses TT_sim_UoN01_max.P.MotorLosses TT_sim_UoN01_max.Rpm])
xlabel('Time (s)')
ylabel('Losses (W)')
legend('Drive (IGBTs)','Motor (total)','Motor RPM')
 
file = [file '_losses.mat']
clear -x TT_sim_UoN01_max file
save -mat7-binary file