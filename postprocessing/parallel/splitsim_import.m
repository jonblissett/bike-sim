files = glob('*_part_*.mat');
S = [];
qu = length(files);
for i = 1:qu
  fname = cell2mat(files(i));
  data = load(fname);
  fn = fieldnames(data.Sims);
  for j = 1:length(fn)
    if !strcmp(cell2mat(fn(j)),'ALL') & !strcmp(cell2mat(fn(j)),'names')
      if i == 1
        eval(['S.' cell2mat(fn(j)) '=cell2mat(data.Sims.' cell2mat(fn(j)) ');']);
      else
        eval(['S.' cell2mat(fn(j)) '=vertcat(S.' cell2mat(fn(j)) ',cell2mat(data.Sims.' cell2mat(fn(j)) '));']);
      end
    end
  end
  %sprintf('Imported %.0f%%', 100*i/qu)
end
sprintf('Import complete')


file = glob('*_variables_*.mat');

for k = 1:length(file)
  fname = cell2mat(file(k));
  load(fname);

  fn = fieldnames(list);
  vars = '';
  for j = 1:length(fn)
    name = cell2mat(fn(j));;
    eval(['temp = list.' name ';']);
    %sprintf([name '=' num2str(min(temp)) ':' num2str(mean(diff(temp))) ':' num2str(max(temp))]) 
    vars = [vars sprintf('\n%s = %d:%d:%d', name, min(temp), mean(diff(temp)), max(temp))];
  end
  sprintf(vars)
end

file = glob('*_constants_*.mat');

fname = cell2mat(file(1));
load(fname);

Ah = settings.battery.cellAh;

clear -x 'S' 'list' 'vars' 'Ah' 'settings'
save -v7 'temp.mat'
a = struct2cell(S.ERROR);
lvc = zeros(size(a));
for i = 1:length(a)
  lvc(i) = 1 - cell2mat(a(i));
end
lvc = logical(lvc);
%struct_levels_to_print(0)
% try different CdA and battery capacities

rotate3d on
hold on
plot3(S.Pmax(lvc),S.Tmax(lvc),S.Eused(lvc),'+')
plot3(S.Pmax(~lvc),S.Tmax(~lvc),S.Eused(~lvc),'r+')

figure
hold on
plot3(S.Pmax(lvc),S.Tmax(lvc),S.vmax(lvc),'+')
plot3(S.Pmax(~lvc),S.Tmax(~lvc),S.vmax(~lvc),'r+')
xlabel('Peak Power (kW)')
ylabel('Peak Torque (Nm)')
%clc

figure()
hold on
for mass = unique(S.N2)'
  srch = (lvc' & (S.N2 == mass) & range_filter(S.dTbat, 0, 45));
  [mi, im] = min(S.tmax(srch));%min(S.Eused(srch));
  list = 1:length(S.Eused);
  list = list(srch);
  im = list(im);
  fprintf('Best t for %.0fkg bike is %.1fs with; max=%.0f mph, %dNm, %dkW, %.0f rpm, N=%d/%d, E=%.0f Wh, Euse/Erat=%.3f,dT=%.1f°, Ns=%.0f\n',S.m(im)-90,S.tmax(im),S.vmax(im)*2.23,S.Tmax(im),S.Pmax(im),S.RPMmax(im),S.N2(im),S.N1(im),S.Eused(im),S.Eused(im)./S.Erated(im),S.dTbat(im),S.Ns(im))
  fprintf('Motor L=%.0f, N=%.1f,D=%d,%dAh\n',S.Lcore(im), S.turns(im), S.drives(im),Ah*S.Np(im))
  %sprintf('%d,%.1f,%.0f,%d,%d,%d,%d/%d,%.0f,%.3f,%.1f',S.m(im)-90,37.73*3600/S.tmax(im),S.vmax(im)*2.23,S.Tmax(im),S.Pmax(im),S.RPMmax(im),S.N2(im),S.N1(im),S.Eused(im),S.Eused(im)./S.Erated(im),S.dTbat(im))
  plot(S.N2(im), S.tmax(im),'o')%, S.Tmax(im), S.Eused(im)/10)
  % race number of laps vs battery size for best lap times and total race time.
  % plot bike distances vs time for each battery size (show if any overtake)
  % 8 laps
end
sprintf('Make battery energy rating based on 5C test')
ylabel('Lap time (s)')
xlabel('Battery Series cells')
hold off