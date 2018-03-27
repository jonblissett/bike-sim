filepath = '';
filestring = '*_part_*.mat';
if (exist ('OCTAVE_VERSION', 'builtin') > 0)
    dfiles = glob([filepath filestring]);
    struct_levels_to_print(0)
    vfile = glob('*_variables_*.mat');
    cfile = glob('*_constants_*.mat');
else
    dfiles = dir([filepath filestring]);
    dfiles = {dfiles.name};
    vfile = dir(strcat(filepath,'*_variables_*.mat'));
    vfile = {vfile.name};
    cfile = dir([filepath,'*_constants_*.mat']);
    cfile = {cfile.name};
end

S = [];
qu = length(dfiles);
for i = 1:qu
  fname = cell2mat(dfiles(i));
  data = load(strcat(filepath,fname));
  fn = fieldnames(data.Sims);
  for j = 1:length(fn)
    if ~strcmp(cell2mat(fn(j)),'ALL') && ~strcmp(cell2mat(fn(j)),'names')
      if i == 1
        eval(['S.' cell2mat(fn(j)) '=cell2mat(data.Sims.' cell2mat(fn(j)) ');']);
      else
        eval(['S.' cell2mat(fn(j)) '=vertcat(S.' cell2mat(fn(j)) ',cell2mat(data.Sims.' cell2mat(fn(j)) '));']);
      end
    end
  end
end
sprintf('Import complete')

for k = 1:length(vfile)
  fname = cell2mat(vfile(k));
  load([filepath fname]);

  fn = fieldnames(list);
  vars = '';
  for j = 1:length(fn)
    name = cell2mat(fn(j));
    eval(['temp = list.' name ';']);
    %sprintf([name '=' num2str(min(temp)) ':' num2str(mean(diff(temp))) ':' num2str(max(temp))]) 
    vars = [vars sprintf('\n%s = %d:%d:%d', name, min(temp), mean(diff(temp)), max(temp))];
  end
  sprintf(vars)
end

fname = cell2mat(cfile(1));
load([filepath fname]);



if (exist ('OCTAVE_VERSION', 'builtin') > 0)
    clear -x 'S' 'list' 'vars' 'settings'
else
    clearvars -except 'S' 'list' 'vars' 'settings'
end

Ah = settings.battery.cellAh;
save -v7 'temp.mat'

a = struct2cell(S.ERROR);
lvc = zeros(size(a));
for i = 1:length(a)
  lvc(i) = 1 - cell2mat(a(i));
end
lvc = logical(lvc);

figure(1)
hold on
for mass = unique(S.Tmax)'
  srch = (lvc' & (S.Tmax == mass) & range_filter(S.dTbat, 0, 45));
  [mi, im] = min(S.tmax(srch));%min(S.Eused(srch));
  list = 1:length(S.Eused);
  list = list(srch);
  im = list(im);
  fprintf('Best t for %.0fkg bike is %.1fs with; max=%.0f mph, %dNm, %dkW, %.0f rpm, N=%d/%d, E=%.0f Wh, Euse/Erat=%.3f,dT=%.1f°, Ns=%.0f\n',S.m(im)-90,S.tmax(im),S.vmax(im)*2.23,S.Tmax(im),S.Pmax(im),S.RPMmax(im),S.N2(im),S.N1(im),S.Eused(im),S.Eused(im)./S.Erated(im),S.dTbat(im),S.Ns(im))
  fprintf('Motor L=%.0f, N=%.1f,D=%d,%dAh\n',S.Lcore(im), S.turns(im), S.drives(im),Ah*S.Np(im))
  %sprintf('%d,%.1f,%.0f,%d,%d,%d,%d/%d,%.0f,%.3f,%.1f',S.m(im)-90,37.73*3600/S.tmax(im),S.vmax(im)*2.23,S.Tmax(im),S.Pmax(im),S.RPMmax(im),S.N2(im),S.N1(im),S.Eused(im),S.Eused(im)./S.Erated(im),S.dTbat(im))
  plot(S.Tmax(im), S.tmax(im),'o')%, S.Tmax(im), S.Eused(im)/10)
end
sprintf('Make battery energy rating based on 5C test')
ylabel('Lap time (s)')
xlabel('Battery Series cells')
hold off