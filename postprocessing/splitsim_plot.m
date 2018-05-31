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
x = [];
y = [];
dTlim = 30:5:60;
for j = 1:length(dTlim)
  i = 0;
  for mass = 162%unique(S.Ns)'
    i = i + 1;
    srch = (lvc' & (S.Ns == mass) & range_filter(S.dTbat, 0, dTlim(j)) & (S.RPMmax < 11e3));
    [mi, im] = min(S.tmax(srch));%min(S.Eused(srch));
    list = 1:length(S.Eused);
    list = list(srch);
    im = list(im);
    fprintf('Best t for %.0fkg bike is %.1fs with; max=%.0f mph, %dNm, %dkW, %.0f rpm, N=%d/%d, E=%.0f Wh, Euse/Erat=%.3f,dT=%.1f°, Ns=%.0f\n',S.m(im)-90,S.tmax(im),S.vmax(im)*2.23,S.Tmax(im),S.Pmax(im),S.RPMmax(im),S.N2(im),S.N1(im),S.Eused(im),S.Eused(im)./S.Erated(im),S.dTbat(im),S.Ns(im))
    fprintf('Motor L=%.0f, N=%.1f,D=%d,%dAh\n',S.Lcore(im), S.turns(im), S.drives(im),Ah*S.Np(im))
    %sprintf('%d,%.1f,%.0f,%d,%d,%d,%d/%d,%.0f,%.3f,%.1f',S.m(im)-90,37.73*3600/S.tmax(im),S.vmax(im)*2.23,S.Tmax(im),S.Pmax(im),S.RPMmax(im),S.N2(im),S.N1(im),S.Eused(im),S.Eused(im)./S.Erated(im),S.dTbat(im))
    %plot(S.Ns(im), S.tmax(im),'o')%, S.Tmax(im), S.Eused(im)/10)
    if isempty(S.tmax(im))
      x(i,j) = NaN;
      y(i,j) = NaN;
    else
      x(i,j) = S.Ns(im);
      y(i,j) = S.tmax(im);
    end
    % race number of laps vs battery size for best lap times and total race time.
    % plot bike distances vs time for each battery size (show if any overtake)
    % 8 laps
  end
end
figure()
plot(x, y)
hold on
sprintf('Make battery energy rating based on 5C test')
ylabel('Lap time (s)')
xlabel('Battery Series cells')
hold off