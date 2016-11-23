load 'PikesPeakMap.mat'
close all
%PikesPeak.t = (0:2:768)'-8;
PikesPeak.h = PikesPeak.z - 2884.6;

pkg load all

locsMin = [4, 22, 26, 31, 33, 40, 43, 45, 59, 64, 77, 88, 95, 101, 104, 111, 119, 123, 127, 133, 137, 143, 146, 151, 160, 180, 192, 195, 205, 207, 214, 222, 226, 236, 238, 247, 258, 262, 272, 276, 296, 299, 308, 313, 316, 330, 337, 340, 355, 359, 366, 379, 384];
locsMin = locsMin + 1;
plot(PikesPeak.v)
hold on
scatter(locsMin, PikesPeak.v(locsMin))
ylabel 'Speed (m/s)'
xlabel 'Data index'
hold off

figure();
plot(PikesPeak.dist/1609.34,[PikesPeak.h PikesPeak.h(1)+cumtrapz(gradient(PikesPeak.h,PikesPeak.dist))]);
%figure();
%plot3(PikesPeak.x,PikesPeak.y,PikesPeak.z,PikesPeak.x,PikesPeak.y,zeros(size(PikesPeak.h)),'--')

%duplicates = diff(PikesPeak.dist) == 0;
%list = 1:384;
%list = list(duplicates) + 1;
%PikesPeak.dist2 = PikesPeak.dist;
%PikesPeak.dist2(duplicates) = (PikesPeak.dist2(list)-PikesPeak.dist2(list+1))/2+PikesPeak.dist2(list);

not_duplicates = diff(PikesPeak.dist) ~= 0;
dist_not_dup = PikesPeak.dist(not_duplicates);
%%
d_spacing = 1/10*length(dist_not_dup);%66;%round(PikesPeak.dist(end)/245);
d_even = linspace(dist_not_dup(1),dist_not_dup(end),d_spacing)';
h_even = interp1(dist_not_dup,PikesPeak.h(not_duplicates),d_even,'cubic');
figure();
plot(PikesPeak.dist,PikesPeak.h,d_even,h_even)

%grad = diff([PikesPeak.h' PikesPeak.h(1)]')./PikesPeak.dDist;
grad = gradient(PikesPeak.h,PikesPeak.dist);
angle = 180/pi*atan(grad);

%grad_even = diff([h_even' h_even(1)]')/d_spacing;
grad_even = gradient(h_even,d_even);%gradient(filtfilt(b,a,h_even));
angle_even = 180/pi*atan(grad_even);

[b, a]=butter(3, 0.1);                  % 5 Hz low-pass filter
filtangle = filtfilt(b,a,angle);
filtgrad = filtfilt(b,a,grad_even);%tan(angle_even*pi/180);

figure();
plot(PikesPeak.dist,angle,d_even,angle_even,PikesPeak.dist,filtangle)

figure();
plot(PikesPeak.dist,PikesPeak.h,d_even,PikesPeak.h(1)+cumtrapz(d_even,grad_even),d_even,PikesPeak.h(1)+cumtrapz(d_even,filtgrad))

PikesPeak.gradient = interp1(d_even,grad_even,PikesPeak.dist,'cubic');
PikesPeak.filtgradient = interp1(d_even,filtgrad,PikesPeak.dist,'cubic');

figure();
plot(PikesPeak.dist,[PikesPeak.h PikesPeak.h(1)+cumtrapz(PikesPeak.dist,PikesPeak.gradient) PikesPeak.h(1)+cumtrapz(PikesPeak.dist,PikesPeak.filtgradient)])

%PikesPeak.Distance = PikesPeak.dist;
PikesPeak.heading = zeros(size(PikesPeak.dist));
PikesPeak.lean = zeros(size(PikesPeak.dist));


% resample at 20Hz
PikesPeakRef.t = linspace(min(PikesPeak.t),max(PikesPeak.t),(max(PikesPeak.t)-min(PikesPeak.t))*20+1);
PikesPeakRef.v = interp1(PikesPeak.t,PikesPeak.v,PikesPeakRef.t);
PikesPeakRef.Distance = interp1(PikesPeak.t,PikesPeak.dist,PikesPeakRef.t);
PikesPeakRef.Iq = zeros(size(PikesPeakRef.t));
PikesPeakRef.vGPS = zeros(size(PikesPeakRef.t));
PikesPeakRef.lean = zeros(size(PikesPeakRef.t));
locsMin2 = int32(round(interp1(PikesPeakRef.t,1:length(PikesPeakRef.t),PikesPeak.t(locsMin))));

%clearvars -except PikesPeak PikesPeakRef locsMin locsMin2
%save 'PikesPeakData.mat'

%clear -x PikesPeak PikesPeakRef locsMin locsMin2
%save -mat7-binary 'PikesPeakData.mat'