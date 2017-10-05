clear all

load('TT_Laps_2016.mat')
load('TT_map.mat')
load('Python_sim.mat')
TT_Sim = TT_sim;

TT_Race.v = TT_Race.Rpm/60*19/83*2.003;%*2.236941;

figure()
plot(TT_Sim.Distance,[TT_Sim.v TT_Sim.Iq*K],TT_Race.Distance,[TT_Race.v TT_Race.Iq*K])

TT_Sim.Rpm = TT_Sim.v/(2*pi*TT_Sim.constants.r)*60*TT_Sim.N(1)/TT_Sim.N(2);

sprintf('Simulated lap time = %.1fs = %d:%.1fs\nSimulated lap speed = %.1f mph',TT_Sim.t(end),floor(TT_Sim.t(end)/60),TT_Sim.t(end)-60*floor(TT_Sim.t(end)/60),37.733/TT_Sim.t(end)*3600)
sprintf('Top speed = %.1f mph at %.0f motor rpm',max(TT_Sim.v)*2.237,max(TT_Sim.Rpm))
%    [w, Tmotor, ~] = MotTorqueSpeed(MotTmax,10500/30*pi,MotPmax,10600/30*pi,20);
%

TT_Race.Rdq = 0.006333*(1+0.004041*(TT_Race.TMot-140)); %Rs = 0.005816/2*2/3;

% add losses in cables to motor and DC cables
% Can get DC cable loss from internal resistance?
% plot losses vs weight in copper?
%plot(TT_Race.t,[resLoss 3*TT_Race.Iq.^2/sqrt(2)*R*0.7])
%[~, R] = skin_effect(1.67e-8,0.008,0)


TT_Race.P.motRLoss=TT_Race.Iq.^2.*TT_Race.Rdq;
TT_Race.P.motwLoss=2.6853e-5*TT_Race.Rpm.^2+0.01528*TT_Race.Rpm;
TT_Race.P.MotorLosses=TT_Race.P.motRLoss+TT_Race.P.motwLoss;
[TT_Race.Vs, TT_Race.Vd, TT_Race.Vq, TT_Race.PF] = Vdq_PMSM(0.34688/6,12,0.007313,53e-6,61e-6,0,TT_Race.Iq,TT_Race.Rpm/30*pi);
[TT_Race.P.DriveLosses, ~, ~, ~, ~, TT_Race.Iripple] = inverter_losses(TT_Race.Vdc,TT_Race.Vs,abs(TT_Race.Iq)/1.4146,TT_Race.PF,82e-6,13e3,0.8,1,0.95e-3,0.54e-3,12e-3,25e-3,9.5e-3);
%plot(TT_Race.t,[TT_Race.Vd TT_Race.Vq TT_Race.Vs 100*cos(PF)])
%[Ptotal, Pct, Pcd, Pst, Psd, i_ripple]
%plot(TT_Race.t,[Pct, Pcd, Pst, Psd])
%plot(TT_Race.t,TT_Race.MotorLosses);
%plot(abs((TT_Race.Power-TT_Race.MotorLosses)./TT_Race.Power))

TT_Race.P.Mech = TT_Race.Rpm/30*pi.*TT_Race.Iq*K;
TT_Race.P.Motor = TT_Race.P.Mech+TT_Race.P.MotorLosses;
TT_Race.P.Drive = TT_Race.P.Motor + TT_Race.P.DriveLosses;
TT_Race.Idc = TT_Race.P.Drive./TT_Race.Vdc;
TT_Race.Energy = cumtrapz(TT_Race.t,TT_Race.P.Drive);
TT_Race.Charge = cumtrapz(TT_Race.t,TT_Race.Idc);
TT_Race.P.DCLoss=TT_Race.Idc.^2*0.055;


sprintf('Race motor losses = %.0f Wh, drive losses = %.0f Wh, Total = %.0f Wh',trapz(TT_Race.t,TT_Race.P.MotorLosses)/3600,trapz(TT_Race.t,TT_Race.P.DriveLosses)/3600,trapz(TT_Race.t,TT_Race.P.DriveLosses+TT_Race.P.MotorLosses)/3600)
sprintf('DC (battery + cable IR) losses = %.0f Wh',trapz(TT_Race.t,TT_Race.P.DCLoss)/3600)
sprintf('Race total energy =%.0f Wh',TT_Race.Energy(end)/3600)

% ########### assumptions ###############
TT_Sim.Vdc = linspace(TT_Race.Vdc(1),TT_Race.Vdc(end),length(TT_Sim.t))';
TT_Sim.TMot = linspace(TT_Race.TMot(1),TT_Race.TMot(end),length(TT_Sim.t))';
% ##############

TT_Sim.Rdq = 0.006333*(1+0.004041*(TT_Sim.TMot-140)); %Rs = 0.005816/2*2/3;

% add losses in cables to motor and DC cables
% Can get DC cable loss from internal resistance?
% plot losses vs weight in copper?
%plot(TT_Sim.t,[resLoss 3*TT_Sim.Iq.^2/sqrt(2)*R*0.7])
%[~, R] = skin_effect(1.67e-8,0.008,0)


TT_Sim.P.motRLoss=TT_Sim.Iq.^2.*TT_Sim.Rdq;
TT_Sim.P.motwLoss=2.6853e-5*TT_Sim.Rpm.^2+0.01528*TT_Sim.Rpm;
TT_Sim.P.MotorLosses=TT_Sim.P.motRLoss+TT_Sim.P.motwLoss;
[TT_Sim.Vs, TT_Sim.Vd, TT_Sim.Vq, TT_Sim.PF] = Vdq_PMSM(0.34688/6,12,0.007313,53e-6,61e-6,0,TT_Sim.Iq,TT_Sim.Rpm/30*pi);
[TT_Sim.P.DriveLosses, ~, ~, ~, ~, TT_Sim.Iripple] = inverter_losses(TT_Sim.Vdc,TT_Sim.Vs,abs(TT_Sim.Iq)/1.4146,TT_Sim.PF,82e-6,13e3,0.8,1,0.95e-3,0.54e-3,12e-3,25e-3,9.5e-3);
%plot(TT_Sim.t,[TT_Sim.Vd TT_Sim.Vq TT_Sim.Vs 100*cos(PF)])
%[Ptotal, Pct, Pcd, Pst, Psd, i_ripple]
%plot(TT_Sim.t,[Pct, Pcd, Pst, Psd])
%plot(TT_Sim.t,TT_Sim.MotorLosses);
%plot(abs((TT_Sim.Power-TT_Sim.MotorLosses)./TT_Sim.Power))

TT_Sim.P.Mech = TT_Sim.Rpm/30*pi.*TT_Sim.Iq*K;
TT_Sim.P.Motor = TT_Sim.P.Mech+TT_Sim.P.MotorLosses;
TT_Sim.P.Drive = TT_Sim.P.Motor + TT_Sim.P.DriveLosses;
TT_Sim.Idc = TT_Sim.P.Drive./TT_Sim.Vdc;
TT_Sim.Energy = cumtrapz(TT_Sim.t,TT_Sim.P.Drive);
TT_Sim.Charge = cumtrapz(TT_Sim.t,TT_Sim.Idc);
TT_Sim.P.DCLoss=TT_Sim.Idc.^2*0.055;

sprintf('Sim motor losses = %.0f Wh, drive losses = %.0f Wh, Total = %.0f Wh',trapz(TT_Sim.t,TT_Sim.P.MotorLosses)/3600,trapz(TT_Sim.t,TT_Sim.P.DriveLosses)/3600,trapz(TT_Sim.t,TT_Sim.P.DriveLosses+TT_Sim.P.MotorLosses)/3600)
sprintf('DC (battery + cable IR) losses = %.0f Wh',trapz(TT_Sim.t,TT_Sim.P.DCLoss)/3600)
sprintf('Sim total energy =%.0f Wh',TT_Sim.Energy(end)/3600)

%

mToM = 0.000621371;
D = TT_Race.Distance;
TTsT = interp1(TT_Sim.Distance,TT_Sim.Iq*K,D);
TTsV = interp1(TT_Sim.Distance,TT_Sim.v,D);

h = fdesign.lowpass('N,Fc',16,100,13000);
d = design(h);
TT_Race.v = filtfilt(d.Numerator,1,TT_Race.v);            %zero-phase filtering
TT_Race.Iq = filtfilt(d.Numerator,1,TT_Race.Iq);              %zero-phase filtering

figure();
hold on
[hAx,hLine1,hLine2] = plotyy(D,[TT_Race.Iq*K TTsT],D,[TT_Race.v TTsV]);
xlabel('Distance(m)')%, 'FontSize', 16)
ylabel(hAx(1),'Motor Torque (Nm)')%, 'FontSize', 16)
ylabel(hAx(2),'Speed (m/s)')
%set(hLine2,'Color',[0,0.8,0])
set(hAx(1),'YTick',-50:25:125,'xlim',[-50 2800],'ylim',[-50 125])%),'LineWidth',2
set(hAx(2),'YTick',0:10:60,'xlim',[-50 2800],'ylim',[0 60])%,'ycolor','black','LineWidth',2)
%set(hLine1,'LineWidth',2)
grid(hAx(1),'on')
h_legend=legend('Torque, Sim','Speed, Sim','Torque, Real','Speed, Real','location','south');
%set(h_legend,'FontSize',16);
hold off
%title 'Simulation and recorded data comparison'

figure();
hold on
plot(TT_Race.Distance*mToM,[TT_Race.Iq*K TT_Race.v])
plot(TT_Sim.Distance*mToM,TT_Sim.Iq*K,'c','LineWidth',2)
plot(TT_Sim.Distance*mToM,TT_Sim.v,'g','LineWidth',2)
hold off
%xlim([16 18])
xlabel 'Distance (miles)'
ylabel 'Torque (Nm) & Speed (mph)'
legend('Power, Sim','Speed, Sim','Power, Race','Speed, Race','location','southeast')
title 'Qualification and Race lap comparison'

%%
Ramsey      =	[-2.861e5 -2.837e5; 6.0370e6 6.0405e6];
BrayHill    =   [-2.935e5 -2.905e5; 6.0215e6 6.0245e6];
BrayCorner  =   [-2.938e5 -2.925e5; 6.0217e6 6.0230e6];
Ballacraine =   [-3.015e5 -2.985e5; 6.0260e6 6.0300e6];
Kirkmichael =   [-2.985e5 -2.965e5; 6.0345e6 6.0368e6];
Sulby       =   [-2.970e5 -2.890e5; 6.0370e6 6.0410e6];
All = [];

% vname=@(x) inputname(1);
% toto=pi;
% s=vname(toto);

%track_plot(TT_map.x,TT_map.y,TT_map.z,TT_map.dist,TT_Sim.Distance,TT_Sim.v,300,[])

track_plotXY(TT_map.x,TT_map.y,TT_map.dist,TT_Sim.Distance,TT_Sim.v,300,[0.2*max(TT_Race.v) 1.2*max(TT_Race.v)],Ramsey)
track_plotXY(TT_map.x,TT_map.y,TT_map.dist,TT_Race.Distance,TT_Race.v,300,[0.2*max(TT_Race.v) 1.2*max(TT_Race.v)],Ramsey)

%%
% Remove duplicate distance values
    [~, ia, ~] = unique(TT_Race.Distance);
figure()
color_line3(TT_map.x,TT_map.y,TT_map.z,interp1(TT_Race.Distance(ia),TT_Race.v(ia),TT_map.dist),'linew',6);

    xlabel('East-West (m)','FontSize',14)
    ylabel('North-South (m)','FontSize',14)
    title('TT Zero - Ramsey corner','FontSize',14)
    
    set(gca,'DataAspectRatio',[1 1 1])
xlimit = Ramsey(1,:);
ylimit = Ramsey(2,:);
xlim(xlimit);
ylim(ylimit);



    [~, ia, ~] = unique(TT_Sim.Distance);
figure()
color_line3(TT_map.x/1000,TT_map.y/1000,TT_map.z,interp1(TT_Sim.Distance(ia),TT_Sim.v(ia),TT_map.dist),'linew',6);

    xlabel('East-West (km)','FontSize',14)
    ylabel('North-South (km)','FontSize',14)
    title('TT Zero - Ramsey corner simulation','FontSize',14)
    
        set(gca,'DataAspectRatio',[1 1 1])

xlimit = Ramsey(1,:)/1000;
ylimit = Ramsey(2,:)/1000;
xlim(xlimit);
ylim(ylimit);
%track_plotXY(TT_map.x,TT_map.y,TT_map.dist,TT_Sim.Distance,TT_Sim.TMot,300,[0.2*max(TT_Sim.TMot) 1.2*max(TT_Sim.TMot)],location)
%track_plotXY(TT_map.x,TT_map.y,TT_map.dist,TT_Race.Distance,TT_Race.TMot,300,[0.2*max(TT_Race.TMot) 1.2*max(TT_Race.TMot)],location)
%%
% This creates the 'background' axes
ha = axes('units','normalized', ...
            'position',[0 0 1 1]);
% Move the background axes to the bottom
uistack(ha,'bottom');
% Load in a background image and display it using the correct colors
% The image used below, is in the Image Processing Toolbox.  If you do not have %access to this toolbox, you can use another image file instead.
I=imread('TT_Googlemap.png');
hi = imagesc(I);
colormap gray
% Turn the handlevisibility off so that we don't inadvertently plot into the axes again
% Also, make the axes invisible
set(ha,'handlevisibility','off', ...
            'visible','off')
% Now we can use the figure, as required.
% For example, we can put a plot in an axes
axes('position',[0.3,0.35,0.4,0.4])
plot(rand(10))