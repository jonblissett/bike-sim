
Ipk = 987.1;
Tpk = 308.64;
kt = 0.366;
x2 = 493.6;   % stall current continuous
y2 = 178.47;  % stall torque continiuous

Ipk = 697;
Tpk = 256.96;
kt = 0.425;
x2 = 348.5;   % stall current continuous
y2 = 146.67;  % stall torque continiuous
% let's try another cubic!
%x1 = Ipk;
%y1 = Tpk;
%b = (y2+x2^3/x1^3*(kt*x1-y1)-kt*x2)/(x2^2-x2^3/x1);
%a = 1/x1^3*(y1-kt*x1-b*x1^2);
b = (y2+x2^3/Ipk^3*(kt*Ipk-Tpk)-kt*x2)/(x2^2-x2^3/Ipk);
a = 1/Ipk^3*(Tpk-kt*Ipk-b*Ipk^2);

I = linspace(0,Ipk,100);
T = a*I.^3+b*I.^2+kt*I;

Ip = I*sqrt(2);
%T2 = a/2/sqrt(2)*Ip.^3+b/2*Ip.^2+kt/sqrt(2)*Ip;
a = a/2/sqrt(2);
b = b/2;
c = kt/sqrt(2);

figure
plot(I,T,I,I.*kt,[0 Ipk], [0 Tpk])
grid on
xlim([0 1100])
xlabel('Current (rms)')
ylabel('Torque (Nm)')


figure
plot(Ip,T,Ip,Ip.*kt/sqrt(2),[0 Ipk*sqrt(2)], [0 Tpk])
grid on
xlim([0 1400])
xlabel('Current (pk-sin Amps)')
ylabel('Torque (Nm)')

%% Assume: T = aI^2 + bx 
%a = (Tpk-kt*Ipk)/Ipk^2;
%b = kt;
%
%I = linspace(0,Ipk,100);
%
%T = a*I.^2 + b*I;
%
%plot(I,T,I,I.*kt)
%grid on
%xlim([0 1100])
%

%
%% let's try a cubic
%
%
%x2 = 420;
%
%a = (2*y1-2*kt*x1)./(2*x1^3-3*x2*x1^2);
%b = -3/2*a*x2;
%
%T = a*I.^3+b*I.^2+kt*I;
%
%figure
%plot(I,T,I,I.*kt)
%grid on
%xlim([0 1100])
%
%
%
%%p = polyfit([0 x1 x2],[0 y1 y2],2);
%%T = polyval(p,I);
%
%figure
%plot(I,T,I,I.*kt)
%grid on
%xlim([0 1100])