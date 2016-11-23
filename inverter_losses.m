% See http:www.iaeng.org./publication./IMECS2010./IMECS2010_pp1401-1406.pdf
% and Infineon 'IGBT Power Losses Calculation Using the Data-Sheet Parameters'
function [Ptotal, Pct, Pcd, Pst, Psd, i_ripple] = inverter_losses(Vbus,Voll,Iorms,PF,L,f_sw,u_ce0,u_d0,r_c,r_d,Eton,Etoff,Ed)
    i_ripple = (Vbus-2^0.5.*Voll).*Voll./(2.*L.*Vbus.*f_sw);
    Iopk = 2^0.5.*Iorms;
    m = 1./(3^0.5).*Voll.*2^0.5./Vbus; % ?????????????
    Pct = u_ce0.*Iopk.*(1./(2.*pi)+m.*cos(PF)./8)+r_c.*Iopk.^2.*(1./8+m.*cos(PF)./(3.*pi));
    Pcd = u_d0 .*Iopk.*(1./(2.*pi)-m.*cos(PF)./8)+r_d.*Iopk.^2.*(1./8-m.*cos(PF)./(3.*pi));
    Idc = Iopk./pi;  % DC equivalent to actual AC output current
    Icon = Idc-i_ripple./2;
    Icoff = Idc+i_ripple./2;
    Pst = (Eton.*Icon+Etoff.*Icoff).*f_sw./300.*Vbus./550;    % Test at Vce=300V, Ic=550A
    Psd = Ed.*f_sw./300.*Vbus./550.*Idc;
    Ptotal = 6.*(Pct+Pcd+Pst+Psd);
end

% e.g.
%[Ptotal, Pct, Pcd, Pst, Psd, i_ripple] = inverter_losses(450,230,350,0.95,75e-6,1e4,0.8,1,0.95e-3,0.54e-3,12e-3,25e-3,9.5e-3)