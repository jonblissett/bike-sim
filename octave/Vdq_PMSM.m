function [Vs, Vd, Vq, PF] = Vdq_PMSM(Ke,P,Rs,Ld,Lq,Id,Iq,w_mech)
% Calculate PMSM voltages with dq model, assuming dI/dt~=0
% Ke in V/electrical rad/s
% P = pole pairs
    w = w_mech*(P/2);
    Vd = Rs*Id-w.*Lq.*Iq;
    Vq = Rs*Iq+w.*(Ld.*Id+Ke);
    Vs = sqrt(Vd.^2 + Vq.^2);
    %  
    Valpha = atan(Vd./Vq);
    Ialpha = atan(Id./Iq);
    PF = Ialpha - Valpha;
   
    PF(isnan(PF)) = 1;
end

%[Vs, Vd, Vq] =
%Vdq_PMSM(0.34688/6,12,0.007313,53e-6,61e-6,-133,921,9000/30*pi);    %
%Example for Parker GVM210-150-V6 motor, Id,Iq are peak of sine values
