function [A,B,C,D] = plot_ADK_p(lam,Duration,Energy,w0,Eion,Eion2,eta_crit)
lam = lam*1E-6;
lam_au  = lam/5.29E-11; % wavelength in au
w_L = 2*pi*137/lam_au;
Duration = Duration*1E-15;
Energy = Energy*1E-3;
%w0 = str2double(get(handles.spot_size,'String'));
w0 = w0*1E-6;
Eion = Eion/27.21;
%Eion2 = Eion2/27.21;
%QN = 1./sqrt(2*Eion);
%QN2 = 2./sqrt(2*Eion2);
Ain0 = 2*Energy/(pi*w0^2*Duration)*1E-4; % Peak intensity in W/cm2
E0 = sqrt(Ain0/3.51E16); % Electric field amplitude
Tpul = Duration*4.1322E16/1.66;
T = linspace(-2*Tpul,2*Tpul,1E4);
HT = mean(diff(T));
ft = sqrt(exp(-T.^2/(Tpul)^2)); % need cos(omega*t) term for short pulses
%E=E0*ft*cos(wo*T);
phi = w_L*T;
Ephotmax= 27.21*Eion + 3.17*9.33E-14*Ain0*(lam*1E6)^2*ft.^2; %change this

theta = phi;
for i = 1:length(phi)
    if phi(i)<pi/2
        while theta(i)<pi/2
            theta(i)=theta(i)+pi;
        end
        theta(i)=theta(i)-pi;
    elseif phi(i)>pi/2
        while theta(i)>pi/2
            theta(i)=theta(i)-pi;
        end
    end
    
end

    gam_c = sqrt(2*Eion*w_L^2./(E0*ft).^2);


    a = 1+gam_c.^2-sin(theta).^2;
    b = sqrt(a.^2+4*gam_c.^2.*sin(theta).^2);
    
    c = sqrt((sqrt((b+a)/2)+gam_c).^2+(sqrt((b-a)/2)+sin(abs(theta))).^2);
    PSI = (gam_c.^2+sin(theta).^2+0.5).*log(c)-(3/(2*sqrt(2)))*sqrt(b-a).*sin(abs(theta))-(1/(2*sqrt(2)))*sqrt(b+a).*gam_c;
    rate_exp = exp(-(E0*ft).^2.*PSI./w_L.^3); % exponential part of ionization rate
    % prefactor
    nstr = 1/sqrt(2*Eion); % l^* =0 and m=0 for noble gases
    A_fac = 2.^(2*nstr)/(nstr*gamma(nstr+1)*gamma(nstr));
    kappa = log(gam_c+sqrt(gam_c.^2+1))-gam_c./sqrt(gam_c.^2+1);
    C_fac = sqrt(3*kappa./gam_c.^3); %assumed C=1 from YI paper, probably
    % ok. Varies from C=1 for gamma<<1 to C = 1.2/sqrt(gamma) for gamma>>1. 
    % where C is the Perlomov-Popov-Terennt'ev (PPT) correction 
    rate = A_fac.*Eion.*C_fac.*(2*(2*Eion).^1.5./(E0*ft)).^(2*nstr-1).*rate_exp;
    %rate_QS = A_fac.*Ip.*(2*(2*Ip).^1.5./(E0*ft.*abs(cos(w_L*t+CEP)))).^(2*nstr-1).*exp(-2*(2*Ip).^1.5./(3*E0*ft.*abs(cos(w_L*t+CEP))));
    
    N(1)=1.;

    Rint=0.;
    
    for i = 1:length(phi)-1
            Rint=Rint+rate(i)*HT;
            N(i+1)=exp(-Rint); % Number of Neutrals
            % Neutrals decrease exponentially with time
            N1_YI(i+1)=1.-N(i); % Number of Single Ions          
    end
Ions=N1_YI; % # ionised atoms

if eta_crit<max(Ions)
    [~,b]=min(abs(Ions-eta_crit)); % the time when the ionization fraction is equal to the critical ionization
else
    [~,b]=max(Ephotmax);
end
D = Ain0*ft(b).^2;
A=b;
B=Ions;
C=Ephotmax;