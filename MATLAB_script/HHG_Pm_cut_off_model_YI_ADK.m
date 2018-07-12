% Calculates HHG phasematching pressure and cut off photon energy
% See "Harmonic_generation_GUI_manual.pdf" for details on the calculation

clear all
% Input parameters
% handles.values left over from GUI format
handles.target_photon = 40; % eV; problems if > cut off
handles.pulse_energy = 0.25; % mJ
handles.pulse_duration = 100; % fs
handles.spot_size = 27.5; % um (1/e^2 rad.)
handles.Laser_wavelength = 1.3; % um

% choose ionization model
Ionize = 'YI'; % 'ADK' or 'YI'
% choose gas species
Gas = 'Kr'; % 'Ar','Kr','Xe','Ne','He' or 'H2'

lam = handles.Laser_wavelength;
lam = lam*1E-6; %convert back to SI

% constants
N_0 = 2.504E25; %Standard number density
r_e = 2.82E-15; %Classical electron radius
c = 3E8; % speed of light
me = 9.11E-31; % electron mass
e = 1.6E-19; % electron charge
ep0 = 8.85E-12; % vacuum permitivity
h = 6.63E-34; % Planck's constant


% gas properties: complex refractive index and Ip
abs_data = load([Gas '_300mb_1mm.csv']);

if strcmp(Gas,'Ar')
    n_gas = 1+ 1E-6*(67.86711+30182.943*(lam*1E6).^2./(144*(lam*1E6).^2-1));
    E1 = 15.76; % Ionization potential in eV
    E2 = 27.65; % 2+ Ionization potential in eV
elseif strcmp(Gas,'Ne')
    n_gas = 1+0.001205*(0.1063*(lam*1E6).^2./(184.661*(lam*1E6).^2-1)+182.90*(lam*1E6).^2/(376.84*(lam*1E6)^2-1));
    E1 = 21.56;
    E2 = 40.10;
elseif strcmp(Gas,'He')
    n_gas = 1+0.01470091/(423.98*(lam*1E6).^2-1);
    E1 = 24.59;
    E2 = 54.42;
elseif strcmp(Gas,'Kr')
    n_gas = 1+0.012055*(0.2104*(lam*1E6)^2/(65.4742*(lam*1E6)^2-1)+0.227*(lam*1E6)^2/(73.698*(lam*1E6)^2-1)+5.14975*(lam*1E6)^2/(181.08*(lam*1E6)^2-1));
    E1 = 13.99;
    E2 = 24.36;
elseif strcmp(Gas,'Xe')
    n_gas = 1+0.012055*(0.26783*(lam*1E6)^2/(46.301*(lam*1E6)^2-1)+0.29841*(lam*1E6)^2/(50.578*(lam*1E6)^2-1)+5.0333*(lam*1E6)^2/(112.74*(lam*1E6)^2-1));
    E1 = 12.13;
    E2 = 21.21;
elseif strcmp(Gas,'H2')
    n_gas = 1+0.012055*(0.26783*(lam*1E6)^2/(46.301*(lam*1E6)^2-1)+0.29841*(lam*1E6)^2/(50.578*(lam*1E6)^2-1)+5.0333*(lam*1E6)^2/(112.74*(lam*1E6)^2-1));
    E1 = 13.6;
    E2 = 1E3; % no IP2 so make it impossible to ionize
end

enrgy = abs_data(:,1);
T = abs_data(:,2);
alpha_300 = -log(T)/1E-3; % convert CXRO absorption data to attenuation coeffient
del_n = n_gas-1;
eta_crt = (1+N_0 * lam.^2 *r_e./(2*pi*del_n)).^-1; % critical ionization fraction

% functions located at bottom of the script
% Calcuate ionization rate and effective harmonic cut off

if strcmp(Ionize,'YI')
    [b,Ions,Ephotmax,Gen_Int] = plot_YI(handles,E1,eta_crt);
elseif strcmp(Ionize,'ADK')
     [b,Ions,Ephotmax,Gen_Int] = plot_ADK(handles,E1,E2,eta_crt);
else
    print('Wrong string')
    return
end

% calculate phasematching pressure
[E_cut_ef,Pm] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt)

E_eff = sqrt(2*Gen_Int*1E4/(c*ep0)); % effective laser electric field at generation time
Keldysh = (2*pi*c/lam)/(e*E_eff)*sqrt(2*E1*e*me) % Keldysh parameter
E_cut_nm = h*c/(E_cut_ef*e)*1E9 % Cut off in nm

% Yudin-Ivanov ionization calcualtion
% PRA, 64, 013409 (2001)
function [A,B,C,D] = plot_YI(handles,Eion,eta_crit)

    lam = handles.Laser_wavelength;
    lam = lam*1E-6;
    lam_au  = lam/5.29E-11; % wavelength in au
    w_L = 2*pi*137/lam_au;
    Duration = handles.pulse_duration;
    Duration = Duration*1E-15;
    Energy = handles.pulse_energy;
    Energy = Energy*1E-3;
    w0 = handles.spot_size;
    w0 = w0*1E-6;
    Eion = Eion/27.21; % covert to atomic units
    %Eion2 = Eion2/27.21; % no 2+ ionization in YI model right now
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
    Ephotmax= 27.21*Eion + 3.17*9.33E-14*Ain0*(lam*1E6)^2*ft.^2; %instantaneous cut off
    
    % transform frame
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
    
    % YI calculation, see paper
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
    % plot ionization fraction and instantaneous cut-off
    figure(1)   
    clf
    set(0,'DefaultAxesFontSize',22)

    axis tight
    [AX1 ,H1 ,H2] = plotyy(1E15*T/(4.1322E16/1.66),Ions,1E15*T/(4.1322E16/1.66),Ephotmax);
    set(H1,'color',[0.8 0.0 0.1])
    set(H2,'color',[0.4 0.0 0.7])
    set(AX1(1),'ycolor',[0.8 0.0 0.1],'linewidth',4,'xlim',1E15*[-Duration*1.2 Duration*1.2],'ylim',[0 1],'ytick',[0:0.25:1],'xtick',1E15*1.5*linspace(-Duration,Duration,5),'yminortick','on','box','off')
    set(AX1(2),'ycolor',[0.4 0.0 0.7],'linewidth',4,'xlim',1E15*[-Duration*1.2 Duration*1.2],'ylim',[0 1.1*max(Ephotmax)],'ytick',round(linspace(0,1.2*max(Ephotmax),4)),'xtick',1E15*1.5*linspace(-Duration,Duration,5),'yminortick','on','box','off')
    set(AX1(2), 'XTickLabel','','XAxisLocation','Top')
    hold(AX1(2),'on')

    if eta_crit<max(Ions)
        [~,ind]=min(abs(Ions-eta_crit));
    else
        [~,ind]=max(Ephotmax);
    end

    if eta_crit<max(Ions)
        plot(AX1(2),1E15*T(ind)*ones(1,20)/(4.1322E16/1.66),linspace(0,max(Ephotmax)*1.2,20),'.','markersize',12,'color','k')
        hold(AX1(2),'off')
    else
        hold(AX1(2),'off')
    end    
    xlabel('Time (fs)')
    ylabel(AX1(1),'Ionization Fraction')
    ylabel(AX1(2),'Cut-off energy (eV)')
    
    % outputs
    D = Ain0*ft(ind).^2; %generating intensity
    A=ind; %index of target photon energy
    B=Ions;
    C=Ephotmax;
end

% ADK ionization calculation
% Sov. Phys. JETP 64(6) (1986)
function [A,B,C,D] = plot_ADK(handles,Eion,Eion2,eta_crit)
    lam = handles.Laser_wavelength;
    lam = lam*1E-6;
    Duration = handles.pulse_duration;
    Duration = Duration*1E-15;
    Energy = handles.pulse_energy;
    Energy = Energy*1E-3;
    w0 = handles.spot_size;
    w0 = w0*1E-6;
    Eion = Eion/27.21; % convert to atomic units
    Eion2 = Eion2/27.21;
    QN = 1./sqrt(2*Eion);
    QN2 = 2./sqrt(2*Eion2);
    Ain0 = 2*Energy/(pi*w0^2*Duration)*1E-4; % Peak intensity in W/cm2
    Ain = sqrt(Ain0/3.51E16); % Electric field amplitude
    Tpul = Duration*1E15; % Pulse FWHM in fs
    Tpul = Tpul*1E-15*4.1322E16;
    Tf = Tpul*0.84932;
    lambda0 = lam*1.89E10;
    wo = 2*pi/((lam/3E8)*41.34*1E15);
    TT = 2.*pi/wo;
    Nx = 1E4; % Number of time steps (512)
    HT = 4.*Tf/Nx;
    T_env = [1:Nx]*HT;
    ft = exp(-(T_env-2.0*Tf).^2/Tf^2); % envelope
    for i=1:Nx
        T = i*HT;
        E(i)=Ain*exp(-(T-2.0*Tf)^2/Tf^2)*cos(wo*(T-2*Tf));
        T2=(T-2*Tf)/TT*2.67;
        % Top photon energy
        Ephotmax(i)= 27.21*Eion + 3.17*9.33E-14*Ain0*(lam*1E6)^2*exp(-T2.^2/((Duration*1E15)/1.66).^2); %change this
    end
    N(1)=1.;
    N1(i)=0.;
    N2(1)=0.;
    % Initially all neutral, no ions
    T(1)=0.;
    Rint=0.;
    for i=1:Nx-1
        T(i+1) = i*HT;
        F = abs(E(i))+1.E-30;
        % why?, Keldysh notation
        C2=2^(2*QN)/(QN*gamma(QN+1)*gamma(QN));
        C22=2^(2*QN2)/(QN2*gamma(QN2+1)*gamma(QN2));
        Rate1(i)=Eion*C2*(4.0*Eion*sqrt(2.0*Eion)/F)^(2.0*QN-1.0)*...
        exp(-4.0*Eion*sqrt(2.0*Eion)/(3.0*F));
        % Probability of ionization of first electron per unit time
        Rate2(i)=Eion2*C22*(4.0*Eion2*sqrt(2.0*Eion2)/F)^(2.0*QN2-1.0)*...
        exp(-4.0*Eion2*sqrt(2.0*Eion2)/(3.0*F));
        % Probability of ionization of second electron per unit time
        Rint=Rint+Rate1(i)*HT;
        N(i+1)=exp(-Rint); % Number of Neutrals
        % Neutrals decrease exponentially with time
        N1(i+1)=1.-N(i)-N2(i); % Number of Single Ions
        % Everthing minus # neutrals and # doubles
        N2(i+1)= N2(i)+Rate2(i)*N1(i)*HT; % Number of Double Ions
        % # Doubles from last time plus # newly ionised singles
        Rate1(i)=Rate1(i)*N(i);
        Rate2(i)=Rate2(i)*N1(i);
    end
    Ions=1-N; % # ionised atoms
    Ion1=N1;
    Ion2=N2;
    T=(T-2*Tf)/TT*2.67;
    % plot ionization fraction and instantaneous cut-off
    figure(1)
    clf
    set(0,'DefaultAxesFontSize',22)

    [AX1 ,H1 ,H2] = plotyy(T,Ions,T,Ephotmax);
    set(H1,'color',[0.8 0.0 0.1])
    set(H2,'color',[0.4 0.0 0.7])
    set(AX1(1),'ycolor',[0.8 0.0 0.1],'linewidth',4,'xlim',1E15*[-Duration*1.2 Duration*1.2],'ylim',[0 1],'ytick',[0:0.25:1],'xtick',[-40:20:40],'yminortick','on','box','off')
    set(AX1(2),'ycolor',[0.4 0.0 0.7],'linewidth',4,'xlim',1E15*[-Duration*1.2 Duration*1.2],'ylim',[0 1.1*max(Ephotmax)],'xtick',[-40:20:40],'ytick',[0:round(1.2*max(Ephotmax)/100)*10:1.2*max(Ephotmax)],'yminortick','on','box','off')
    set(AX1(2), 'XTickLabel','','XAxisLocation','Top')
    hold(AX1(2),'on')
    if eta_crit<max(Ions)
        [~,ind]=min(abs(Ions-eta_crit));
    else
        [~,ind]=max(Ephotmax);
    end

    if eta_crit<max(Ions)
        plot(AX1(2),T(ind)*ones(1,20),linspace(0,max(Ephotmax)*1.2,20),'.','color','k')
        hold(AX1(2),'off')
    else
        hold(AX1(2),'off')
    end    
    xlabel('Time (fs)')
    ylabel(AX1(1),'Ionization Fraction')
    ylabel(AX1(2),'Cut-off energy (eV)')
    
    %outputs
    D = Ain0*ft(ind).^2; % effective intensity at generation time
    A=ind; %index of target photon energy
    B=Ions; % ionization fraction
    C=Ephotmax;
end

% calculate phasematching pressure and equivalent number density
function [E_cut_eff, Pm] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt)
        grey = 0.5*ones(1,3);
        w0 = handles.spot_size;
        w0 = w0*1E-6;
        lam = handles.Laser_wavelength;
        lam = lam*1E-6;
        E_cut_eff = Ephotmax(b);
        eta_up = Ions(1:b-1);
        Ephot_up = Ephotmax(1:b-1);
        Phs_mtch_prsr = 1013*lam^2./(2*pi^2*w0.^2*del_n*(1-eta_up/eta_crt));
        
        [~,ind] = min(abs(Ephotmax(1:(length(Ephotmax)/2))-handles.target_photon));
        Pm = Phs_mtch_prsr(ind);
        
        
        V_foc = pi^2 *w0.^4/(2*lam);
        nmbr_dnsty = 100*(Phs_mtch_prsr)/(1.38E-23*294); % 100 factor to convert to Pa
        nmbr_emttrs = V_foc .*nmbr_dnsty;
        figure(2)
        clf
        set(0,'DefaultAxesFontSize',22)
        scale = 1E-13;
        [AX,H1,H2] = plotyy(Ephot_up,nmbr_emttrs*scale,Ephot_up,Phs_mtch_prsr);
        set(H1,'color',grey)
        set(H2,'color',grey)
        set(AX(2),'ycolor','k','linewidth',4,'xlim',[min(Ephot_up) 1.1*max(Ephot_up)],'ylim',[min(Phs_mtch_prsr)*0.8 max(Phs_mtch_prsr)*1.1],'xtick',linspace(0,1.1*max(Ephot_up),5),'xminortick','on','yminortick','on','box','off')
        %set(AX(1),'ycolor','k','linewidth',4,'xlim',[min(Ephot_up) 1.1*max(Ephot_up)],'ylim',[min(nmbr_emttrs)*0.8 1.1*max(nmbr_emttrs)]*scale,'xtick',round(linspace(0,1.1*max(Ephot_up),5)),'ytick',round(scale*linspace(min(nmbr_emttrs)*0.8,1.1*max(nmbr_emttrs),5),0),'xminortick','on','yminortick','on','box','off')
        set(AX(1),'ycolor','k','linewidth',4,'xlim',[min(Ephot_up) 1.1*max(Ephot_up)],'ylim',[min(nmbr_emttrs)*0.8 1.1*max(nmbr_emttrs)]*scale,'xtick',round(linspace(0,1.1*max(Ephot_up),5)),'xminortick','on','yminortick','on','box','off')
        set(AX(2), 'XTickLabel','','XAxisLocation','Top')
        ylabel(AX(2),'P_m (mbar)')
        ylabel(AX(1),'# emitters (\times10^{13})')
        xlabel('Cut-off photon energy (eV)')
end