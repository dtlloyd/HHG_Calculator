% WHY ERRORS?

function varargout = main_YI4_tidied(varargin)
 set(0,'defaultlinelinewidth',6)
 % LATEST VERSION
 % notes: MATLAB 2016a compatible 
 % focal volume length = zr/2
 % Ionization by YI model up to +1  charge state
 % Cut-off calc ignores carrier wave
 % now with H2
 % new FoM includes driver intensity scaling, lambda^-9 scaling and w_0^2
 % scaling of pulse energy.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_YI3_OpeningFcn, ...
                   'gui_OutputFcn',  @main_YI3_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main_YI3 is made visible.
function main_YI3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_YI3 (see VARARGIN)

% Choose default command line output for main_YI3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = main_YI3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Gas_type.
function Gas_type_Callback(hObject, eventdata, handles)

contents=get(hObject,'Value');
targ = str2double(get(handles.target_photon,'String'));
lam = str2double(get(handles.Laser_wavelength,'String'));
lam = lam*1E-6; %convert back to SI
N_0 = 2.504E25; %Standard number density
r_e = 2.82E-15; %Classical electron radius

switch contents
    case 1 % Argon
        abs_data = load('Ar_300mb_1mm.csv');
        enrgy = abs_data(:,1);
        T = abs_data(:,2);
        alpha_300 = -log(T)/1E-3;
        n_gas = 1+ 1E-6*(67.86711+30182.943*(lam*1E6).^2./(144*(lam*1E6).^2-1));
        E1 = 15.76;
        E2 = 27.65;
        del_n = n_gas-1;
        eta_crt = (1+N_0 * lam.^2 *r_e./(2*pi*del_n)).^-1;
        [b,Ions,Ephotmax,Gen_Int] = plot_ADK(handles,E1,E2,eta_crt);
        [E_cut_ef,nmbr_ems, E_phot_up] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt);
        %targ = str2double(get(handles.target_photon,'String'));
        if E_cut_ef>=targ
            yes_no = 1;
        else
            yes_no = 0;
        end
        data = struct('yn',yes_no,'ems',nmbr_ems,'Ion',Ions,...
            'Top_phot',Ephotmax,'E_phot_up',E_phot_up,'T_enrgy',enrgy,...
            'alpha_gas',alpha_300,'Gen_Int',Gen_Int);
        hObject.UserData = data;
    case 2 % Neon
        abs_data = load('Ne_300mb_1mm.csv');
        enrgy = abs_data(:,1);
        T = abs_data(:,2);
        alpha_300 = -log(T)/1E-3;
        n_gas = 1+0.001205*(0.1063*(lam*1E6).^2./(184.661*(lam*1E6).^2-1)+182.90*(lam*1E6).^2/(376.84*(lam*1E6)^2-1));
        E1 = 21.56;
        E2 = 40.10;
        del_n = n_gas-1;
        eta_crt = (1+N_0 * lam.^2 *r_e./(2*pi*del_n)).^-1;
        [b,Ions,Ephotmax,Gen_Int]=plot_ADK(handles,E1,E2,eta_crt);
        [E_cut_ef,nmbr_ems,E_phot_up] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt);
        targ = str2double(get(handles.target_photon,'String'));
        if E_cut_ef>=targ
            yes_no = 1;
        else
            yes_no = 0;
        end
        data = struct('yn',yes_no,'ems',nmbr_ems,'Ion',Ions,'Top_phot',Ephotmax,'E_phot_up',E_phot_up,'T_enrgy',enrgy,'alpha_gas',alpha_300,'Gen_Int',Gen_Int);
        hObject.UserData = data;
    case 3 % Helium
        abs_data = load('He_300mb_1mm.csv');
        enrgy = abs_data(:,1);
        T = abs_data(:,2);
        alpha_300 = -log(T)/1E-3;
        n_gas = 1+0.01470091/(423.98*(lam*1E6).^2-1);
        E1 = 24.59;
        E2 = 54.42;
        del_n = n_gas-1;
        eta_crt = (1+N_0 * lam.^2 *r_e./(2*pi*del_n)).^-1;
        [b,Ions,Ephotmax,Gen_Int]=plot_ADK(handles,E1,E2,eta_crt);
        [E_cut_ef,nmbr_ems,E_phot_up] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt);
        targ = str2double(get(handles.target_photon,'String'));
        if E_cut_ef>=targ
            yes_no = 1;
        else
            yes_no = 0;
        end
        data = struct('yn',yes_no,'ems',nmbr_ems,'Ion',Ions,'Top_phot',Ephotmax,'E_phot_up',E_phot_up,'T_enrgy',enrgy,'alpha_gas',alpha_300,'Gen_Int',Gen_Int);
        hObject.UserData = data;
    case 4 % Krypton
        abs_data = load('Kr_300mb_1mm.csv');
        enrgy = abs_data(:,1);
        T = abs_data(:,2);
        alpha_300 = -log(T)/1E-3;
        n_gas = 1+0.012055*(0.2104*(lam*1E6)^2/(65.4742*(lam*1E6)^2-1)+0.227*(lam*1E6)^2/(73.698*(lam*1E6)^2-1)+5.14975*(lam*1E6)^2/(181.08*(lam*1E6)^2-1));
        E1 = 13.99;
        E2 = 24.36;
        del_n = n_gas-1;
        eta_crt = (1+N_0 * lam.^2 *r_e./(2*pi*del_n)).^-1;
        [b,Ions,Ephotmax,Gen_Int]=plot_ADK(handles,E1,E2,eta_crt);
        [E_cut_ef,nmbr_ems,E_phot_up] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt);
        targ = str2double(get(handles.target_photon,'String'));
        if E_cut_ef>=targ
            yes_no = 1;
        else
            yes_no = 0;
        end
        data = struct('yn',yes_no,'ems',nmbr_ems,'Ion',Ions,'Top_phot',Ephotmax,'E_phot_up',E_phot_up,'T_enrgy',enrgy,'alpha_gas',alpha_300,'Gen_Int',Gen_Int);
        hObject.UserData = data;
    case 5 % Xenon
        abs_data = load('Xe_300mb_1mm.csv');
        enrgy = abs_data(:,1);
        T = abs_data(:,2);
        alpha_300 = -log(T)/1E-3;
        n_gas = 1+0.012055*(0.26783*(lam*1E6)^2/(46.301*(lam*1E6)^2-1)+0.29841*(lam*1E6)^2/(50.578*(lam*1E6)^2-1)+5.0333*(lam*1E6)^2/(112.74*(lam*1E6)^2-1));
        E1 = 12.13;
        E2 = 21.21;
        del_n = n_gas-1;
        eta_crt = (1+N_0 * lam.^2 *r_e./(2*pi*del_n)).^-1;
        [b,Ions,Ephotmax,Gen_Int]=plot_ADK(handles,E1,E2,eta_crt);
        [E_cut_ef,nmbr_ems,E_phot_up] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt);
        targ = str2double(get(handles.target_photon,'String'));
        if E_cut_ef>=targ
            yes_no = 1;
        else
            yes_no = 0;
        end
        data = struct('yn',yes_no,'ems',nmbr_ems,'Ion',Ions,'Top_phot',Ephotmax,'E_phot_up',E_phot_up,'T_enrgy',enrgy,'alpha_gas',alpha_300,'Gen_Int',Gen_Int);
        hObject.UserData = data;
    case 6 % Hydrogen
        abs_data = load('H2_300mb_1mm.csv');
        enrgy = abs_data(:,1);
        T = abs_data(:,2);
        alpha_300 = -log(T)/1E-3;
        n_gas = 1+0.012055*(0.26783*(lam*1E6)^2/(46.301*(lam*1E6)^2-1)+0.29841*(lam*1E6)^2/(50.578*(lam*1E6)^2-1)+5.0333*(lam*1E6)^2/(112.74*(lam*1E6)^2-1));
        E1 = 13.6;
        E2 = 13.6;
        del_n = n_gas-1;
        eta_crt = (1+N_0 * lam.^2 *r_e./(2*pi*del_n)).^-1;
        [b,Ions,Ephotmax,Gen_Int]=plot_ADK(handles,E1,E2,eta_crt);
        [E_cut_ef,nmbr_ems,E_phot_up] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt);
        targ = str2double(get(handles.target_photon,'String'));
        if E_cut_ef>=targ
            yes_no = 1;
        else
            yes_no = 0;
        end
        data = struct('yn',yes_no,'ems',nmbr_ems,'Ion',Ions,'Top_phot',Ephotmax,'E_phot_up',E_phot_up,'T_enrgy',enrgy,'alpha_gas',alpha_300,'Gen_Int',Gen_Int);
        hObject.UserData = data;
    
    otherwise    
end
% Hints: contents = cellstr(get(hObject,'String')) returns Gas_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Gas_type


% --- Executes during object creation, after setting all properties.
function Gas_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gas_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pulse_energy_Callback(hObject, eventdata, handles)
% hObject    handle to pulse_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulse_energy as text
%        str2double(get(hObject,'String')) returns contents of pulse_energy as a double


% --- Executes during object creation, after setting all properties.
function pulse_energy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulse_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pulse_duration_Callback(hObject, eventdata, handles)
% hObject    handle to pulse_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulse_duration as text
%        str2double(get(hObject,'String')) returns contents of pulse_duration as a double


% --- Executes during object creation, after setting all properties.
function pulse_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulse_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output_power_Callback(hObject, eventdata, handles)
% hObject    handle to output_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','Gas_type');
data = h.UserData;
yes_no = data.yn;
if yes_no ==1;
    Out = 'Yes';
else
    Out = 'No - not generated';
end
set(hObject,'String',Out)
% Hints: get(hObject,'String') returns contents of output_power as text
%        str2double(get(hObject,'String')) returns contents of output_power as a double


% --- Executes during object creation, after setting all properties.
function output_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Laser_wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to Laser_wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Laser_wavelength as text
%        str2double(get(hObject,'String')) returns contents of Laser_wavelength as a double


% --- Executes during object creation, after setting all properties.
function Laser_wavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Laser_wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spot_size_Callback(hObject, eventdata, handles)
% hObject    handle to spot_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spot_size as text
%        str2double(get(hObject,'String')) returns contents of spot_size as a double


% --- Executes during object creation, after setting all properties.
function spot_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spot_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [A,B,C,D] = plot_ADK(handles,Eion,Eion2,eta_crit)
lam = str2double(get(handles.Laser_wavelength,'String'));
lam = lam*1E-6;
lam_au  = lam/5.29E-11; % wavelength in au
w_L = 2*pi*137/lam_au;
Duration = str2double(get(handles.pulse_duration,'String'));
Duration = Duration*1E-15;
Energy = str2double(get(handles.pulse_energy,'String'));
Energy = Energy*1E-3;
w0 = str2double(get(handles.spot_size,'String'));
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
set(0,'DefaultAxesFontSize',22)
axes(handles.axes1)
axis tight
[AX1 ,H1 ,H2] = plotyy(1E15*T/(4.1322E16/1.66),Ions,1E15*T/(4.1322E16/1.66),Ephotmax);
set(H1,'color',[0.8 0.0 0.1])
set(H2,'color',[0.4 0.0 0.7])
set(AX1(1),'ycolor',[0.8 0.0 0.1],'linewidth',4,'xlim',1E15*[-Duration*1.2 Duration*1.2],'ylim',[0 1],'ytick',[0:0.25:1],'xtick',1E15*1.5*linspace(-Duration,Duration,5),'yminortick','on','box','off')
set(AX1(2),'ycolor',[0.4 0.0 0.7],'linewidth',4,'xlim',1E15*[-Duration*1.2 Duration*1.2],'ylim',[0 1.1*max(Ephotmax)],'ytick',round(linspace(0,1.2*max(Ephotmax),4)),'xtick',1E15*1.5*linspace(-Duration,Duration,5),'yminortick','on','box','off')
set(AX1(2), 'XTickLabel','','XAxisLocation','Top')
hold(AX1(2),'on')
if eta_crit<max(Ions)
    [~,b]=min(abs(Ions-eta_crit));
else
    [~,b]=max(Ephotmax);
end

if eta_crit<max(Ions)
    plot(AX1(2),1E15*T(b)*ones(1,20)/(4.1322E16/1.66),linspace(0,max(Ephotmax)*1.2,20),'.','markersize',12,'color','k')
    hold(AX1(2),'off')
else
    hold(AX1(2),'off')
end    
xlabel('Time (fs)')
ylabel(AX1(1),'Ionization Fraction')
ylabel(AX1(2),'Cut-off energy (eV)')

D = Ain0*ft(b).^2; %generating intensity
A=b;
B=Ions;
C=Ephotmax;



function [E_cut_eff, nmbr_emttrs,Ephot_up] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt)
        w0 = str2double(get(handles.spot_size,'String'));
        w0 = w0*1E-6;
        lam = str2double(get(handles.Laser_wavelength,'String'));
        lam = lam*1E-6;
        E_cut_eff = Ephotmax(b);
        eta_up = Ions(1:b-1);
        Ephot_up = Ephotmax(1:b-1);
        Phs_mtch_prsr = 1013*lam^2./(2*pi^2*w0.^2*del_n*(1-eta_up/eta_crt));
        V_foc = pi^2 *w0.^4/(2*lam);
        nmbr_dnsty = 100*(Phs_mtch_prsr)/(1.38E-23*294); % 100 factor to convert to Pa
        nmbr_emttrs = V_foc .*nmbr_dnsty;
        axes(handles.axes2)
        set(0,'DefaultAxesFontSize',22)
        scale = 1E-13;
        [AX,H1,H2] = plotyy(Ephot_up,nmbr_emttrs*scale,Ephot_up,Phs_mtch_prsr);
        set(H1,'color','k')
        set(H2,'color','k')
        set(AX(2),'ycolor','k','linewidth',4,'xlim',[min(Ephot_up) 1.1*max(Ephot_up)],'ylim',[min(Phs_mtch_prsr)*0.8 max(Phs_mtch_prsr)*1.1],'xtick',linspace(0,1.1*max(Ephot_up),5),'xminortick','on','yminortick','on','box','off')
        %set(AX(1),'ycolor','k','linewidth',4,'xlim',[min(Ephot_up) 1.1*max(Ephot_up)],'ylim',[min(nmbr_emttrs)*0.8 1.1*max(nmbr_emttrs)]*scale,'xtick',round(linspace(0,1.1*max(Ephot_up),5)),'ytick',round(scale*linspace(min(nmbr_emttrs)*0.8,1.1*max(nmbr_emttrs),5),0),'xminortick','on','yminortick','on','box','off')
        set(AX(1),'ycolor','k','linewidth',4,'xlim',[min(Ephot_up) 1.1*max(Ephot_up)],'ylim',[min(nmbr_emttrs)*0.8 1.1*max(nmbr_emttrs)]*scale,'xtick',round(linspace(0,1.1*max(Ephot_up),5)),'xminortick','on','yminortick','on','box','off')
        set(AX(2), 'XTickLabel','','XAxisLocation','Top')
        ylabel(AX(2),'P_m (mbar)')
        ylabel(AX(1),'# emitters (\times10^{13})')
        xlabel('Cut-off photon energy (eV)')

        
function target_photon_Callback(hObject, eventdata, handles)
% hObject    handle to target_photon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of target_photon as text
%        str2double(get(hObject,'String')) returns contents of target_photon as a double


% --- Executes during object creation, after setting all properties.
function target_photon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_photon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%data = struct('yn',yes_no,'ems',nmbr_ems,'Ion',Ions,'Top_phot',Ephotmax,'E_phot_up',E_phot_up,'T_enrgy',enrgy,'alpha_gas',alpha_300);        
h = findobj('Tag','Gas_type');
data = h.UserData;
lam = str2double(get(handles.Laser_wavelength,'String'));
targ = str2double(get(handles.target_photon,'String'));
w0 = str2double(get(handles.spot_size,'String'));
w0 = w0*1E-6;
Ephotmax = data.Top_phot;
eta_q = data.Ion;
em = data.ems;
T_en = data.T_enrgy;
alpha_300 = data.alpha_gas;
E_phot_X = data.E_phot_up;
GI = data.Gen_Int;
[~,c1] = min(abs(Ephotmax(1:length(Ephotmax)/2)-targ));
[~,c2] = min(abs(E_phot_X-targ));
[~,c3] = min(abs(T_en-targ));
num_ems = em(c2);
V_foc = pi^2*w0.^4/(2*lam*1E-6); % Update here and throughout
zr = 2*V_foc/(pi*w0.^2);
Pm = num_ems*1.38E-23*294/(V_foc*100);
alpha = alpha_300(c3)*Pm/300;
mean_T = (2/(alpha*zr)).*(1-exp(-alpha*zr/2));

% figure of merit for harmonic efficiency
FoM = (GI/1E14)^3*((0.8/lam).^9)*(em(c2)*1E-14*(eta_q(c1))).^2./w0^2*mean_T; % for rationale behind w0^-2 factor see Heyl or Midorikawa

set(hObject,'String',FoM)


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [n_gas] = Sellmeier(lam,cfs)
if length(cfs)<7
    cfs(length(cfs)+1:7) = 0;
end
n_gas = 1+cfs(1)*(cfs(2)*lam.^2./(cfs(3)*lam.^2-1)+...
        cfs(4)*lam.^2./(cfs(4)*lam.^2-1)+cfs(6)*lam.^2./(cfs(7)*lam.^2-1));

    %Sells = load('Sellmeir_coeffs_nobles.csv'); % top line should be names
    
    
%     function[data] = Main_calc(abs_data,handles,E1,E2)       
%         
%         enrgy = abs_data(:,1);
%         T = abs_data(:,2);
%         alpha_300 = -log(T)/1E-3;
%         cfs = [1E-6 
%         n_gas = 1+ 1E-6*(67.86711+30182.943*(lam*1E6).^2./(144*(lam*1E6).^2-1));
%         E1 = 15.76;
%         E2 = 27.65;
%         del_n = n_gas-1;
%         eta_crt = (1+N_0 * lam.^2 *r_e./(2*pi*del_n)).^-1;
%         [b,Ions,Ephotmax,Gen_Int] = plot_ADK(handles,E1,E2,eta_crt);
%         [E_cut_ef,nmbr_ems, E_phot_up] = plot_phs_mtch_prsr(handles,b,Ions,Ephotmax,del_n,eta_crt);
%         %targ = str2double(get(handles.target_photon,'String'));
%         if E_cut_ef>=targ
%             yes_no = 1;
%         else
%             yes_no = 0;
%         end
%         data = struct('yn',yes_no,'ems',nmbr_ems,'Ion',Ions,...
%             'Top_phot',Ephotmax,'E_phot_up',E_phot_up,'T_enrgy',enrgy,...
%             'alpha_gas',alpha_300,'Gen_Int',Gen_Int);
%         hObject.UserData = data;
