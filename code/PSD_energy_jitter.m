function PSD_energy_jitter(filename)
%% Loading File
load(filename);
%% Measurement and Projections;

h_measure = [a.NTout.Uout;conj(a.NTout.Uout)];  % UNIT: SI, For energy jitter

VR = [a.ev(1:a.N,:)+1i*a.ev((1+a.N):end,:); a.ev(1:a.N,:)-1i*a.ev((1+a.N):end,:)]*sqrt(a.N);  % Say, unit: 1e6
VL = conj([V_left(1:a.N,:)-1i*V_left((1+a.N):end,:); ...
    V_left(1:a.N,:)+1i*V_left((1+a.N):end,:)]/2)/sqrt(a.N);  % So, unit: 1e6

h_vec = conj(VR'*h_measure*a.dt);  % unit: 1e6*1e-12 = 1e-6 SI

Lambda = a.ew;

if 1
    ind_big = find(abs(h_vec)/max(abs(h_vec)) > 1e-4);
    ind = max(ind_big):length(Lambda);   
%     ind = 512:length(Lambda);   
    VL(:,ind) = [];
    Lambda(ind) = [];
    h_vec(ind) = [];
end

%% Dynamical Method
wavelength0 = 1564e-9;
nu0 = 3e8/wavelength0;  % unit: SI
h_planck = 6.62607004e-34; % SI
TR = 1/300e6;
S_auto = a.NTout.gsat*h_planck*nu0*TR;  % unit: Watt * s
% S_auto = 1e10;

Hm = (h_vec*h_vec'); % unit: 1e-12 SI
Dm = (VL'*VL)*a.dt*S_auto;  %  unit: 1e12*1e-12 = SI

Nm = length(Lambda);
Lm = Lambda*ones(1,Nm);
Cross_Add_Lambda = Lm + Lm';
Iden_zero = (abs(Cross_Add_Lambda)>1e-8);
M = Iden_zero.*Hm.*Dm;
Sh_f = @(f) sum(sum(M./((Lm-2i*pi*f).*(Lm'+2i*pi*f))));

faxis  = [0:0.01:100]*1e6;
f_comp = faxis*TR;

parfor ifreq = 1:length(faxis)
    PSD_h(ifreq) = Sh_f(f_comp(ifreq));
end
% toc
PSD_plot = abs(PSD_h*1e-12); % unit: SI

figure(1); hold on;
w0 = norm(a.NTout.Uout)^2*a.dt*1e-12;
plot((faxis),10*log10(2*PSD_plot/w0^2),'--');

%% Haus Mecozzi Model

if 1
    Parameters_SESAM300_paper2;
    wavelength0 = 1564e-9; % unit: SI
    nu0 = 3e8/wavelength0;  % unit: SI
    %TR = 1/300e6; % unit: SI
    h_planck = 6.62607004e-34; % SI
    S_auto = a.NTout.gsat*h_planck*nu0*TR;  % unit: Watt * s
    
    A0 = norm(a.NTout.Uout,inf); % unit: SI
    w0 = norm(a.NTout.Uout)^2*a.dt*1e-12;  % unit: SI
    tau0 = PlsWidth_FWHM(a.t,a.NTout.Uout)/1.763*1e-12;   % unit: SI
    Omega = Omega*1e12;
    
    var_w = 2*w0*S_auto;    
    
    %Decay Coefficients
    del = -gam*gain_sat/2/Omega^2/beta2;
    g1 = a.NTout.gsat^2/(g0*PsatTR*1e-12);  % unit: SI
    rw = -(2*del*A0^2 - g1*w0 + 2*g1*A0^2/(6*Omega^2*tau0));
    Swf_analytical = var_w./(rw^2+(2*pi*f_comp).^2);
    
    figure(1); 
    plot(faxis,10*log10(2*Swf_analytical/w0^2),'r-');
    hold on;
    
end
%% Figure Customization
figure(1);
legend('Dynamical', 'Haus-Mecozzi', 'location', 'northeast');
xlabel('Frequency (MHz)');
ylabel('S_{\omega}(f) (dBc/Hz)');
fig = gca;
fig.XRuler.Exponent = 6;
title('Energy Jitter Power Spectral Density');
xlim([0,6]*1e7);
grid on;
end