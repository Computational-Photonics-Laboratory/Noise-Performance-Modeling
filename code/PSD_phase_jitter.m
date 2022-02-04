function PSD_phase_jitter(filename)
%% Loading passed in File
load(filename);

%% Measurement and Projections;

w0 = norm(a.NTout.Uout)^2*a.dt*1e-12;  % unit: SI
h_measure = [a.t.*a.NTout.Uout;conj(a.t.*a.NTout.Uout)]*1e-12/w0;  % UNIT: SI, For timing jitter
% h_measure = [a.t.*a.NTout.Uout;conj(a.t.*a.NTout.Uout)]/w0/100;  % UNIT: SI, For timing jitter

VR = [a.ev(1:a.N,:)+1i*a.ev((1+a.N):end,:); a.ev(1:a.N,:)-1i*a.ev((1+a.N):end,:)];  % Say, unit: 1e6
VL = conj([V_left(1:a.N,:)-1i*V_left((1+a.N):end,:); V_left(1:a.N,:)+1i*V_left((1+a.N):end,:)]/2);   % So, unit: 1e6
Lambda = a.ew;
h_vec = conj(VR'*h_measure*a.dt);  % unit: 1e6*1e-12 = 1e-6 

if 1
    ind_big = find(abs(h_vec)/max(abs(h_vec)) > 1e-5);
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
M = Iden_zero.*Hm.*Dm;  % unit: 1e-36
Sh_f = @(f) sum(sum(M./((Lm-2i*pi*f).*(Lm'+2i*pi*f))));

faxis  = [0:0.05:50]*1e6;
f_comp = faxis*TR;

parfor ifreq = 1:length(faxis)
    PSD_h(ifreq) = Sh_f(f_comp(ifreq))*1e-12;
end
% toc

%% Timing Jitter Calculation
RandWalk = sum(sum((1-Iden_zero).*Hm.*(Dm+Dm')))/2*1e-12;  % unit: SI
PSD_plot = abs((2*pi*f_comp).^2.*PSD_h + RandWalk); % unit: SI

figure(3); hold on;
PSD_timing = PSD_plot./(f_comp*TR).^2;
plot((faxis),10*log10(2*PSD_timing),'--');
axis tight;
xlim([0,35]*1e6);

%% Haus Mecozzi
if 1
    Parameters_SESAM300_paper2;
    wavelength0 = 1564e-9; % unit: SI
    nu0 = 3e8/wavelength0;  % unit: SI
    %TR = 1/300e6; % unit: SI
    h_planck = 6.62607004e-34; % SI
    S_auto = a.NTout.gsat*h_planck*nu0*TR;  % unit: Watt * s
%     S_auto = 1e10*S_auto;  % this is just an arbituary scale
    
    w0 = norm(a.NTout.Uout)^2*a.dt*1e-12;  % unit: SI
    tau0 = PlsWidth_FWHM(a.t,a.NTout.Uout)/1.763*1e-12;   % unit: SI
    Omega = Omega*1e12;
    beta2 = beta2*1e-24;
    
    var_f = 2/(3*w0*tau0^2)*S_auto;  % unit: SI
    var_t = (pi*tau0)^2/(6*w0)*S_auto; % unit: SI
    
    %% Decay Coefficients
    rf = 2*a.NTout.gsat/3/(Omega*sqrt(2)*tau0)^2; % unitless
    rt = beta2;  % unit: SI
    
    Stf_analytical = rt^2*var_f./(rf^2+(2*pi*f_comp).^2) + var_t;
    figure(3);
    hold on;
    Sphif_analytical = Stf_analytical./(TR*f_comp).^2;
    plot(faxis,10*log10(2*Sphif_analytical),'r-');
end
    figure(3);
    hold on;
    legend('Dynamical', 'Haus-Mecozzi', 'location', 'northeast');
    xlabel('Frequency (MHz)');
    ylabel('S_{\Psi}(f) (dBc/Hz)');
    title('Phase Jitter Power Spectral Density');
    fig = gca;
    fig.XRuler.Exponent = 6;
    grid on;
end