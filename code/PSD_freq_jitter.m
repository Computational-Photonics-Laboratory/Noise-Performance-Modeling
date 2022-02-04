function PSD_freq_jitter(filename)
%% Loading passed in File
load(filename);

%% Measurement and Projections;

w0 = norm(a.NTout.Uout)^2*a.dt*1e-12;
u0t = ifft(-1i*a.w.*fft(a.NTout.Uout));
h_measure = [u0t;-conj(u0t)]/(-1i*w0)*1e12;  % UNIT: SI, For freq jitter

VR = [a.ev(1:a.N,:)+1i*a.ev((1+a.N):end,:); a.ev(1:a.N,:)-1i*a.ev((1+a.N):end,:)];  % Say, unit: 1e6
VL = conj([V_left(1:a.N,:)-1i*V_left((1+a.N):end,:); V_left(1:a.N,:)+1i*V_left((1+a.N):end,:)]/2);   % So, unit: 1e6
h_vec = conj(VR'*h_measure*a.dt);  % unit: 1e6*1e-12 = 1e-6 SI
h_vec(7:end) = 0;

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

Hm = (h_vec*h_vec'); % unit: 1e-12 SI
Dm = (VL'*VL)*a.dt*S_auto;  %  unit: 1e12*1e-12 = SI
Nm = length(Lambda);
Lm = Lambda*ones(1,Nm);

Cross_Add_Lambda = Lm + Lm';
Iden_zero = (abs(Cross_Add_Lambda)>1e-8);

Term1 = 2*Iden_zero.*Hm;
Term2 = @(f) real(Dm./((Lm+Lm').*(Lm'+2i*pi*f)));
Sh_f = @(f) sum(sum(Term1.*Term2(f))); % unit: 1e-12 SI

faxis  = [0:0.01:50]*1e6;
f_comp = faxis*TR;

parfor ifreq = 1:length(faxis)
    PSD_h(ifreq) = Sh_f(f_comp(ifreq))*1e-12;  % unit: SI
end 
% toc
PSD_plot = abs(PSD_h)/nu0^2; % unit: SI

figure(2); hold on;
plot((faxis),10*log10(2*PSD_plot),'b--');


%% Haus Mecozzi Model
if 1
    Parameters_SESAM300_paper2;
    A0 = norm(a.NTout.Uout,inf); % unit: SI
    w0 = norm(a.NTout.Uout)^2*a.dt*1e-12;  % unit: SI
    tau0 = PlsWidth_FWHM(a.t,a.NTout.Uout)/1.763*1e-12;   % unit: SI
    Omega = Omega*1e12;
    var_p = 2/(3*w0*tau0^2)*S_auto;  % unit: SI

    %% Decay Coefficients
    rp = -2*a.NTout.gsat/3/(Omega*sqrt(2)*tau0)^2; % unitless
    
    Sff_analytical = var_p./(rp^2+(2*pi*f_comp).^2);
    
    figure(2);
    plot(faxis,10*log10(2*Sff_analytical/nu0^2),'r-');
    hold on;
end

%% Figure Customization
figure(2);
xlabel('Frequency (MHz)');
ylabel('S_{f_c}(f) (dBc/Hz)');
title('Frequency Jitter Power Spectral Density');
legend('Dynamical', 'Haus-Mecozzi', 'location', 'northeast');
fig = gca;
fig.XRuler.Exponent = 6;
grid on;
end