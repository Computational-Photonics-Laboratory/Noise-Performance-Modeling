%% Set to 1 to generate new pulse file, 0 otherwise
%% Parameters
a = hme_eigen();

TS.mode = 'tw2tau';
TS.N = 512*2;
TS.tw2tau = 50*2;
TS.switch_FourierFD = 'off';

%% Object setup parameters
Parameters_SESAM300_paper2;

%%
T = 200; %Time window
N = 4096; %Num points in domain
% a.l = loss;
% a.g0 = g0;
% a.g0 = 7.3;
% a.b2 = beta2;
a.setup(T,N,g0,Omega,PsatTR,loss,beta2,gam,rho,TA,100);
a.g0 = 7.74;

%% Other System Parameters
w0 = 190;  %unit: pJ
tau0 = 0.133; %unit: ps

U0 = w0/2;   % pJ
A0 = sqrt(U0/tau0); %W^(1/2)

phi = 1/2*A0^2;
taus = 1;

us = A0*(sech(a.t/(tau0)));

[u1, phi1] = a.newton_solver(us,taus,phi,TS);

if a.NTout.res > 1e-5
    error('No Convergence!!!');
end

%% Eigen-analysis
switch_FourierFD_Eigen = 'on';
switch_LeftEigen = 'on';
[V_left] = a.CompleteEigen(switch_FourierFD_Eigen, switch_LeftEigen);    


%% Saves workspace
filename = sprintf('Pulse_g0_%.5g__b2_%.5g__loss_%.5g.mat', a.g0,-a.b2,a.l(1));
save(filename);

%% Calls scripts for creating PSD figures
PSD_energy_jitter(filename);
PSD_freq_jitter(filename);
PSD_phase_jitter(filename);
