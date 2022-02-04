Omega = 30;
PsatTR = 30;   % W*ps = pJ
g0 = 7.74;
loss = 1.05;
beta2 = -0.0144;
gam = 1.11*1e-3; % W^-1
rho = 0.0726;

TR = 1/300e6;

% the saturable absorber
TA = 2; % ps
wA = 157; % 157 pJ

wavelength0 = 1560e-9; % m

gain_sat = loss+rho-0.025;  % just an estimation

w0 = 1.94;  % unit: pJ
tau0 = 0.139;
phi0 = 0.8;
A0 = 26.5;


