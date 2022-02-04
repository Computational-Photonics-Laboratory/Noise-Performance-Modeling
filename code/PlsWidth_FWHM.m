function tau = PlsWidth_FWHM(t,u)

u = (abs(u)/norm(u,inf)).^2;
[~,ind_max] = max(u);
tl = spline(u(1:ind_max),t(1:ind_max),1/2);
tr = spline(u(ind_max:end),t(ind_max:end),1/2);
tau = abs(tr-tl);
end