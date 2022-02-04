classdef hme_newton_3N < hme_parameters
    properties
        Psat
        PsatF
    end
    
    % variables for eigenanalysis
    properties         
        options      % for solving the stationary solution
        NTout    % the output of newton solver
        TW2tau;  % the ratio of time window to pulse width
    end
    
    % Boundary Absorber
    properties (GetAccess='public',SetAccess='public')
        BAw
        BAh
        BAw2tau;
        BAf
    end
    % Internal matrices
    properties (GetAccess='public',SetAccess='public')        
        Df1
        Df1_fun
        Df2
        Df2_fun
        Lm
        Zv        
        dispD2
        diagsp
    end        
    
    methods
        function obj=hme_newton_3N()
        end
        
        function setup(self,T,N,g0,Omega,PsatTR,l,beta2,gam,rho,TA,wA)
            
            self.N = N;
            self.naxis = (-self.N/2:1:self.N/2-1).';
            self.T = T;
            self.dt = self.T/self.N;
            self.t = self.dt*self.naxis;
            self.dw = 2*pi/self.T;                     % Freq resolution
            self.w = fftshift(self.naxis*self.dw);                         % Freq axis
                                   
            if nargin > 3
                self.l = l;
                self.Omega = Omega;
                self.g0 = g0;
                self.PsatTR = PsatTR;
                self.b2 = beta2;
                self.gam = gam;
                self.rho = rho;
                self.TA = TA;
                self.wA = wA;
            end
            
            self.options = optimset('Display','iter','Jacobian','on','TolFun',1e-16,'TolX',...
                1e-10,'Algorithm','levenberg-marquardt','ScaleProblem','Jacobian',...
                'MaxIter',600,'LargeScale','on');
        end
        
        %% Solver
        
        function [u0, phi0] = newton_solver(self,u_in,taus_in,phi_in,TS)
            
            phi0 = phi_in;
            
%             figure(101);
%             plot(self.t,abs(u_in),'k.'); hold on;
            
            if 1
                u0 = self.normalization(u_in,TS);
            else
                u0 = u_in;
                phi0 = phi_in;
                self.Df1_fun = @(x) (ifft(1i*self.w.*fft(x)));
                self.Df2_fun = @(x) (ifft(-self.w.^2.*fft(x)));
                if strcmp(TS.switch_FourierFD,'on')
                    [~, self.Df1] = fourdif(self.N,1);
                    self.Df1 = self.dw*self.Df1;
                    [~, self.Df2] = fourdif(self.N,2);
                    self.Df2 = self.dw^2*self.Df2;
                else
                    self.Df1 = sparseFD(5,1,self.dt,self.N);
                    self.Df2 = sparseFD(7,2,self.dt,self.N);
                end
                
                self.dispD2 = self.b2*self.Df2;       % dispersion operator
                self.gainD2 = self.I + (1/2/self.Omega^2)*self.Df2;  %  Gain multiplier
                self.Psat = self.PsatTR/self.dt;       % Saturated power
            end
            
%             figure(101);
%             plot(self.t,abs(u0),'b--');
            
            gi0 = self.g0/(1 + norm(u0)^2/self.Psat);            
            n0 = self.absorber_cumtrapz(u0);
            taus0 = taus_in;
            x0 = [real(u0);imag(u0);n0;gi0;taus0;phi0];  
            
%             [F J] = self.HMEsetup(x0);
            
            f = @(x0) self.HMEsetup(x0);
            
%             figure(101); hold on;
%             eps = 2.^(-10:-7);
%             Dx = [exp(-(self.t).^2);exp(-(self.t).^2);exp(-(self.t).^2);rand(self.N*0+3,1)];
%             f0 = f(x0);
%             [~,J0] = self.HMEsetup(x0);
%             for ie = 1:length(eps)
%                 epsilon = eps(ie);
%                 x = x0 + epsilon*Dx;
%                 fval = f(x);
%                 err = log2(abs(fval - f0 - epsilon*J0*Dx));
%                 plot(err); hold on; axis tight; grid on;
%             end

            [xout,F,exitflag,output] = fsolve (f,x0,self.options);
            
            out.xout = xout;
            
            Uout = xout(1:self.N) + 1i*xout(self.N+1:2*self.N);            
            out.Uout = Uout; % pulse shape
            out.n = xout(self.N*2+1:self.N*3);
            out.gsat = xout(end-2);
            out.taus = xout(end-1);   % the time-shift per roundtrip;
            out.phi = xout(end);   % under normalization: the real 'phi' = out.phi*out.pu^4;                        
            out.res = norm(F); % Residual of the function [ ideally, F = \vec(0) ]
            out.exitflag = exitflag; % the status of output.
            out.t = self.t;
            out.T = self.T;
            out.output = output;
            self.NTout = out;
        end        
        
        function [F, J] = HMEsetup(self,x)
            vi = x(1:self.N);
            wi = x(self.N+1:self.N*2);
            ni = x(self.N*2+1:self.N*3);
            gi = x(self.N*3+1);
            ti = x(self.N*3+2);  % the "taus" in the equations
            phi = x(end);            
            
            Ei = norm(vi+1i*wi)^2;
           
            vd1 = self.Df1_fun(vi);
            wd1 = self.Df1_fun(wi);
            vd2 = self.Df2_fun(vi);
            wd2 = self.Df2_fun(wi);
            
            vi2 = vi.^2;
            wi2 = wi.^2;
            ui2 = vi2 + wi2;            
            viwi = vi.*wi;
            
            w0 = cumtrapz(self.t,ui2); % the cumulative integration of the intensity of the pulse
            nf1 = exp(self.t/self.TA + w0/self.wA);
            nf1m = 1./nf1;
            nf2 = exp(self.t(1)/self.TA) + cumsimps(self.t,nf1)/self.TA;     
            
            F1 = (gi-self.l-self.rho*ni)/2.*vi + gi/(4*self.Omega^2)*vd2 + ti*vd1 + (phi*wi+self.b2*wd2)/2 - self.gam*wi.*ui2;
            F2 = -(phi*vi+self.b2*vd2)/2 + (gi-self.l-self.rho*ni)/2.*wi + gi/(4*self.Omega^2)*wd2 + ti*wd1 ...
                + self.gam*vi.*ui2;
            F3 = nf1m.*nf2 - ni;
            F4 = gi*(self.Psat+Ei) - self.g0*self.Psat;
                                    
%             F = [F1];                       
            F = real([F1;F2;F3;F4]);                       
       
            F1v = self.diagsp((gi-self.l-self.rho*ni)/2-2*self.gam*viwi) + gi/(4*self.Omega^2)*self.Df2 + ti*self.Df1;
            F1w = self.diagsp(phi/2-self.gam*(vi2+3*wi2)) + self.dispD2/2;
            F1n = self.diagsp(-self.rho/2*vi);
            F1g = vi/2+1/(4*self.Omega^2)*vd2;
            F1tau = vd1;
            F1phi = wi/2;            
                       
            F2v = self.diagsp(-phi/2+self.gam*(3*vi2+wi2)) - self.dispD2/2;
            F2w = self.diagsp((gi-self.l-self.rho*ni)/2+2*self.gam*viwi) + gi/(4*self.Omega^2)*self.Df2 + ti*self.Df1;
            F2n = self.diagsp(-self.rho/2*wi);
            F2g = wi/2+1/(4*self.Omega^2)*wd2;
            F2tau = wd1;
            F2phi = -vi/2;
                                
            Mnf1 = self.diagsp(nf1);
            Mnf1m = self.diagsp(nf1m);
            Mnf2 = self.diagsp(nf2);

%             The original version, this is Correct!!!
%             commen_factor = 2/self.wA*Mnf1m*(self.Lm*Mnf1*self.Lm/self.TA-Mnf2*self.Lm);

            commen_factor = 2/self.wA*Mnf1m*(self.Lm*Mnf1/self.TA-Mnf2)*self.Lm;
            
            F3v = commen_factor*self.diagsp(vi);
            F3w = commen_factor*self.diagsp(wi);
            F3n = self.diagsp(-ones(self.N,1));
            F3g = self.Zv;
            F3tau = self.Zv;
            F3phi = self.Zv;
             
            F4v = 2*gi*vi.';
            F4w = 2*gi*wi.';
            F4n = self.Zv.';
            F4g = self.Psat + Ei;
            F4tau = 0;
            F4phi = 0;
            
             
            J = real([F1v, F1w, F1n, F1g, F1tau, F1phi; F2v, F2w, F2n, F2g, F2tau, F2phi;...
                F3v, F3w, F3n, F3g, F3tau, F3phi; F4v, F4w, F4n, F4g, F4tau, F4phi]); 
            
%             J = [zeros(self.N,self.N), F1w, zeros(self.N,self.N), zeros(self.N,1), zeros(self.N,1), zeros(self.N,1)];
%                 F2v, zeros(self.N,self.N), zeros(self.N,self.N), zeros(self.N,1), zeros(self.N,1), zeros(self.N,1)]; 
            
        end
    end % Methods
    
    %% utilities
    methods
        function tau = PlsWidth (self,u)
            % To compute the pulse width of "u"
            ua = abs(u).^2;
            E = norm(u)^2;
            tc = sum(self.t.*ua)/E;
            tau = 2*sqrt(sum((self.t-tc).^2.*ua)/E);
        end
        
        function out = centering(self,u,t,w) % to shift the center of pulse
            u = full(u);
            if (nargin<3)
                t = self.t;
                w = self.w;
            end            
            ua = abs(u).^2;
            E = norm(u)^2;
            tc = sum(t.*ua)/E;
            out = ifft(exp(1i*w*tc).*fft(u));
        end                
                        
        function [U0,TU_in,TU_out] = normalization(self,u0,TS)                       
            Tin.mode = 'tw'; Tin.tw = self.T;
%             Tout.mode = TS.mode; Tout.N = TS.N; Tout.tw2tau = TS.tw2tau;
            Tout = TS;
            [TU_in, TU_out] = Smp_Rate_Convert(u0,Tin,Tout,'on');            
            if 1
%                 tau0 = PlsWidth_FWHM(TU_out.t,TU_out.U)/1.763;
%                 TU_out.U = time_shift(TU_out.U,TU_out.t,(TU_out.t(1)+35*tau0));
                TU_out.U = time_shift(TU_out.U,TU_out.t,0);
            end
            IS.Uout = TU_out.U;
            IS.T = TU_out.tw;
            IS.t = TU_out.t;
            self.setup(IS.T,length(IS.t));
            U0 = TU_out.U;            
            
            self.Df1_fun = @(x) (ifft(1i*self.w.*fft(x)));
            self.Df2_fun = @(x) (ifft(-self.w.^2.*fft(x)));
            if strcmp(TS.switch_FourierFD,'on')
                [~, self.Df2] = fourdif(self.N,2);
                self.Df2 = self.dw^2*self.Df2;
                [~, self.Df1] = fourdif(self.N,1);
                self.Df1 = self.dw*self.Df1;
            else
                self.Df1 = sparseFD(5,1,self.dt,self.N);
                self.Df2 = sparseFD(7,2,self.dt,self.N);
            end
            
            BA = 0*fftshift(TU_out.U);
            self.l = self.l + BA;
                        
            self.dispD2 = self.b2*self.Df2;       % dispersion operator
            self.Psat = self.PsatTR/self.dt;       % Saturated power
            self.Lm = self.Trap_Quand_matrix(self.N)*self.dt;
            self.Zv = sparse(zeros(self.N,1));
            self.diagsp = @(x) spdiags(x,0,self.N,self.N);
        end
        
        function n = absorber_cumsimps(self,u0)   
            w = cumsimps(self.t,abs(u0).^2);
            exp_tw = exp(self.t/self.TA+w/self.wA);
            n = (exp(self.t(1)/self.TA) + cumsimps(self.t,exp_tw)/self.TA)./(exp_tw);
        end
        
        function q = absorber_cumtrapz(self,u0)            
            w = cumtrapz(self.t,abs(u0).^2);
            exp_tw = exp(self.t/self.TA+w/self.wA);
            q = (exp(self.t(1)/self.TA) + cumtrapz(self.t,exp_tw)/self.TA)./(exp_tw);               
        end
    end
    
    methods (Static)
        function Lm = Trap_Quand_matrix(N)
            Lm = tril(ones(N),-1)+eye(N)*1/2;
            Lm(:,1) = 1/2*ones(N,1);
            Lm(1,1) = 0;
            Lm = sparse(Lm);
        end
    end
    
end % Classdef
