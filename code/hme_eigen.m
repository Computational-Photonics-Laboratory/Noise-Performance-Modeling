classdef hme_eigen < hme_newton_3N
    
    properties        
        ev
        ew        
    end
    
    methods
        function obj=hme_eigen()
        end
        
        function [V_left] = CompleteEigen(self,switch_FourierFD,switch_LeftEigen)
            
            u0 = self.NTout.Uout;
            x0 = [real(u0);imag(u0)];
            
            % %             x00 = [x0;self.NTout.n;self.NTout.gsat;self.NTout.taus;self.NTout.phi];
            % %             F00 = self.HMEsetup(x00);
            %             [J0,F0] = self.Jacobian('on',x0);
            %             figure(102); hold on;
            %             eps = 2.^(-10:-7);
            %             Dx = [rand*exp(-(self.t).^2);rand*exp(-(self.t).^2)];
            %             for ie = 1:length(eps)
            %                 epsilon = eps(ie);
            %                 x = x0 + epsilon*Dx;
            %                 [~,fval] = self.Jacobian('on',x);
            %                 err = log2(abs(fval - epsilon*J0*Dx));
            %                 plot(err); hold on; axis tight; grid on;
            %             end
            
            J0 = self.Jacobian(switch_FourierFD,x0);
            [X, D] = eig(full(J0));
            D = diag(D);
            [self.ev, D] = self.sort_EigenOutput_complex(X,D);
            self.ew = D;
            
            V_left = nan;
            if strcmp(switch_LeftEigen,'on')
                [X, D] = eig(full(J0).');
                % the complex conjugacy of D is due to the difference
                % between physical splitting and computational realization
                % when splitting the system:
                % Real and Imag --- OR --- u and u-bar
                % Left and Right eigenvectors --- OR --- Adjoint eigenvectors
                % Inner product without --- OR --- with complex conjugacy
                
                [V_left,D] = self.sort_EigenOutput_complex(X,conj(diag(D)));
                if abs(D(1)) < 1e-8 && abs(D(1)) < 1e-8
                    % First, regulate the eigenvectors of zero
                    % eigenvalues
%                     ut = self.Df1*self.NTout.Uout;
%                     self.ev(:,1) = [real(ut);imag(ut)];
%                     self.ev(:,2) = [real(1i*self.NTout.Uout);imag(1i*self.NTout.Uout)];
                    % Next, find the left eigenvectors of zero eigenvalues
                    
                    v1 = V_left(:,1);
                    v2 = V_left(:,2);
                    c1 = inner_product(v1,self.ev(:,1),self.dt);
                    c2 = inner_product(v1,self.ev(:,2),self.dt);
                    d1 = inner_product(v2,self.ev(:,1),self.dt);
                    d2 = inner_product(v2,self.ev(:,2),self.dt);
                    
                    V_left(:,1) = (v1*d2-v2*c2)/(c1*d2-c2*d1);
                    V_left(:,2) = (v1*d1-v2*c1)/(c2*d1-c1*d2);
                                       
%                     figure(335);
%                     crossproduct = V_left.'*self.ev;
%                     surf(abs(crossproduct(1:20,1:20))); colorbar; axis tight; view(2);
                end
                coeffs = (1./(self.dt*dot(conj(V_left),self.ev))).';
                V_left = V_left*diag(coeffs);
            end
        end
    end
    
           
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Some Utilities
    methods    
        function [J, F] = Jacobian (self,switch_FourierFD,x0) 
            if strcmp(switch_FourierFD,'on')
                [~, self.Df1] = fourdif(self.N,1);
                self.Df1 = self.dw*self.Df1;
                [~, self.Df2] = fourdif(self.N,2);
                self.Df2 = self.dw^2*self.Df2;
                self.dispD2 = self.b2*self.Df2;       % dispersion operator
            end
            
            vi = x0(1:self.N);
            wi = x0(1+self.N:end);
            Ei = norm(vi+1i*wi)^2;
            ti = self.NTout.taus;  % the "taus" in the equations
            gi = self.g0/(1 + Ei/self.Psat);
            g2 = -gi^2/(self.g0*self.Psat);
            phi = self.NTout.phi;    
            
            vi2 = vi.^2;
            wi2 = wi.^2;
            ui2 = vi2 + wi2;            
            viwi = vi.*wi;            
            w0 = cumtrapz(self.t,ui2); % the cumulative integration of the intensity of the pulse
            nf1 = exp(self.t/self.TA + w0/self.wA);
            nf1m = 1./nf1;
            nf2 = exp(self.t(1)/self.TA) + cumsimps(self.t,nf1)/self.TA;
            ni = nf1m.*nf2;                       
            
            I = speye(self.N);
            Mnf1 = self.diagsp(nf1);
            Mnf1m = self.diagsp(nf1m);
            Mnf2 = self.diagsp(nf2);
            commen_factor = 2/self.wA*Mnf1m*(self.Lm*Mnf1*self.Lm/self.TA-Mnf2*self.Lm);
       
            F1v = self.diagsp((gi-self.l-self.rho*ni)/2-2*self.gam*viwi)...
                - self.rho/2*self.diagsp(vi)*commen_factor*self.diagsp(vi)...
                + gi/(4*self.Omega^2)*self.Df2 + g2*(I + (1/2/self.Omega^2)*self.Df2)*vi*(vi)'...
                + ti*self.Df1;
            F1w = self.diagsp(phi/2-self.gam*(vi2+3*wi2))...
                - self.rho/2*self.diagsp(vi)*commen_factor*self.diagsp(wi)...
                + self.dispD2/2 ...
                + g2*(I + (1/2/self.Omega^2)*self.Df2)*vi*(wi)';
                                   
            F2v = self.diagsp(-phi/2+self.gam*(3*vi2+wi2))...
                - self.rho/2*self.diagsp(wi)*commen_factor*self.diagsp(vi)...
                - self.dispD2/2 ...
                + g2*(I + (1/2/self.Omega^2)*self.Df2)*wi*(vi)';
            F2w = self.diagsp((gi-self.l-self.rho*ni)/2+2*self.gam*viwi) ...
                - self.rho/2*self.diagsp(wi)*commen_factor*self.diagsp(wi)...
                + gi/(4*self.Omega^2)*self.Df2 + g2*(I + (1/2/self.Omega^2)*self.Df2)*wi*(wi)' ...
                + ti*self.Df1;
             
            J = [F1v, F1w; F2v, F2w]; 
            
            if nargout == 2
                vd1 = self.Df1_fun(vi);
                wd1 = self.Df1_fun(wi);
                vd2 = self.Df2_fun(vi);
                wd2 = self.Df2_fun(wi);
                
                F1 = (gi-self.l-self.rho*ni)/2.*vi + gi/(4*self.Omega^2)*vd2 + ti*vd1 + (phi*wi+self.b2*wd2)/2 - self.gam*wi.*ui2;
                F2 = -(phi*vi+self.b2*vd2)/2 + (gi-self.l-self.rho*ni)/2.*wi + gi/(4*self.Omega^2)*wd2 + ti*wd1 + self.gam*vi.*ui2;
                F = [F1;F2];
            end
            
        end
        
        function [J, F] = Jacobian1 (self,switch_FourierFD,x0) 
            if strcmp(switch_FourierFD,'on')
                [~, self.Df1] = fourdif(self.N,1);
                self.Df1 = self.dw*self.Df1;
                [~, self.Df2] = fourdif(self.N,2);
                self.Df2 = self.dw^2*self.Df2;
                self.dispD2 = self.b2*self.Df2;       % dispersion operator
            end
            
            vi = x0(1:self.N);
            wi = x0(1+self.N:end);
            Ei = norm(vi+1i*wi)^2;
            ti = self.NTout.taus;  % the "taus" in the equations
            gi = self.g0/(1 + Ei/self.Psat);
            g2 = -gi^2/(self.g0*self.Psat);
            phi = self.NTout.phi;    
            
            vi2 = vi.^2;
            wi2 = wi.^2;
            ui2 = vi2 + wi2;            
            viwi = vi.*wi;            
            w0 = cumtrapz(self.t,ui2); % the cumulative integration of the intensity of the pulse
            nf1 = exp(self.t/self.TA + w0/self.wA);
            nf1m = 1./nf1;
            nf2 = exp(self.t(1)/self.TA) + cumsimps(self.t,nf1)/self.TA;
            ni = nf1m.*nf2;                       
            
            I = speye(self.N);
            Mnf1 = self.diagsp(nf1);
            Mnf1m = self.diagsp(nf1m);
            Mnf2 = self.diagsp(nf2);
            commen_factor = 2/self.wA*Mnf1m*(self.Lm*Mnf1*self.Lm/self.TA-Mnf2*self.Lm);
       
            F1v = self.diagsp((gi-self.l-self.rho*ni)/2-2*self.gam*viwi)...
                - self.rho/2*self.diagsp(vi)*commen_factor*self.diagsp(vi);
            F1w = self.diagsp(phi/2-self.gam*(vi2+3*wi2))...
                - self.rho/2*self.diagsp(vi)*commen_factor*self.diagsp(wi);
                                   
            F2v = self.diagsp(-phi/2+self.gam*(3*vi2+wi2))...
                - self.rho/2*self.diagsp(wi)*commen_factor*self.diagsp(vi);
            F2w = self.diagsp((gi-self.l-self.rho*ni)/2+2*self.gam*viwi) ...
                - self.rho/2*self.diagsp(wi)*commen_factor*self.diagsp(wi);
             
            J = [F1v, F1w; F2v, F2w]; 
            
            if nargout == 2
                vd1 = self.Df1_fun(vi);
                wd1 = self.Df1_fun(wi);
                vd2 = self.Df2_fun(vi);
                wd2 = self.Df2_fun(wi);
                
                F1 = (gi-self.l-self.rho*ni)/2.*vi + gi/(4*self.Omega^2)*vd2 + ti*vd1 + (phi*wi+self.b2*wd2)/2 - self.gam*wi.*ui2;
                F2 = -(phi*vi+self.b2*vd2)/2 + (gi-self.l-self.rho*ni)/2.*wi + gi/(4*self.Omega^2)*wd2 + ti*wd1 + self.gam*vi.*ui2;
                F = [F1;F2];
            end
            
        end
                                  
    end
    
    methods (Static)
        function [V_out, D_out,orders] = sort_EigenOutput_complex(V,D)
            %% Orders:
            % 1. all real eigenvalues, descending
            % 2. all complex eigenvalues, descending by the real part
            
            [~,orders1] = sort(abs(imag(D)),'ascend');
            if length(D) > 4    
                order_discrete = orders1(1:4);
                order_others = orders1(5:end);

                D_tmp = D(order_discrete);          % the four real eigenvalues
                [~, order_d] = sort(real(D_tmp),'descend');  % sort the 4 real eigenvalues
                
                D_tmp = D(order_others);  % take out other eigenvalues 
                [~, order_o] = sort(real(D_tmp),'descend');  % sort the other eigenvalues
                
                orders = [order_discrete(order_d);order_others(order_o)];  % the sorting for 4 discrete eigenvalues is done!
                D_out = D(orders);
                V_out = V(:,orders);
%                 pause;
            else
                [D_out,order_d] = sort(real(D(orders1(1:4))),'descend');
                V_out = V(:,order_d);
            end
        end                 
    end
    
end
