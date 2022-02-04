classdef hme_parameters < handle
    
    properties        
         % 300 MHz       
        Omega = 30;
        PsatTR = 30;   % W*ps = pJ
        g0 = 7.74;
        l = 1.05;
        b2 = -0.0144;   %unit: ps^2
        gam = 1.11*1e-3; % W^-1
        rho = 0.0726;
        TA = 2;  % unit: ps
        wA = 100;  % unit: pJ
    end
    
    % Internal variables
    properties (GetAccess='public',SetAccess='public')
        N
        T
        dt
        dw
        naxis
        t
        w
    end
          
    %% setup function
    methods
        function obj=hme_parameters()            
        end
        
        function setup(self,N,T,delta,sgm)
            self.N = N; %Number of points in time domain
            self.T = T; %Computational time window
            self.delta=delta; 
            self.sgm=sgm;
            self.naxis = (-N/2:1:N/2-1)';
            
            self.dt = self.T/self.N;                   % temporal resolution
            self.t = self.naxis*self.dt;                         % time axis
            
            self.dw = 2*pi/self.T;                     % Freq resolution
            self.w = fftshift(self.naxis*self.dw);                         % Freq axis
        end
        
    end
end