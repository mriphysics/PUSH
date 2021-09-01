%%% [Mss,Mt] = SPGR_BSB(b1pulse,dt,Delta_Hz,TR,pulse_sep,tissuepars,Nt)
%
%   Steady-state MT SPGR sequence with eigenvector based time integration method
%
%   INPUTS:         
%           b1pulse    = RF pulse, Mx1 array (M=#timepoints). Units are uT
%           dt         = dwell time, sec
%           Delta_Hz   = frequency of the MT sat pulse. Units Hz
%           TR         = repetition time, sec
%           pulse_seq  = time interval between saturation and readout
%                        pulses, sec
%           tissuepars = structure containing all tissue parameters. See
%                        init_tissue()
%           Nt         = number of pulses to simulate
%
%  OUTPUTS:
%           Mss        = Steady-state Mxy (after excitation pulse)
%           Mt         = Transient response for Nt pulses starting at
%                        equilibrium
%
% (c) Shaihan Malik 2019. King's College London
%     Adapted by David Leitao 2020 to single MT saturation pulse

function [Mss,Mt] = SPGR_BSB(b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars, Nt, varargin)

% unpack tissue parameters
M0f = tissuepars.M0 * (1-tissuepars.f); 
M0s = tissuepars.M0 * tissuepars.f;     

R1  = [tissuepars.free.R1 tissuepars.semi.R1];
R2f = tissuepars.free.R2;

% semisolid T2, used in lineshape calculation
T2s = tissuepars.semi.T2;

% overall exchange rate for free and both semisolid pools
K = tissuepars.K;

% lineshape at Delta_Hz and 0
switch tissuepars.lineshape
    case 'SL'
        [G,~] = SuperLorentzian_lineshape_7T(T2s,[0 Delta_Hz]); % seconds
end

% gamma for RF calculation
gam = 267.5221; %< rad /s /uT


%% Lambda matrix is time invariant

% no off-resonance here
Lambda = [-R2f      0       0               0               0;
           0       -R2f     0               0               0;
           0        0      -K*M0s-R1(1)     K*M0f           R1(1)*M0f;
           0        0       K*M0s          -K*M0f-R1(2)     R1(2)*M0s;
           0        0       0               0               0];

% evolution in time between pulses
Nt_exc  = size(b1pulse.exc,1); 
tau_exc = dt*Nt_exc;
Nt_sat  = size(b1pulse.sat,1); 
tau_sat = dt*Nt_sat;
% pulse_sep supplied is interval between extremeties of two pulses; 
% in hard pulse approximation add half duration of each pulse to calculate
% interval between the centre of the pulses
pulse_sep = pulse_sep + tau_sat/2 + tau_exc/2; 

Xtilde_norf           = expm(Lambda*(TR - pulse_sep));    
Xtilde_norf_pulse_sep = expm(Lambda*pulse_sep);


%% RF matrices and integrate over mean (1) excitation and (2) saturation pulses

%(1) excitation pulse
b1x = sum(real(b1pulse.exc))*dt;
b1y = sum(imag(b1pulse.exc))*dt;
% calculate saturation
w1 = gam * sqrt(sum(abs(b1pulse.exc).^2)*dt);
w  = pi*w1^2*G(1); 

Omega = [0          0          -gam*b1y    0   0;
         0          0           gam*b1x    0   0;
         gam*b1y   -gam*b1x     0          0   0;
         0          0           0         -w   0;
         0          0           0          0   0];
     
% instantaneous RF pulse
Xtilde_rf_exc = expm(Omega);


%(2) saturation pulse
b1x = sum(real(b1pulse.sat))*dt;
b1y = sum(imag(b1pulse.sat))*dt;
% calculate saturation
w1 = gam*sqrt(sum(abs(b1pulse.sat).^2)*dt);
w  = pi*w1^2*G(2);

Omega = [0          0          -gam*b1y    0   0;
         0          0           gam*b1x    0   0;
         gam*b1y   -gam*b1x     0          0   0;
         0          0           0         -w   0;
         0          0           0          0   0];
    
% instantaneous RF pulse
Xtilde_rf_sat = expm(Omega);   
    

%% Now compile these into a sequence prediction

% perfect spoiling; apply at end of TR immediately before next RF pulse
Phi = diag([0 0 1 1 1]);

% combine all operators in the following order:
% 1. evolution during period TR-pulse_sep
% 2. spoiling
% 3. saturation pulse
% 4. evolution during period pulse_sep 
% 5. spoiling
% 6. excitation pulse
% the steady state is calculated in the time point right after the 
% excitation pulse 
X = Xtilde_rf_exc * Phi * Xtilde_norf_pulse_sep * Xtilde_rf_sat * Phi * Xtilde_norf;

% eigenvector decomposition
[v,~,~] = eig(X);
Mss = v(1:end-1,end)/v(end,end); % normalise to the last component so it stays as 1 (rounding errors from eig function)

% simulate dynamics of transient state for Nt TRs (starting at equilibrium)
if nargout>1
    Mt = zeros(5, Nt+1);
    Mt(:,1) = [0; 0; M0f; M0s; 1];
    for nn=1:Nt
        Mt(:,nn+1) = X * Mt(:,nn);
    end
    Mt(5,:) = [];
end

end