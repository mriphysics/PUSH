%%% function [b1_pulse] = gen_B1_pulse(theta, tau, dt, pulse_shape)
%
%   Calculates the B1 waveform for a specific pulse shape and flip angle
%   
%   INPUTS:         
%           theta       = nutation angle (radians)
%           tau         = duration of the RF pulse (s)
%           dt          = time discretization step (s)
%           pulse_shape = string with pulse shape: 'Rect' or 'Gaussian'
%
%   OUTPUTS:
%           b1_pulse    = B1 waveform of the pulse (uT)
%
% (c) David Leitao 2019. King's College London

function [b1_pulse, teff] = gen_B1_pulse(theta, tau, dt, pulse_shape)

gam = 267.5221; %< rad /s /uT
nt = round(tau/dt);

switch pulse_shape
    case 'Rect'
        pulse = ones(nt, 1);
    case 'Gaussian'
        gausswin_alpha = 3; %<-- shape of basic pulse
        pulse = gausswin(nt, gausswin_alpha);
        pulse = pulse - min(pulse);
        pulse = pulse / max(pulse);
    otherwise
        error('Pulse shape not identified')
end

teff =  (gam*sum(pulse)*dt)/(gam*max(pulse)*tau); %<-- teff of shape - normalised to duration

b1_max = (theta)/(gam*teff*tau); % uT units

%%% Scale pulse to achieve desired nutation angle
b1_pulse = pulse * b1_max;

end

