%%% function [b1sqrd] = b1sqrd_integral(b1,TR,tau,dt)
%
%   Calculates the mean square B1+ for a given B1 waveform
%   
%   INPUTS:         
%           b1     = discretized B1 waveform (uT)
%           TR     = repetition time of the sequence (s)
%           tau    = duration of the RF pulse (s)
%           dt     = time discretization step (s)
%
%   OUTPUTS:
%           b1sqrd = mean square B1+ (uT^2)
%
% (c) David Leitao 2019. King's College London

function [b1sqrd] = b1sqrd_integral(b1,TR,tau,dt)

%Check if discretization of B1 matches the number of samples of the RF
if numel(b1)~=round(tau/dt)
    error('#B1 samples != #RF samples')
end

b1sqrd = (1/TR) * sum(dt*abs(b1).^2);

end

