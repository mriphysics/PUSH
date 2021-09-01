function pars = init_tissue(name)
%pars = init_tissue(name)
% Function to generate structure containing tissue properties. Possible
% name inputs are:
%   - 'wm_7t'  = white matter tissue parameters at 7T
%
%
%  Shaihan Malik (Feb 2019)
%  David Leitao  (Jan 2021)

pars = struct;

% Replace all letters with lower case
name = lower(name);

% initialise the empty structure
% Rates are all s^-1, times on seconds

%%% Free water pool
pars.free.R1 = [];
pars.free.R2 = [];
pars.free.M0 = [];  
%%% Semisolid pool
pars.semi.M0 = [];   
pars.semi.R1 = [];   
pars.semi.T2 = [];  
%%% General parameters
pars.K  = [];
pars.f  = [];
pars.M0 = 1; % by definition sum of all magnetizations equals 1
pars.lineshape = [];

switch name
    
    case 'wm_7t'
  
        pars.M0 = 1;
        pars.f  = 0.1357;   
        pars.K  = 32.79;    
        pars.lineshape = 'SL'; 
        
        pars.free.R1 = 0.4;
        pars.free.R2 = 1/81e-3; 
        pars.free.M0 = pars.M0 * (1-pars.f); 
        
        pars.semi.M0 = pars.M0 * pars.f;    
        pars.semi.R1 = 1.85;
        pars.semi.T2 = 9.6e-6; 

end

end