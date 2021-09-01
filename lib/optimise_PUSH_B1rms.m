function [wopt, fval] = optimise_PUSH_B1rms(tx, target, Nsp, opt_options)
%[wopt, fval] = optimise_PUSH_B1rms(tx, target, Nsp, Sys, Sar, opt_options)
%   Optimises PUSH pulses to achieve target B1rms 
%
%           INPUTS:
%   tx          = transmit sensitivity maps for all coils (units of uT/V); 
%                 dimenions should be (Nr x Nch), where Nr are the number 
%                 of spatial positions and Nch is the number of channels
%   target      = target peak B1rms for all spatial positions (units of
%                 uT); dimensions should be (Nr x 1)
%   Nsp         = number of sub-pulses in the PUSH pulse
%   Sys         = Sys structure used in IDEA's Matlab framework
%   Sar         = Sar structure used in IDEA's Matlab framework
%   opt_options = structure containing optimisation parameters:
%       *.rng_seed          = seed for Matlab's pseudo random number
%                             generator (optional)
%       *.weights           = vector with spatial weights; dimensions
%                             should be (1 x Nr)
%       *.multistart        = {0,1} not/use multiple initialisations for
%                             the optimisation
%       *.multistart_trials = number of initialisations for the
%                             optimisation
%       *.verbose           = {'off','final','iter'} output from fmincon 
%       *.Niter             = maximum number of iterations for the
%                             optimisation
%       *.quality           = quality factor used in stopping condition for
%                             the optimisation
%       *.power_factor_s    = integral of the squared normalised waveform 
%                             of one sub-pulse (units of seconds)
%       *.constraints       = substructure containing limits for 
%                             optimisation constraints and 'VOPs' used to 
%                             evaluate local SAR, global SAR and average 
%                             power
%
%           OUTPUTS:
%   wopt         = optimised complex weights for PUSH pulse (units of V);
%                  returns complex matrix with dimensions (Nch x Nsp)
%   fval         = cost function value of the optimisation solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) David Leitao & Raphael Tomi-Tricot 2019. King's College London
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(opt_options,'rng_seed') && ~isempty(opt_options.rng_seed)
    rng(opt_options.rng_seed)
end

[Nr, Nch] = size(tx);

%%% Transformation from complex input to real/imaginary input
cp2ri = [eye(Nch), 1i*eye(Nch)];

%%% Transform tx to take real/imaginary input
txri = tx * cp2ri;

%%% Equal spatial weights (if not provided)
if isempty(opt_options.weights)
    W = ones(1, Nr); 
else
    W = opt_options.weights(:).';
end

%%% Intialisation
if isempty(opt_options.w0)
    if opt_options.multistart
        w0 = cell(opt_options.multistart_trials,1);
        for tt=1:opt_options.multistart_trials
            w0{tt} = 2*(rand(Nch*Nsp,1) + 1i*rand(Nch*Nsp,1) - (1+1i)/2) / sqrt(2) * opt_options.base_voltage; 
            w0{tt} = [real(w0{tt}); imag(w0{tt})];
        end
    else
        w0{1} = Sys.txScaleFactor * opt_options.base_voltage;
        w0{1} = repmat([real(w0{1}); imag(w0{1})],[Nsp 1]);
    end
else
    w0 = opt_options.w0;
    if iscell(w0) && numel(w0)>1
        opt_options.multistart        = true;
        opt_options.multistart_trials = numel(w0);
    else
        opt_options.multistart        = false;
    end
end


%% Set optimization parameters

A = []; 
b = [];

cost_normalising = mean(target)^2;
cost =@(x) cost_function(x, target);

costHistory = [];

nonlincon_limits = cat(1,opt_options.constraints.limits.lSAR,...
                         opt_options.constraints.limits.gSAR,...
                         opt_options.constraints.limits.Pmax,...
                         opt_options.constraints.limits.Vmax*ones(Nsp*Nch,1));

options = optimoptions('fmincon',...
                       'Display',opt_options.verbose,...
                       'MaxIter',opt_options.Niter,...
                       'MaxFunEvals',Inf,...
                       'Algorithm','interior-point',... 
                       'SpecifyObjectiveGradient', true,...
                       'SpecifyConstraintGradient',false,...                                
                       'CheckGradients', false,...
                       'InitBarrierParam',1e-1,...
                       'UseParallel', false,...
                       'FiniteDifferenceType','central',...
                       'OutputFcn',@outfun_fmin);


%% Execute fmincon

all_wopt = cell(opt_options.multistart_trials,1); 
all_fval = zeros(opt_options.multistart_trials,1);
for mm=1:opt_options.multistart_trials
    [all_wopt{mm},all_fval(mm)] = fmincon(cost,w0{mm},A,b,[],[],[],[],@(x) nonlincon(x),options);
    costHistory = [];
end
idx_best = find(all_fval==min(all_fval),1,'first');
wopt = cp2ri * reshape(all_wopt{idx_best},[2*Nch Nsp]);
fval = all_fval(idx_best);


%% Cost function and its derivative

    function [f, grad_f] = cost_function(x, target)
        
        x = reshape(x,[2*Nch Nsp]);
        b1 = txri * x;
        b1sq = conj(b1) .* b1;
        b1rms = sqrt(sum(b1sq,2));
        diff_b1rms = b1rms - target;
        SQdiff_b1rms = diff_b1rms.^2;
        f = W * SQdiff_b1rms / sum(W) / cost_normalising;   
                
        if nargout>1
            for pp=1:Nsp
                grad_b1rms(:,(pp-1)*2*Nch+1:pp*2*Nch) = bsxfun(@times,real(b1(:,pp))./b1rms,real(txri)) + bsxfun(@times,imag(b1(:,pp))./b1rms,imag(txri)); % == real( conj(b) .* txri ./ b1rms )
            end
            grad_f = 2 * (W * bsxfun(@times,diff_b1rms,grad_b1rms) / cost_normalising / sum(W)).';    
        end
        
    end


%% Non-linear constraints

	function [c,ceq] = nonlincon(x)    
        
        x = reshape(x,[2*Nch Nsp]);
        b_coeff = transpose(x(1:Nch,:)+1i*x(Nch+1:end,:));
        b_coeff = b_coeff(:);
        
        cval = optimisation_constraints(opt_options.constraints.VOP,b_coeff,Nch,opt_options.power_factor_s);
        
        c = cval - nonlincon_limits;
        ceq = [];
         
    end

%% Own stop condition

    function stop = outfun_fmin(x, optimValues, state) 
        stop = false; 
        switch state
            case 'init'
            case 'iter'    
                
                % get violations of constraints:
                [c,ceq] = nonlincon(x);
                
                violations = length(find(cat(1, [c,ceq]) > 0));
                
                value = optimValues.fval;

                costHistory  = [costHistory,  value];  
                
                if numel(costHistory)>=4
                    s = std(costHistory((end-3):end));
                    q = costHistory(4)* 10^(-opt_options.quality);
                    if s <= q
                        stop = true;
                    end
                end
                
                if violations > 0    % number of constraint violations
                    stop = false;    % don't stop unless all constraints are valid
                end
            otherwise
        end
    end

end