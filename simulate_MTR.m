%%% Simulate MTR maps from optimised solutions

clearvars; close all; clc;

% if isempty(gcp('nocreate'))
%     c = parcluster('local');
%     c.NumWorkers = 10;
%     parpool(c, c.NumWorkers);
% end

%% Load optimised solutions to use in MTR simulation

% filename = './bin/PUSH_optimisation_3D.mat'; 
filename = './bin/PUSH_optimisation_2D_slice_12.mat';

load(filename)

%% Define acquisition settings and tissue parameters

dt = 10e-6;
alpha_exc = 5 * pi/180; % flip angle of the excitation pulse
tau_exc = 2e-3; % duration of the excitation pulse 
tau_sat = tau;  % duration of the saturation pulse
Delta_Hz = 2e3; % offset frequency of the saturation pulse
pulse_sep = 3e-3; % time between end of saturation pulse and beginning of excitation pulse

%%% Create nominal RF objects
b1pulse.exc = gen_B1_pulse(alpha_exc, tau_exc, dt, 'Gaussian');
b1pulse.sat = gen_B1_pulse(1, tau_sat, dt, 'Gaussian') .* exp(1i*2*pi*Delta_Hz*(0:round(tau_sat/dt)-1)'*dt);
b1pulse.sat = b1pulse.sat ./ max(abs(b1pulse.sat)); % normalise saturation pulse

%%% Create strucutre with tissue parameters
tissuepars = init_tissue('WM_7T'); 


%% Simulate MTR maps

MTRimg =@(x,x0) reshape(100*(x0 - x)./x0, dims);

MTR_PUSH   = cell(Nsp,Ntg); Msat_PUSH   = cell(Nsp,Ntg); Mref_PUSH = cell(Nsp,Ntg);
MTR_CPmode = cell(Ntg,1);   Msat_CPmode = cell(Ntg,1);   Mref_CPmode = cell(Ntg,1);

% calculate normalised peak B1+ for CP mode (used in excitation pulse)
CPmode_norm_peakB1 = tx*Sys.txScaleFactor; CPmode_norm_peakB1 = CPmode_norm_peakB1 ./ mean(abs(CPmode_norm_peakB1(mask(:))));

parfor tt=1:Ntg
    for ss=1:Nsp
        B1peak = tx*all_wopt{ss,tt};

        Msat_PUSH{ss,tt} = zeros(dims); Mref_PUSH{ss,tt} = zeros(dims);
        for rr=1:Nr
            if ~mask(rr)
                continue;
            end
            
            % create RF pulses waveforms from the nominal waveforms
            aux_b1pulse = b1pulse;
            % scale excitation pulse waveform
            aux_b1pulse.exc = CPmode_norm_peakB1(rr) .* aux_b1pulse.exc;
            % create saturation pulse waveform with all sub-pulses scaled
            aux_b1pulse.sat = []; 
            for pp=1:subpulses(ss)
                aux_b1pulse.sat = cat(1, aux_b1pulse.sat, B1peak(rr,pp).*b1pulse.sat);
            end
            
            % calculate SPGR steady-state signal with saturation pulse
            Msat_PUSH{ss,tt}(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars));
            
            % calculate SPGR steady-state signal without saturation pulse
            aux_b1pulse.sat = 0;
            Mref_PUSH{ss,tt}(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars));
        end
        
        MTR_PUSH{ss,tt} = MTRimg(Msat_PUSH{ss,tt},Mref_PUSH{ss,tt});
        
    end
    
    % Repeat calculation for CP mode:
    B1peak = tx*CPmode_wopt{tt};

    Msat_CPmode{tt} = zeros(dims); Mref_CPmode{tt} = zeros(dims);
    for rr=1:Nr
        if ~mask(rr)
            continue;
        end
            
        % create RF pulses waveforms from the nominal waveforms
        aux_b1pulse = b1pulse;
        % scale excitation pulse waveform
        aux_b1pulse.exc = CPmode_norm_peakB1(rr) .* aux_b1pulse.exc;
        % scale saturation pulse waveform
        aux_b1pulse.sat = B1peak(rr) .* aux_b1pulse.sat;
            
        % calculate SPGR steady-state signal with saturation pulse
        Msat_CPmode{tt}(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars));

        % calculate SPGR steady-state signal without saturation pulse
        aux_b1pulse.sat = 0;
        Mref_CPmode{tt}(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars));
    end
    
    MTR_CPmode{tt} = MTRimg(Msat_CPmode{tt},Mref_CPmode{tt});
        
end

if is3D
    save('./bin/MTR_simulation_3D','MTR_PUSH','MTR_CPmode')    
else
    save(['./bin/MTR_simulation_2D_slice_',num2str(slices)],'MTR_PUSH','MTR_CPmode')    
end

