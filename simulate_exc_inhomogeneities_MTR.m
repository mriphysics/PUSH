%%% Simulate MTR with different excitation flip angle and B1rms profiles

clearvars; close all; clc;

% if isempty(gcp('nocreate'))
%     c = parcluster('local');
%     c.NumWorkers = 10;
%     parpool(c, c.NumWorkers);
% end

%% Load transmit maps from optimized solution file

filename = './bin/PUSH_optimisation_2D_slice_12.mat';

aux = load(filename);
tx = aux.tx;
txScaleFactor = aux.Sys.txScaleFactor;
mask = aux.mask;
dims = aux.dims;
Nr = aux.Nr;

%% Define acquisition settings and tissue parameters

dt = 1e-6; % dwell time [s]
TR = 22e-3; % repetition time [s]
alpha_exc = [5 15] * pi/180; % flip angle of the excitation pulse [rad]
tau_exc = 100e-6; % duration of the excitation pulse [s]
tau_sat = 4e-3; % duration of the saturation pulse [s]
Delta_Hz = 2e3; % offset frequency of the saturation pulse [Hz]
pulse_sep = 3e-3; % time between end of saturation pulse and beginning of excitation pulse [s]

target_avgB1rms = 1; % target beta [uT]

Nfa = numel(alpha_exc);

%%% Create nominal RF objects
for aa = 1:Nfa
    b1pulse{aa}.exc = gen_B1_pulse(alpha_exc(aa), tau_exc, dt, 'Rect');
    b1pulse{aa}.sat = gen_B1_pulse(1, tau_sat, dt, 'Gaussian') .* exp(1i*2*pi*Delta_Hz*(0:round(tau_sat/dt)-1)'*dt);
    b1pulse{aa}.sat = b1pulse{aa}.sat ./ sqrt(b1sqrd_integral(b1pulse{aa}.sat,TR,tau_sat,dt)); % normalise saturation pulse to units B1rms
end

%%% Create strucutre with tissue parameters
tissuepars = init_tissue('WM_7T'); 


%% Simulate MTR maps

MTRimg =@(x,x0) reshape(100*(x0 - x)./x0, dims);

MTR_maps = cell(4,Nfa); Msat_imgs = cell(4,Nfa); Mref_imgs = cell(4,Nfa);

% calculate normalised peak B1+ for CP mode (used in excitation pulse)
CPmode_norm_peakB1 = tx*txScaleFactor; CPmode_norm_peakB1 = CPmode_norm_peakB1 ./ mean(abs(CPmode_norm_peakB1(mask(:))));

noise_std = 0;

for aa=1:Nfa
    aux_sat00 = zeros(dims);  aux_ref00 = zeros(dims);
    aux_sat10 = zeros(dims);  aux_ref10 = zeros(dims);
    aux_sat01 = zeros(dims);  aux_ref01 = zeros(dims);  
    aux_sat11 = zeros(dims);  aux_ref11 = zeros(dims);
    for rr=1:Nr
        if ~mask(rr)
            continue;
        end           
        %%% CP flip angle, CP saturation
        aux_b1pulse = b1pulse{aa};
        aux_b1pulse.sat = aux_b1pulse.sat * target_avgB1rms;
        aux_b1pulse.exc = CPmode_norm_peakB1(rr) .* b1pulse{aa}.exc;
        aux_sat00(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars));
        aux_b1pulse.sat = 0;
        aux_ref00(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars));

        %%% CP flip angle, CP saturation
        aux_b1pulse = b1pulse{aa};
        aux_b1pulse.sat = aux_b1pulse.sat * target_avgB1rms;
        aux_b1pulse.exc = CPmode_norm_peakB1(rr) .* b1pulse{aa}.exc;
        aux_sat10(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars, 'FlipAngUni', b1pulse{aa}.exc));
        aux_b1pulse.sat = 0;
        aux_ref10(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars, 'FlipAngUni', b1pulse{aa}.exc));

        %%% uniform flip angle, CP saturation
        aux_b1pulse = b1pulse{aa};
        aux_b1pulse.sat = aux_b1pulse.sat * target_avgB1rms;
        aux_b1pulse.exc = CPmode_norm_peakB1(rr) .* b1pulse{aa}.exc;
        aux_sat10(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars, 'FlipAngUni', b1pulse{aa}.exc));
        aux_b1pulse.sat = 0;
        aux_ref10(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars, 'FlipAngUni', b1pulse{aa}.exc));

        %%% CP flip angle, uniform saturation
        aux_b1pulse = b1pulse{aa};
        aux_b1pulse.sat = aux_b1pulse.sat * target_avgB1rms;
        aux_b1pulse.exc = CPmode_norm_peakB1(rr) .* b1pulse{aa}.exc;
        aux_sat01(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars, 'AvgSatUni', b1pulse{aa}.exc));
        aux_b1pulse.sat = 0;
        aux_ref01(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars, 'AvgSatUni', b1pulse{aa}.exc));

        %%% uniform flip angle, uniform saturation
        aux_b1pulse = b1pulse{aa};
        aux_b1pulse.sat = aux_b1pulse.sat * target_avgB1rms;
        aux_b1pulse.exc = CPmode_norm_peakB1(rr) .* b1pulse{aa}.exc;
        aux_sat11(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars, 'FlipAngUni', b1pulse{aa}.exc, 'AvgSatUni', b1pulse{aa}.exc));
        aux_b1pulse.sat = 0;
        aux_ref11(rr) = abs([1 1i 0 0] * SPGR_BSB(aux_b1pulse, dt, Delta_Hz, TR, pulse_sep, tissuepars, 'FlipAngUni', b1pulse{aa}.exc, 'AvgSatUni', b1pulse{aa}.exc));
    end
    Msat_imgs{1,aa} = aux_sat00; Mref_imgs{1,aa} = aux_ref00;
    Msat_imgs{2,aa} = aux_sat01; Mref_imgs{2,aa} = aux_ref01;
    Msat_imgs{3,aa} = aux_sat10; Mref_imgs{3,aa} = aux_ref10;
    Msat_imgs{4,aa} = aux_sat11; Mref_imgs{4,aa} = aux_ref11;

    MTR_maps{1,aa} = MTRimg(Msat_imgs{1,aa}+randn(dims)*noise_std,Mref_imgs{1,aa}+randn(dims)*noise_std);
    MTR_maps{2,aa} = MTRimg(Msat_imgs{2,aa}+randn(dims)*noise_std,Mref_imgs{2,aa}+randn(dims)*noise_std);
    MTR_maps{3,aa} = MTRimg(Msat_imgs{3,aa}+randn(dims)*noise_std,Mref_imgs{3,aa}+randn(dims)*noise_std);
    MTR_maps{4,aa} = MTRimg(Msat_imgs{4,aa}+randn(dims)*noise_std,Mref_imgs{4,aa}+randn(dims)*noise_std);
end

% save('./bin/MTR_simulation_excitation_inhomogeneity')    

%% Plot results

Nzr1 = 4;
Nzr2 = 10;
Nzr3 = 10;

auxMTR_1 = [zeros(Nzr1,88+Nzr3); rot90(mask.*MTR_maps{1,1},-1),  zeros(50,Nzr3), rot90(mask.*MTR_maps{3,1},-1);  zeros(Nzr2,88+Nzr3); rot90(mask.*MTR_maps{2,1},-1),  zeros(50,Nzr3),  rot90(mask.*MTR_maps{4,1},-1); zeros(Nzr2,88+Nzr3)];
auxMTR_2 = [zeros(Nzr1,88+Nzr3); rot90(mask.*MTR_maps{1,2},-1), zeros(50,Nzr3), rot90(mask.*MTR_maps{3,2},-1); zeros(Nzr2,88+Nzr3); rot90(mask.*MTR_maps{2,2},-1), zeros(50,Nzr3), rot90(mask.*MTR_maps{4,2},-1); zeros(Nzr2,88+Nzr3)];

figure; set(gcf,'color',[1 1 1],'position',[500 400 900 450])
hsp = subplot(1,2,1); 
imagesc(auxMTR_1); colorbar; axis image; colormap('inferno'); hcb = colorbar; xticks([]); yticks([]); 
caxis([0 30]); hcb.Ticks = 0:10:30; hcb.TickLabels = {'0%','10%','20%','30%'}; hcb.FontSize = 14; hcb.Location = 'southoutside';
hsp.Position = [0.1300    0.060    0.372    0.8150]; hcb.Position(2) = 0.06;
text(44+Nzr3/2,-15,'Excitation pulse \alpha','fontsize',14,'fontweight','bold','horizontalalignment','center')
text(22,-5,'Inhomogeneous','fontsize',12,'fontweight','bold','horizontalalignment','center')
text(66+Nzr3,-5,'Homogeneous','fontsize',12,'fontweight','bold','horizontalalignment','center')
text(-15,Nzr1+50+Nzr2/2,'Excitation pulse \beta','fontsize',14.5,'fontweight','bold','horizontalalignment','center','rotation',90)
text(-5,Nzr1+25,'Inhomogeneous','fontsize',12.5,'fontweight','bold','horizontalalignment','center','rotation',90)
text(-5,Nzr1+Nzr2+75,'Homogeneous','fontsize',12.5,'fontweight','bold','horizontalalignment','center','rotation',90)
text(-25,-15,'(A)','fontsize',18,'fontweight','bold','horizontalalignment','center')

hsp = subplot(1,2,2); 
imagesc(auxMTR_2); colorbar; axis image; colormap('inferno'); hcb = colorbar; xticks([]); yticks([]); 
caxis([0 15]); hcb.Ticks = 0:5:15; hcb.TickLabels = {'0%','5%','10%','15%'}; hcb.FontSize = 14; hcb.Location = 'southoutside';
hsp.Position = [0.6300    0.060    0.372    0.8150]; hcb.Position(2) = 0.06;
text(44+Nzr3/2,-15,'Excitation pulse \alpha','fontsize',14,'fontweight','bold','horizontalalignment','center')
text(22,-5,'Inhomogeneous','fontsize',12,'fontweight','bold','horizontalalignment','center')
text(66+Nzr3,-5,'Homogeneous','fontsize',12,'fontweight','bold','horizontalalignment','center')
text(-15,Nzr1+50+Nzr2/2,'Excitation pulse \beta','fontsize',14.5,'fontweight','bold','horizontalalignment','center','rotation',90)
text(-5,Nzr1+25,'Inhomogeneous','fontsize',12.5,'fontweight','bold','horizontalalignment','center','rotation',90)
text(-5,Nzr1+Nzr2+75,'Homogeneous','fontsize',12.5,'fontweight','bold','horizontalalignment','center','rotation',90)
text(-25,-15,'(B)','fontsize',18,'fontweight','bold','horizontalalignment','center')
