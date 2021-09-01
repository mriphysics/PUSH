%%% Optimisation of PUSH pulses for 2D or 3D imaging

clearvars; close all; clc; 

% if isempty(gcp('nocreate'))
%     c = parcluster('local');
%     c.NumWorkers = 10;
%     parpool(c, c.NumWorkers);
% end

%% Load TX maps and BET mask; select ROI

load('.\maps\TXmaps.mat')
dims = size(txmaps);
Nch = dims(4); dims = dims(1:3); Nr = prod(dims);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >> select here slices for optimisation 
% numel(slices)==1 -> 2D optimisation
% numel(slices)>1  -> 3D optimisation
slices = 1:dims(3);
dims(3) = numel(slices); Nr = prod(dims);
if dims(3)>1
    is3D = true;
else
    is3D = false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tx = reshape(txmaps(:,:,slices,:),[Nr Nch]);

load('.\maps\BETmask.mat');
mask = mask(:,:,slices);


%% Set optimisation options and define (SAR & hardware) constraints

opt_options.w0                = [];
opt_options.Niter             = 1e3;
opt_options.quality           = 6;
opt_options.verbose           = 'off'; %{'off', 'final', 'iter}
opt_options.multistart        = true;
opt_options.multistart_trials = 1e1;
opt_options.weights           = [];
opt_options.rng_seed          = 'default';

TR  = 22e-3; %[s]
tau = 4e-3;  %[s]
dt  = 10e-6; %[s] dwell time
b1waveform = gen_B1_pulse(1, tau, dt, 'Gaussian'); b1waveform = b1waveform ./ max(abs(b1waveform));
avg2peak = sqrt( TR / (sum(abs(b1waveform).^2) * dt) ); %scaling from average B1rms to peak B1rms


load('./maps/ConstraintData.mat')

opt_options.power_factor_s = sum(abs(b1waveform).^2)*dt; %for SAR calculations
opt_options.constraints = ConstraintData;

CPmode_limits = cat(1,opt_options.constraints.limits.lSAR,...
                      opt_options.constraints.limits.gSAR,...
                      opt_options.constraints.limits.Pmax,...
                      opt_options.constraints.limits.Vmax*ones(Nch,1));
                  
CPmode_w0 = exp(1i*2*pi*(0:Nch-1)/Nch)' ./ sqrt(Nch);

                  
%% Optimisation

if is3D
    filename = './bin/PUSH_optimisation_3D.mat'; 
else
    filename = ['./bin/PUSH_optimisation_2D_slice_',num2str(slices),'.mat'];
end

if ~exist(filename,'file')
    target_avgB1rms = 0.1:0.1:2.0;      Ntg = numel(target_avgB1rms); 
    subpulses = 1:3;                    Nsp = numel(subpulses);

    all_wopt = cell(Nsp,Ntg); 
    all_fval = zeros(opt_options.multistart_trials,Nsp,Ntg); 
    all_runtime = zeros(Nsp,Ntg);

    CPmode_peakB1rms = abs(tx * CPmode_w0); 
    CPmode_wopt = cell(Ntg,1);

    parfor tt=1:Ntg

        % create copy of opt_options to allow multi-threading 
        % /!\ keep opt_options unchanged inside parfor /!\
        cp_opt_options = opt_options;

        % optimal voltage for CP mode (minimises the cost function):
        cp_opt_options.base_voltage = (target_avgB1rms(tt) * avg2peak) * sum(CPmode_peakB1rms(mask)) / sum(CPmode_peakB1rms(mask).^2);
        CPmode_wopt{tt} = cp_opt_options.base_voltage * CPmode_w0;
        % check if optimal CP voltage is within SED limits
        cval = optimisation_constraints(ConstraintData.VOP,CPmode_wopt{tt},Nch,opt_options.power_factor_s);
        relative_c = cval ./ CPmode_limits; 
        % if above SED limits, downweight CP voltage until it's feasible
        while max(relative_c)>1
            CPmode_wopt{tt} = 0.9999 * CPmode_wopt{tt};
            cval = optimisation_constraints(ConstraintData.VOP,CPmode_wopt{tt},Nch,cp_opt_options.power_factor_s);
            relative_c = cval ./ CPmode_limits; 
        end
        
        % optimise PUSH pulses
        target_peakB1rms = target_avgB1rms(tt) * avg2peak * ones(dims);
        for ss=1:Nsp
            tic
            [all_wopt{ss,tt}, all_fval(:,ss,tt)] = optimise_PUSH_B1rms(tx(mask(:),:), target_peakB1rms(mask(:)), subpulses(ss), cp_opt_options);
            all_runtime(ss,tt) = toc;                                           
        end
        fprintf(1,'Finished optimisations using B1rms target = %.2f uT.\n',target_avgB1rms(tt))
    end
    
    save(filename)
    
else
    
    load(filename)
    
end


%% Compare NRMSE and SAR as a function of Nsp and target B1rms

NRMSE =@(x,x0) sqrt( mean((x-x0).^2) ) / mean(x0);

NRMSE_B1rms_PUSH   = zeros(Nsp,Ntg);
NRMSE_B1rms_CPmode = zeros(Ntg,1);
parfor tt=1:Ntg
    for ss=1:Nsp
        aux_B1rms = reshape(sqrt(sum(abs(tx * all_wopt{ss,tt}).^2,2)),dims) / avg2peak;
        NRMSE_B1rms_PUSH(ss,tt) = NRMSE(aux_B1rms(mask(:)), target_avgB1rms(tt));
    end
    
    aux_B1rms = reshape(sqrt(sum(abs(tx * CPmode_wopt{tt}).^2,2)),dims) / avg2peak;
    NRMSE_B1rms_CPmode(tt) = NRMSE(aux_B1rms(mask(:)), target_avgB1rms(tt));
end

if is3D
    filename2 = './bin/NRMSE_PUSH_optimisation_3D.mat'; 
else
    filename2 = ['./bin/NRMSE_PUSH_optimisation_2D_slice_',num2str(slices),'.mat'];
end

if ~exist(filename2,'file')
    save(filename2,'NRMSE_B1rms_PUSH','NRMSE_B1rms_CPmode','target_avgB1rms','subpulses')
end


%% Plot B1rms maps for some target beta

idx_target_ROI = 4:3:13; 

avgB1rms = zeros(numel(idx_target_ROI), Nsp+1);
stdB1rms = zeros(numel(idx_target_ROI), Nsp+1);

figure; set(gcf,'Units','normalized','Color','w');
if is3D; set(gcf,'Outerposition',[0.05 0.15 0.8 0.6]); else; set(gcf,'Outerposition',[0.3 0.15 0.375 0.8]); end
for tt=1:numel(idx_target_ROI)
    all_b1rms = [];
    aux = mask .* reshape(abs(tx * CPmode_wopt{idx_target_ROI(tt)}),dims) / avg2peak;
    avgB1rms(tt,1) = mean(aux(mask(:))); stdB1rms(tt,1) = std(aux(mask(:)));
    if is3D
        tra = flipud(rot90(aux(:,:,dims(3)/2)));
        cor = [zeros((dims(2)-dims(3))/2, dims(1)); rot90(squeeze(aux(:,dims(2)/2,:))); zeros((dims(2)-dims(3))/2, dims(1))];
        sag = [zeros((dims(2)-dims(3))/2, dims(2)); rot90(squeeze(aux(dims(1)/2,:,:))); zeros((dims(2)-dims(3))/2, dims(2))]; 
        aux = cat(2, tra, cor, sag);
        all_b1rms = cat(1, all_b1rms, aux, zeros(size(aux,1)/10, size(aux,2)));
    else
        all_b1rms = cat(1, all_b1rms, flipud(rot90(aux)), zeros(round(dims(1)/5), dims(1)));
    end
    for ss=1:Nsp
        aux = mask .* reshape(sqrt(sum(abs(tx * all_wopt{ss,idx_target_ROI(tt)}).^2,2)),dims) / avg2peak;
        avgB1rms(tt,ss+1) = mean(aux(mask(:))); stdB1rms(tt,ss+1) = std(aux(mask(:)));
        if is3D
            tra = flipud(rot90(aux(:,:,dims(3)/2)));
            cor = [zeros((dims(2)-dims(3))/2, dims(1)); rot90(squeeze(aux(:,dims(2)/2,:))); zeros((dims(2)-dims(3))/2, dims(1))];
            sag = [zeros((dims(2)-dims(3))/2, dims(2)); rot90(squeeze(aux(dims(1)/2,:,:))); zeros((dims(2)-dims(3))/2, dims(2))]; 
            aux = cat(2, tra, cor, sag);
            all_b1rms = cat(1, all_b1rms, aux, zeros(size(aux,1)/10, size(aux,2)));
        else
            all_b1rms = cat(1, all_b1rms, flipud(rot90(aux)), zeros(round(dims(1)/5), dims(1)));
        end
    end
    
    hsp(tt) = subplot(1, numel(idx_target_ROI),tt);
    imagesc(all_b1rms, round([100/3 400/3]*target_avgB1rms(idx_target_ROI(tt)))/100)
    axis image; colorcet('L3'); set(gca,'Visible','off'); hcb(tt) = colorbar;
    hcb(tt).Location = 'southoutside'; 
    
    if is3D
        hcb(tt).FontSize = 12; 
        text(0.5,1.05,['\beta = ',num2str(target_avgB1rms(idx_target_ROI(tt))),'\mu{}T'],'Units','normalized','Fontsize',16,'Fontweight','bold','HorizontalAlignment','center')
        hcb(tt).YTick = round([115/3 500/6 385/3]*target_avgB1rms(idx_target_ROI(tt)))/100;
        hcb(tt).YTickLabel = [num2str([round([100/3 500/6 400/3]*target_avgB1rms(idx_target_ROI(tt)))/100]'),repmat('\muT',[3 1])];
        hcb(tt).YLim = round([100/3 400/3]*target_avgB1rms(idx_target_ROI(tt)))/100;
    else
        hcb(tt).FontSize = 11;
        text(0.5,1.04,['\beta = ',num2str(target_avgB1rms(idx_target_ROI(tt))),'\mu{}T'],'Units','normalized','Fontsize',16,'Fontweight','bold','HorizontalAlignment','center')
        hcb(tt).YTick = round([130/3 500/6 370/3]*target_avgB1rms(idx_target_ROI(tt)))/100;
        hcb(tt).YTickLabel = round([100/3 500/6 400/3]*target_avgB1rms(idx_target_ROI(tt)))/100;
        hcb(tt).YLim = round([100/3 400/3]*target_avgB1rms(idx_target_ROI(tt)))/100;
    end
    
    if is3D
        hsp(tt).Position(2) = 0.1; hsp(tt).Position(4) = 0.84;
    else
        hsp(tt).Position(2:4) = [0.085, 0.18, 0.87];
    end
    
    if tt==1
        if is3D; xshift = -0.22; else; xshift = -0.5; end
        text(xshift,21.5/24,'CP mode',   'Units','normalized','Fontsize',16,'Fontweight','bold','HorizontalAlignment','center')
        text(xshift,15.5/24,{'PUSH-1';'(static shim)'},'Units','normalized','Fontsize',16,'Fontweight','bold','HorizontalAlignment','center')
        text(xshift,9.5/24, 'PUSH-2','Units','normalized','Fontsize',16,'Fontweight','bold','HorizontalAlignment','center')
        text(xshift,3.5/24, 'PUSH-3','Units','normalized','Fontsize',16,'Fontweight','bold','HorizontalAlignment','center')
    end
    
    for ss=1:Nsp+1
        text(0.5,(24.5-ss*6)/24,[num2str(avgB1rms(tt,ss),'%.2f'),char(177),num2str(stdB1rms(tt,ss),'%.2f'),'\mu{}T'],'Units','normalized','Fontsize',13,'HorizontalAlignment','center','Color','w')
    end
    
    if ~is3D
        text(0.5,-0.082,'\mu{}T','Units','normalized','Fontsize',11,'HorizontalAlignment','center')
    end
end

for tt=numel(idx_target_ROI):-1:1
    if is3D
        hsp(tt).Position(1) = 0.31 + (tt-2)*0.23;       
        hsp(tt).Position(3) = 0.22;
        hcb(tt).Position = [hsp(tt).Position(1)+0.0135 0.055 0.193 0.035];
    else
        hsp(tt).Position(1) = 0.39 + (tt-2)*0.21;    
        hcb(tt).Position = [hsp(tt).Position(1) 0.055 0.18 0.0259];
    end
end
