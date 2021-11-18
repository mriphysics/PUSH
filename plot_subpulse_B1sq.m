%%% Plot mean squared B1+ for each sub-pulse in PUSH solutions

clearvars; close all; clc;

%% Load 2D optimization solutions for one slice (e.g. slice 12)

load('PUSH_optimisation_2D_slice_12.mat');

% select one target (index 10 corresponds to beta = 1uT)
idx_target_ROI = 10;
tt = 1;

% concatenate mean squared B1+ maps in big matrix to plot
Nzrs = 10;
Nspc = 10;
all_b1sq = cell(3,4);

all_b1sq{1,1} =  [[zeros(Nzrs/2,dims(1)); flipud(rot90(mask .* reshape(abs(tx * all_wopt{1,idx_target_ROI(tt)}(:,1)).^2,dims) / avg2peak^2)); zeros(Nzrs/2,dims(1))], zeros(dims(2)+Nzrs,Nspc)];
all_b1sq{1,2} =  zeros(dims(2)+Nzrs,dims(1)+Nspc);
all_b1sq{1,3} =  zeros(dims(2)+Nzrs,dims(1)+Nspc);
all_b1sq{1,4} =  [zeros(Nzrs/2,dims(1)); flipud(rot90(mask .* reshape(abs(tx * all_wopt{1,idx_target_ROI(tt)}(:,1)).^2,dims) / avg2peak^2)); zeros(Nzrs/2,dims(1))];

all_b1sq{2,1} =  [[zeros(Nzrs/2,dims(1)); flipud(rot90(mask .* reshape(abs(tx * all_wopt{2,idx_target_ROI(tt)}(:,1)).^2,dims) / avg2peak^2)); zeros(Nzrs/2,dims(1))], zeros(dims(2)+Nzrs,Nspc)];
all_b1sq{2,2} =  [[zeros(Nzrs/2,dims(1)); flipud(rot90(mask .* reshape(abs(tx * all_wopt{2,idx_target_ROI(tt)}(:,2)).^2,dims) / avg2peak^2)); zeros(Nzrs/2,dims(1))], zeros(dims(2)+Nzrs,Nspc)];
all_b1sq{2,3} =  zeros(dims(2)+Nzrs,dims(1)+Nspc);
all_b1sq{2,4} =  [zeros(Nzrs/2,dims(1)); flipud(rot90(mask .* reshape(sum(abs(tx * all_wopt{2,idx_target_ROI(tt)}).^2,2),dims) / avg2peak^2)); zeros(Nzrs/2,dims(1))];

all_b1sq{3,1} =  [[zeros(Nzrs/2,dims(1)); flipud(rot90(mask .* reshape(abs(tx * all_wopt{3,idx_target_ROI(tt)}(:,1)).^2,dims) / avg2peak^2)); zeros(Nzrs/2,dims(1))], zeros(dims(2)+Nzrs,Nspc)];
all_b1sq{3,2} =  [[zeros(Nzrs/2,dims(1)); flipud(rot90(mask .* reshape(abs(tx * all_wopt{3,idx_target_ROI(tt)}(:,2)).^2,dims) / avg2peak^2)); zeros(Nzrs/2,dims(1))], zeros(dims(2)+Nzrs,Nspc)];
all_b1sq{3,3} =  [[zeros(Nzrs/2,dims(1)); flipud(rot90(mask .* reshape(abs(tx * all_wopt{3,idx_target_ROI(tt)}(:,3)).^2,dims) / avg2peak^2)); zeros(Nzrs/2,dims(1))], zeros(dims(2)+Nzrs,Nspc)];
all_b1sq{3,4} =  [zeros(Nzrs/2,dims(1)); flipud(rot90(mask .* reshape(sum(abs(tx * all_wopt{3,idx_target_ROI(tt)}).^2,2),dims) / avg2peak^2)); zeros(Nzrs/2,dims(1))];

%% Plot figure

figure; set(gcf,'color',[1 1 1],'position',[400 100 850 650])
hsp = subplot(1,1,1); imagesc(cell2mat(all_b1sq),[0 1.5]); axis image; colormap('hot'); 
xticks([]); yticks([]); hcb = colorbar; title(hcb,'\mu{}T^2'); hcb.Ticks = 0:0.5:1.5; hcb.FontSize = 16;

text(1.5*dims(1)+Nspc,-16,'Sub-pulse \langle{}|B_1^+|^2\rangle{}','HorizontalAlignment','center','fontsize',16,'fontweight','bold')
text(0.5*dims(1),-5,'#1','HorizontalAlignment','center','fontsize',16,'fontweight','bold')
text(1.5*dims(1)+Nspc,-5,'#2','HorizontalAlignment','center','fontsize',16,'fontweight','bold')
text(2.5*dims(1)+2*Nspc,-5,'#3','HorizontalAlignment','center','fontsize',16,'fontweight','bold')
text(3.5*dims(1)+3*Nspc,-16,'Total \langle{}|B_1^+|^2\rangle{}','HorizontalAlignment','center','fontsize',16,'fontweight','bold')

text(-20,0.5*(dims(2)+Nzrs),{['PUSH-1'];['(static shim)']},'HorizontalAlignment','center','fontsize',16,'fontweight','bold')
text(3*dims(1)+2.5*Nspc,0.5*(dims(2)+Nzrs),'=','HorizontalAlignment','center','fontsize',24,'fontweight','bold','color',[1 1 1])
text(-20,1.5*(dims(2)+Nzrs),'PUSH-2','HorizontalAlignment','center','fontsize',16,'fontweight','bold')
text(1*dims(1)+0.5*Nspc,1.5*(dims(2)+Nzrs),'+','HorizontalAlignment','center','fontsize',24,'fontweight','bold','color',[1 1 1])
text(3*dims(1)+2.5*Nspc,1.5*(dims(2)+Nzrs),'=','HorizontalAlignment','center','fontsize',24,'fontweight','bold','color',[1 1 1])
text(-20,2.5*(dims(2)+Nzrs),'PUSH-3','HorizontalAlignment','center','fontsize',16,'fontweight','bold')
text(1*dims(1)+0.5*Nspc,2.5*(dims(2)+Nzrs),'+','HorizontalAlignment','center','fontsize',24,'fontweight','bold','color',[1 1 1])
text(2*dims(1)+1.5*Nspc,2.5*(dims(2)+Nzrs),'+','HorizontalAlignment','center','fontsize',24,'fontweight','bold','color',[1 1 1])
text(3*dims(1)+2.5*Nspc,2.5*(dims(2)+Nzrs),'=','HorizontalAlignment','center','fontsize',24,'fontweight','bold','color',[1 1 1])

hsp.Position = [0.150    0.02    0.75    0.88];

% print img/subpule_B1sq.png -dpng -r300
