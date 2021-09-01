%%% Plot MTR of solutions for 2D and 3D optimisations in simulation

clearvars; close all; clc; 

%% Load simulated MTR maps

slices = [6 10 12 14 18];

MTRmaps_2D = cell(numel(slices),1);
MTRmaps_3D = load('./bin/MTR_simulation_3D.mat');

for ii=1:numel(slices)
    filename = ['./bin/MTR_simulation_2D_slice_',num2str(slices(ii)),'.mat'];
    MTRmaps_2D{ii} = load(filename);
end

load('./maps/BETmask.mat');

dims = size(mask);

Nsp = size(MTRmaps_2D{1}.MTR_PUSH,1);

%% Plot MTR maps

tt = 10; % index of B1rms target to plot

avgMTR_2D = zeros(numel(slices),Nsp+1);
stdMTR_2D = zeros(numel(slices),Nsp+1);
avgMTR_3D = zeros(Nsp+1,1);
stdMTR_3D = zeros(Nsp+1,1);

N = numel(slices)*dims(1) + round(dims(2)/2) + 2*dims(1)+dims(2);

%%% concatenate all MTR maps in a single matrix; first CP mode 
aux_fill = [zeros(round(dims(2)/8),numel(slices)*dims(1)), 100*ones(round(dims(2)/8),round(dims(2)/2)), zeros(round(dims(2)/8),2*dims(1)+dims(2))];
tra2D_MTR = cell(Nsp+1,1);
for ii=1:numel(slices)
    tra2D_MTR{ii} = flipud(rot90(MTRmaps_2D{ii}.MTR_CPmode{tt}));
    avgMTR_2D(ii,1) = mean(tra2D_MTR{ii}(flipud(rot90(mask(:,:,slices(ii))))));
    stdMTR_2D(ii,1) = std(tra2D_MTR{ii}(flipud(rot90(mask(:,:,slices(ii))))));
end
tra3D_MTR = flipud(rot90(MTRmaps_3D.MTR_CPmode{tt}(:,:,dims(3)/2)));
cor3D_MTR = [zeros((dims(2)-dims(3))/2, dims(1)); rot90(squeeze(MTRmaps_3D.MTR_CPmode{tt}(:,dims(2)/2,:))); zeros((dims(2)-dims(3))/2, dims(1))];
sag3D_MTR = [zeros((dims(2)-dims(3))/2, dims(2)); rot90(squeeze(MTRmaps_3D.MTR_CPmode{tt}(dims(1)/2,:,:))); zeros((dims(2)-dims(3))/2, dims(2))]; 
avgMTR_3D(tt,1) = mean(MTRmaps_3D.MTR_CPmode{tt}(mask)); 
stdMTR_3D(tt,1) = std(MTRmaps_3D.MTR_CPmode{tt}(mask));
aux_MTR = [];
for ii=1:numel(slices)
    aux_MTR = cat(2, aux_MTR, tra2D_MTR{ii});
end
aux_MTR = cat(2, aux_MTR, 100*ones(dims(2), round(dims(2)/2)), tra3D_MTR, cor3D_MTR, sag3D_MTR);
all_MTR = cat(1, aux_fill, aux_MTR);

% now concatenate PUSH
for ss=1:Nsp
    for ii=1:numel(slices)
        tra2D_MTR{ii} = flipud(rot90(MTRmaps_2D{ii}.MTR_PUSH{ss,tt}));
        avgMTR_2D(ii,ss+1) = mean(tra2D_MTR{ii}(flipud(rot90(mask(:,:,slices(ii))))));
        stdMTR_2D(ii,ss+1) = std(tra2D_MTR{ii}(flipud(rot90(mask(:,:,slices(ii))))));
    end
    tra3D_MTR = flipud(rot90(MTRmaps_3D.MTR_PUSH{ss,tt}(:,:,dims(3)/2)));
    cor3D_MTR = [zeros((dims(2)-dims(3))/2, dims(1)); rot90(squeeze(MTRmaps_3D.MTR_PUSH{ss,tt}(:,dims(2)/2,:))); zeros((dims(2)-dims(3))/2, dims(1))];
    sag3D_MTR = [zeros((dims(2)-dims(3))/2, dims(2)); rot90(squeeze(MTRmaps_3D.MTR_PUSH{ss,tt}(dims(1)/2,:,:))); zeros((dims(2)-dims(3))/2, dims(2))]; 
    avgMTR_3D(tt,ss+1) = mean(MTRmaps_3D.MTR_PUSH{ss,tt}(mask)); 
    stdMTR_3D(tt,ss+1) = std(MTRmaps_3D.MTR_PUSH{ss,tt}(mask));
    aux_MTR = [];
    for ii=1:numel(slices)
        aux_MTR = cat(2, aux_MTR, tra2D_MTR{ii});
    end
    aux_MTR = cat(2, aux_MTR, 100*ones(dims(2), round(dims(2)/2)), tra3D_MTR, cor3D_MTR, sag3D_MTR);        
    all_MTR = cat(1, all_MTR, aux_fill, aux_MTR);
end 
all_MTR = cat(1, all_MTR, aux_fill);

%Note that coronal and sagittal views here appeared squeezed in HF
%direction because acquired voxel size is 4x4x6mm and it is not 
%interpolated here to isotropic resolution for display
figure;
set(gcf,'Units','normalized','Color','w','Outerposition',[0.2 0.2 0.575 0.675])    
hnd1 = subplot(1,1,1);
imagesc(all_MTR,[0 30.01])
axis image; colormap([inferno(1000);[1 1 1]]); set(gca,'Visible','off'); hcb = colorbar;
hcb.Location = 'southoutside'; hcb.FontSize = 13; hcb.Ticks = 0:10:40; hcb.TickLabels = {'0','10','20','30'};
hnd1.Position = [0.125 0.1 0.85 0.86]; hcb.Position(1:3) = [0.1305 0.05 0.84];

text(0.02,0.985,'(A)','Units','normalized','Fontsize',16,'FontWeight','bold','HorizontalAlignment','center','Color','w')
text(0.657,0.985,'(B)','Units','normalized','Fontsize',16,'FontWeight','bold','HorizontalAlignment','center','FontWeight','bold','Color','w')

text(-0.08,28/33,'CP mode','Units','normalized','Fontsize',14,'FontWeight','bold','HorizontalAlignment','center')
text(-0.08,20/33, {'PUSH-1';'(static shim)'},'Units','normalized','Fontsize',14,'FontWeight','bold','HorizontalAlignment','center')
text(-0.08,12/33, 'PUSH-2','Units','normalized','Fontsize',14,'FontWeight','bold','HorizontalAlignment','center')
text(-0.08,4/33,  'PUSH-3','Units','normalized','Fontsize',14,'FontWeight','bold','HorizontalAlignment','center')

text(110/N,1.03,'2D imaging','Units','normalized','Fontsize',14,'FontWeight','bold','HorizontalAlignment','center')
text(314/N,1.03,'3D imaging','Units','normalized','Fontsize',14,'FontWeight','bold','HorizontalAlignment','center')
