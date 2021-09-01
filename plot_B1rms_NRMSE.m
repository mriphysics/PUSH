%%% Plot NRMSE of 2D and 3D optimisations' of solutions

clearvars; close all; clc; 

load('./bin/NRMSE_PUSH_optimisation_2D_slice_12.mat')
NRMSE_B1rms_CPmode_2D = NRMSE_B1rms_CPmode; 
NRMSE_B1rms_PUSH_2D   = NRMSE_B1rms_PUSH;
load('./bin/NRMSE_PUSH_optimisation_3D.mat')
NRMSE_B1rms_CPmode_3D = NRMSE_B1rms_CPmode; 
NRMSE_B1rms_PUSH_3D   = NRMSE_B1rms_PUSH;

% find last target B1rms where CP mode voltage was not capped
idx_SEDlimits_2D = find(diff(round(1e4*NRMSE_B1rms_CPmode_2D))>0,1,'first');
idx_SEDlimits_3D = find(diff(round(1e4*NRMSE_B1rms_CPmode_3D))>0,1,'first');

figure; set(gcf,'Color','w','Units','normalized','Outerposition',[0.325 0.05 0.325 0.9]);
c = lines(4);

subplot(2,1,1)
hold on;
plot(target_avgB1rms(1:idx_SEDlimits_2D),NRMSE_B1rms_CPmode_2D(1:idx_SEDlimits_2D),'-','Linewidth',2,'MarkerFaceColor','k','Color','k')
area([target_avgB1rms(idx_SEDlimits_2D), target_avgB1rms(end)],[0.6 0.6],'FaceColor','k','FaceAlpha',0.2,'EdgeColor','k','EdgeAlpha',0.2)
plot(target_avgB1rms(idx_SEDlimits_2D:end),NRMSE_B1rms_CPmode_2D(idx_SEDlimits_2D:end),'--','Linewidth',2,'MarkerFaceColor','k','Color','k')
plot(target_avgB1rms,NRMSE_B1rms_PUSH_2D(1,:),'-o','Linewidth',2,'MarkerFaceColor',c(1,:),'Color',c(1,:),'MarkerSize',6)
plot(target_avgB1rms,NRMSE_B1rms_PUSH_2D(2,:),'-s','Linewidth',2,'MarkerFaceColor',c(2,:),'Color',c(2,:),'MarkerSize',6)
plot(target_avgB1rms,NRMSE_B1rms_PUSH_2D(3,:),'-^','Linewidth',2,'MarkerFaceColor',c(3,:),'Color',c(3,:),'MarkerSize',5)
ylabel('NRMSE of B_1^{rms}')
xlabel('\beta (\mu{}T)')
set(gca,'Fontsize',14)
box on; grid on; ylim([0 0.6])
ax_obj = findobj(gca,'Type','Line');
legend([ax_obj(5),ax_obj(3),ax_obj(2),ax_obj(1)],'CP mode','PUSH-1','PUSH-2','PUSH-3','Location','northwest') %
title('2D imaging')
text(-0.17,1.075,'(A)','Units','normalized','Fontsize',20,'Fontweight','bold')

subplot(2,1,2)
hold on;
plot(target_avgB1rms(1:idx_SEDlimits_3D),NRMSE_B1rms_CPmode_3D(1:idx_SEDlimits_3D),'-','Linewidth',2,'MarkerFaceColor','k','Color','k')
area([target_avgB1rms(idx_SEDlimits_3D), target_avgB1rms(end)],[0.6 0.6],'FaceColor','k','FaceAlpha',0.2,'EdgeColor','k','EdgeAlpha',0.2)
plot(target_avgB1rms(idx_SEDlimits_3D:end),NRMSE_B1rms_CPmode_3D(idx_SEDlimits_3D:end),'--','Linewidth',2,'MarkerFaceColor','k','Color','k')
plot(target_avgB1rms,NRMSE_B1rms_PUSH_3D(1,:),'-o','Linewidth',2,'MarkerFaceColor',c(1,:),'Color',c(1,:),'MarkerSize',6)
plot(target_avgB1rms,NRMSE_B1rms_PUSH_3D(2,:),'-s','Linewidth',2,'MarkerFaceColor',c(2,:),'Color',c(2,:),'MarkerSize',6)
plot(target_avgB1rms,NRMSE_B1rms_PUSH_3D(3,:),'-^','Linewidth',2,'MarkerFaceColor',c(3,:),'Color',c(3,:),'MarkerSize',5)
ylabel('NRMSE of B_1^{rms}')
xlabel('\beta (\mu{}T)')
set(gca,'Fontsize',14)
box on; grid on; ylim([0 0.6])
ax_obj = findobj(gca,'Type','Line');
legend([ax_obj(5),ax_obj(3),ax_obj(2),ax_obj(1)],'CP mode','PUSH-1','PUSH-2','PUSH-3','Location','northwest')
title('3D imaging')
text(-0.17,1.075,'(B)','Units','normalized','Fontsize',20,'Fontweight','bold')

