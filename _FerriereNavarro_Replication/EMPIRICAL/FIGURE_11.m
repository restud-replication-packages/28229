clear; close all; clc

%% directories

cdir = pwd;
odir = ('output');
fdir = ('figures');

%% load data

cd(odir)
IRF_g_BP_LINEAR  = readtable('LPM_FIGURE_11_lngov.csv','Delimiter',',');
IRFG_BP.m    = table2array(IRF_g_BP_LINEAR(:,1));
IRFG_BP.m_se = table2array(IRF_g_BP_LINEAR(:,2));
IRFG_BP.m_lb = table2array(IRF_g_BP_LINEAR(:,3));
IRFG_BP.m_ub = table2array(IRF_g_BP_LINEAR(:,4));
HORIZ        = table2array(IRF_g_BP_LINEAR(:,6));
H            = max(HORIZ);

MULT_DEF_LINEAR  = readtable('LPM_FIGURE_11_rdef.csv','Delimiter',',');
MULTDF.m     = table2array(MULT_DEF_LINEAR(:,1));
MULTDF.m_se  = table2array(MULT_DEF_LINEAR(:,2));
MULTDF.m_lb  = table2array(MULT_DEF_LINEAR(:,3));
MULTDF.m_ub  = table2array(MULT_DEF_LINEAR(:,4));
MULTDF.HORIZ = table2array(MULT_DEF_LINEAR(:,6));
MULTDF.H     = max(HORIZ);

cd(cdir)

%% Plots

save_graph = 1;

xlb = 0; xub = 15;

figure(13)
subplot(1,2,1)
fill([0:H, H:-1:0],[IRFG_BP.m_lb(1:H+1)', IRFG_BP.m_ub(H+1:-1:1)'],[.4 .6 .6],'edgecolor',[.4 .6 .6])
hold on
yy = plot(0:H,IRFG_BP.m,'k','LineWidth',3.15);
% plot(0:H,gmodel,'color',[0.4 0.4 0.4],'LineWidth',3.15)
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Spending Response','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
% leg = legend(yy,'Linear');
% set(leg,'Interpreter','LaTex','Fontsize',17)
% legend boxoff
xlim([xlb xub])
ylim([-0.1 1.55])
hold off

subplot(1,2,2)
fill([0:H, H:-1:0],[MULTDF.m_lb(1:H+1)', MULTDF.m_ub(H+1:-1:1)'],[.4 .6 .6],'edgecolor',[.4 .6 .6])
hold on
yy = plot(0:H,MULTDF.m(1:H+1),'k','LineWidth',3.15);
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Deficit Multiplier','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
% leg = legend(yy,'Linear');
% set(leg,'Interpreter','LaTex','Fontsize',17)
% legend boxoff
xlim([xlb xub])
ylim([0.35 0.80])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_11','-dpng','-r0')
    cd(cdir)
end


