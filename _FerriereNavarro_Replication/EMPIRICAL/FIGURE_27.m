clear; close all; clc

%% directories

cdir = pwd;
odir = ('output');
fdir = ('figures');


%% Load data

save_graph = 1;

cd(odir)
IRF_g_BP_LINEAR  = readtable('LPM_FIGURE_27_gov_bp.csv','Delimiter',',');
IRFG_BP.m    = table2array(IRF_g_BP_LINEAR(:,1));
IRFG_BP.m_se = table2array(IRF_g_BP_LINEAR(:,2));
IRFG_BP.m_lb = table2array(IRF_g_BP_LINEAR(:,3));
IRFG_BP.m_ub = table2array(IRF_g_BP_LINEAR(:,4));
HORIZ        = table2array(IRF_g_BP_LINEAR(:,6));
H            = max(HORIZ);


IRF_y_BP_LINEAR  = readtable('LPM_FIGURE_27_gdp_bp.csv','Delimiter',',');
IRFY_BP.m    = table2array(IRF_y_BP_LINEAR(:,1));
IRFY_BP.m_se = table2array(IRF_y_BP_LINEAR(:,2));
IRFY_BP.m_lb = table2array(IRF_y_BP_LINEAR(:,3));
IRFY_BP.m_ub = table2array(IRF_y_BP_LINEAR(:,4));
% HORIZ        = table2array(IRF_g_BP_LINEAR(:,10));
% H            = max(HORIZ);

IRF_g_RZ_LINEAR  = readtable('LPM_FIGURE_27_gov_news.csv','Delimiter',',');
IRFG_RZ.m    = table2array(IRF_g_RZ_LINEAR(:,1));
IRFG_RZ.m_se = table2array(IRF_g_RZ_LINEAR(:,2));
IRFG_RZ.m_lb = table2array(IRF_g_RZ_LINEAR(:,3));
IRFG_RZ.m_ub = table2array(IRF_g_RZ_LINEAR(:,4));

IRF_y_RZ_LINEAR  = readtable('LPM_FIGURE_27_gdp_news.csv','Delimiter',',');
IRFY_RZ.m    = table2array(IRF_y_RZ_LINEAR(:,1));
IRFY_RZ.m_se = table2array(IRF_y_RZ_LINEAR(:,2));
IRFY_RZ.m_lb = table2array(IRF_y_RZ_LINEAR(:,3));
IRFY_RZ.m_ub = table2array(IRF_y_RZ_LINEAR(:,4));

cd(cdir)

zz = [IRFG_BP.m, IRFY_BP.m, IRFG_RZ.m, IRFY_RZ.m];

%% normalizing constant

%---normalize bp
bmax = max(IRFG_BP.m);
IRFG_BP.m = IRFG_BP.m/bmax;  IRFG_BP.m_lb = IRFG_BP.m_lb/bmax;  IRFG_BP.m_ub = IRFG_BP.m_ub/bmax;
IRFY_BP.m = IRFY_BP.m/bmax;  IRFY_BP.m_lb = IRFY_BP.m_lb/bmax;  IRFY_BP.m_ub = IRFY_BP.m_ub/bmax;

%---normalize rz
bmax = max(IRFG_RZ.m);
IRFG_RZ.m = IRFG_RZ.m/bmax;  IRFG_RZ.m_lb = IRFG_RZ.m_lb/bmax;  IRFG_RZ.m_ub = IRFG_RZ.m_ub/bmax;
IRFY_RZ.m = IRFY_RZ.m/bmax;  IRFY_RZ.m_lb = IRFY_RZ.m_lb/bmax;  IRFY_RZ.m_ub = IRFY_RZ.m_ub/bmax;

%% plot spending comparisons 

cBP = [0.75 0.75 0.75];  cBP_l = [0.7 0.5 0.5];
cRZ = [0.75 0.75 0.90];  cRZ_l = [0.5 0.5 0.7];

xlb = 0; xub = 15;

fig=figure(1301); clf;
subplot(1,2,1)
fill([0:H, H:-1:0],[IRFG_BP.m_lb(1:H+1)', IRFG_BP.m_ub(H+1:-1:1)'],cBP,'edgecolor',cBP,'FaceAlpha',0.5); hold on
fill([0:H, H:-1:0],[IRFG_RZ.m_lb(1:H+1)', IRFG_RZ.m_ub(H+1:-1:1)'],cRZ,'edgecolor',cRZ,'FaceAlpha',0.5); 
% plot(0:H,gmodel,'color',[0.10 0.10 0.10],'LineWidth',3.25)
yy = plot(0:H,IRFG_BP.m,'color',cBP_l,'LineWidth',3.15);
yy = plot(0:H,IRFG_RZ.m,'color',cRZ_l,'LineWidth',3.15);

set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Spending Response','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
% leg = legend('BP shock','RZ shock'); % leg = legend('BP shock','RZ shock','Model');
% set(leg,'Interpreter','LaTex','Fontsize',27)
% legend boxoff
xlim([xlb xub])
ylim([-0.1 1.55])
hold off

subplot(1,2,2)
fill([0:H, H:-1:0],[IRFY_BP.m_lb(1:H+1)', IRFY_BP.m_ub(H+1:-1:1)'],cBP,'edgecolor',cBP,'FaceAlpha',0.5); hold on
fill([0:H, H:-1:0],[IRFY_RZ.m_lb(1:H+1)', IRFY_RZ.m_ub(H+1:-1:1)'],cRZ,'edgecolor',cRZ,'FaceAlpha',0.5); 
yy = plot(0:H,IRFY_BP.m,'color',cBP_l,'LineWidth',3.15);
yy = plot(0:H,IRFY_RZ.m,'color',cRZ_l,'LineWidth',3.15);

set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Output Response','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
leg = legend('BP shock','RZ shock');
set(leg,'Interpreter','LaTex','Fontsize',27,'Location','NorthWest')
legend boxoff
xlim([xlb xub])
% ylim([-0.05 0.55])
ylim([-0.10 1.55])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_27','-dpng','-r0')
    cd(cdir)    
end