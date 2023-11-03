clear; close all; clc;

%% load data

cdir = pwd;
odir = ('output');
fdir = ('figures');

cd(odir)
%---1) lngov
DATA = readtable('LPM_FIGURE_28_29_lngov.csv','Delimiter',',');
IRF_g_m_P    = table2array(DATA(:,1));
IRF_g_m_se_P = table2array(DATA(:,2));
IRF_g_m_lb_P = table2array(DATA(:,3));
IRF_g_m_ub_P = table2array(DATA(:,4));
IRF_g_m_N    = table2array(DATA(:,5));
IRF_g_m_se_N = table2array(DATA(:,6));
IRF_g_m_lb_N = table2array(DATA(:,7));
IRF_g_m_ub_N = table2array(DATA(:,8));

%---2) tb3
DATA = readtable('LPM_FIGURE_28_29_tb3.csv','Delimiter',',');
IRF_tb3_m_P    = table2array(DATA(:,1));
IRF_tb3_m_se_P = table2array(DATA(:,2));
IRF_tb3_m_lb_P = table2array(DATA(:,3));
IRF_tb3_m_ub_P = table2array(DATA(:,4));
IRF_tb3_m_N    = table2array(DATA(:,5));
IRF_tb3_m_se_N = table2array(DATA(:,6));
IRF_tb3_m_lb_N = table2array(DATA(:,7));
IRF_tb3_m_ub_N = table2array(DATA(:,8));

%---3) rdef
DATA = readtable('LPM_FIGURE_28_29_rdef.csv','Delimiter',',');
IRF_rdef_m_P    = table2array(DATA(:,1));
IRF_rdef_m_se_P = table2array(DATA(:,2));
IRF_rdef_m_lb_P = table2array(DATA(:,3));
IRF_rdef_m_ub_P = table2array(DATA(:,4));
IRF_rdef_m_N    = table2array(DATA(:,5));
IRF_rdef_m_se_N = table2array(DATA(:,6));
IRF_rdef_m_lb_N = table2array(DATA(:,7));
IRF_rdef_m_ub_N = table2array(DATA(:,8));

%---4) lninv
DATA = readtable('LPM_FIGURE_28_29_lninv.csv','Delimiter',',');
IRF_inv_m_P    = table2array(DATA(:,1));
IRF_inv_m_se_P = table2array(DATA(:,2));
IRF_inv_m_lb_P = table2array(DATA(:,3));
IRF_inv_m_ub_P = table2array(DATA(:,4));
IRF_inv_m_N    = table2array(DATA(:,5));
IRF_inv_m_se_N = table2array(DATA(:,6));
IRF_inv_m_lb_N = table2array(DATA(:,7));
IRF_inv_m_ub_N = table2array(DATA(:,8));

%---5) lnwgenf
DATA = readtable('LPM_FIGURE_28_29_lnwgenf.csv','Delimiter',',');
IRF_lnwge_m_P    = table2array(DATA(:,1));
IRF_lnwge_m_se_P = table2array(DATA(:,2));
IRF_lnwge_m_lb_P = table2array(DATA(:,3));
IRF_lnwge_m_ub_P = table2array(DATA(:,4));
IRF_lnwge_m_N    = table2array(DATA(:,5));
IRF_lnwge_m_se_N = table2array(DATA(:,6));
IRF_lnwge_m_lb_N = table2array(DATA(:,7));
IRF_lnwge_m_ub_N = table2array(DATA(:,8));

%---6) lnhours
DATA = readtable('LPM_FIGURE_28_29_lnhours.csv','Delimiter',',');
IRF_lnh_m_P    = table2array(DATA(:,1));
IRF_lnh_m_se_P = table2array(DATA(:,2));
IRF_lnh_m_lb_P = table2array(DATA(:,3));
IRF_lnh_m_ub_P = table2array(DATA(:,4));
IRF_lnh_m_N    = table2array(DATA(:,5));
IRF_lnh_m_se_N = table2array(DATA(:,6));
IRF_lnh_m_lb_N = table2array(DATA(:,7));
IRF_lnh_m_ub_N = table2array(DATA(:,8));

cd(cdir)


HORIZ = table2array(DATA(:,14)) ;
H     = max(HORIZ);

%% plot

save_graph = 1;


colorLci  = [0.40 0.60 0.60];
colorPci  = [0.40 0.60 1.00];
colorNci  = [0.60 0.30 0.50];

colorL   = 'k';
colorP   = 'b';
colorN   = 'r';

czero   = [0.30 0.30 0.30];

xlb = 0; xub = 11;

fig = figure(1301); clf;
subplot(2,2,1)
fill([0:H, H:-1:0],[100*IRF_g_m_lb_P(1:H+1)', 100*IRF_g_m_ub_P(H+1:-1:1)'],colorPci,'edgecolor',colorPci); hold on
fill([0:H, H:-1:0],[100*IRF_g_m_lb_N(1:H+1)', 100*IRF_g_m_ub_N(H+1:-1:1)'],colorNci,'edgecolor',colorNci,'FaceAlpha',.5,'EdgeAlpha',.5); % hold on
yy = plot(0:H,100*IRF_g_m_P(1:H+1),colorP,0:H,100*IRF_g_m_N(1:H+1),colorN,'LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Spending','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('\%','Interpreter','LaTex','Fontsize',27)
leg = legend(yy,'Progressive','Non-Progressive');
set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northeast')
legend boxoff
xlim([xlb xub])
%ylim([ylb yub])
hold off
subplot(2,2,2)
fill([0:H, H:-1:0],[IRF_tb3_m_lb_P(1:H+1)', IRF_tb3_m_ub_P(H+1:-1:1)'],colorPci,'edgecolor',colorPci); hold on
fill([0:H, H:-1:0],[IRF_tb3_m_lb_N(1:H+1)', IRF_tb3_m_ub_N(H+1:-1:1)'],colorNci,'edgecolor',colorNci,'FaceAlpha',.5,'EdgeAlpha',.5); % hold on
yy = plot(0:H,IRF_tb3_m_P(1:H+1),colorP,0:H,IRF_tb3_m_N(1:H+1),colorN,'LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('TB3','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('percentage points','Interpreter','LaTex','Fontsize',27)
xlim([xlb xub])
%ylim([ylb yub])
hold off
subplot(2,2,3)
fill([0:H, H:-1:0],[100*IRF_rdef_m_lb_P(1:H+1)', 100*IRF_rdef_m_ub_P(H+1:-1:1)'],colorPci,'edgecolor',colorPci); hold on
fill([0:H, H:-1:0],[100*IRF_rdef_m_lb_N(1:H+1)', 100*IRF_rdef_m_ub_N(H+1:-1:1)'],colorNci,'edgecolor',colorNci,'FaceAlpha',.5,'EdgeAlpha',.5); % hold on
yy = plot(0:H,100*IRF_rdef_m_P(1:H+1),colorP,0:H,100*IRF_rdef_m_N(1:H+1),colorN,'LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Deficit','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('\%','Interpreter','LaTex','Fontsize',27)
xlim([xlb xub])
%ylim([ylb yub])

subplot(2,2,4)
fill([0:H, H:-1:0],[100*IRF_inv_m_lb_P(1:H+1)', 100*IRF_inv_m_ub_P(H+1:-1:1)'],colorPci,'edgecolor',colorPci); hold on
fill([0:H, H:-1:0],[100*IRF_inv_m_lb_N(1:H+1)', 100*IRF_inv_m_ub_N(H+1:-1:1)'],colorNci,'edgecolor',colorNci,'FaceAlpha',.5,'EdgeAlpha',.5); % hold on
yy = plot(0:H,100*IRF_inv_m_P(1:H+1),colorP,0:H,100*IRF_inv_m_N(1:H+1),colorN,'LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Investment','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('\%','Interpreter','LaTex','Fontsize',27)
%leg = legend(yy,'Progressive','Non-Progressive');
%set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northwest')
%legend boxoff
xlim([xlb xub])
%ylim([ylb yub])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_28','-dpng','-r0')    
    cd(cdir)
end



fig = figure(1302); clf;
subplot(1,2,1)
fill([0:H, H:-1:0],[100*IRF_lnwge_m_lb_P(1:H+1)', 100*IRF_lnwge_m_ub_P(H+1:-1:1)'],colorPci,'edgecolor',colorPci); hold on
fill([0:H, H:-1:0],[100*IRF_lnwge_m_lb_N(1:H+1)', 100*IRF_lnwge_m_ub_N(H+1:-1:1)'],colorNci,'edgecolor',colorNci,'FaceAlpha',.5,'EdgeAlpha',.5); % hold on
yy = plot(0:H,100*IRF_lnwge_m_P(1:H+1),colorP,0:H,100*IRF_lnwge_m_N(1:H+1),colorN,'LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Wages','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('\%','Interpreter','LaTex','Fontsize',27)
xlim([xlb xub])
%ylim([ylb yub])
leg = legend(yy,'Progressive','Non-Progressive');
set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northwest')
legend boxoff
hold off

subplot(1,2,2)
fill([0:H, H:-1:0],[100*IRF_lnh_m_lb_P(1:H+1)', 100*IRF_lnh_m_ub_P(H+1:-1:1)'],colorPci,'edgecolor',colorPci); hold on
fill([0:H, H:-1:0],[100*IRF_lnh_m_lb_N(1:H+1)', 100*IRF_lnh_m_ub_N(H+1:-1:1)'],colorNci,'edgecolor',colorNci,'FaceAlpha',.5,'EdgeAlpha',.5); % hold on
yy = plot(0:H,100*IRF_lnh_m_P(1:H+1),colorP,0:H,100*IRF_lnh_m_N(1:H+1),colorN,'LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Hours','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('\%','Interpreter','LaTex','Fontsize',27)

xlim([xlb xub])
%ylim([ylb yub])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_29','-dpng','-r0')    
    cd(cdir)
end
