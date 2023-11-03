clear; close all; clc

%% directories

cdir = pwd;
odir = ('output');
fdir = ('figures');

%% Load results

save_graph = 1;

adj = 1; % 1/100;

cd(odir)
%----------------------------------------------------------------------------
%---by state
IRF_tb3_allt  = readtable('LPM_FIGURE_24_1913_2006.csv','Delimiter',',');
HORIZ           = table2array(IRF_tb3_allt(:,14));
tb3_allt.m_P    = table2array(IRF_tb3_allt(:,1))*adj;
tb3_allt.m_se_P = table2array(IRF_tb3_allt(:,2))*adj;
tb3_allt.m_lb_P = table2array(IRF_tb3_allt(:,3))*adj;
tb3_allt.m_ub_P = table2array(IRF_tb3_allt(:,4))*adj;

tb3_allt.m_N    = table2array(IRF_tb3_allt(:,5))*adj;
tb3_allt.m_se_N = table2array(IRF_tb3_allt(:,6))*adj;
tb3_allt.m_lb_N = table2array(IRF_tb3_allt(:,7))*adj;
tb3_allt.m_ub_N = table2array(IRF_tb3_allt(:,8))*adj;

H               = max(HORIZ);

IRF_tb3_post51  = readtable('LPM_FIGURE_24_1951_2006.csv','Delimiter',',');
tb3_post51.m_P    = table2array(IRF_tb3_post51(:,1))*adj;
tb3_post51.m_se_P = table2array(IRF_tb3_post51(:,2))*adj;
tb3_post51.m_lb_P = table2array(IRF_tb3_post51(:,3))*adj;
tb3_post51.m_ub_P = table2array(IRF_tb3_post51(:,4))*adj;

tb3_post51.m_N    = table2array(IRF_tb3_post51(:,5))*adj;
tb3_post51.m_se_N = table2array(IRF_tb3_post51(:,6))*adj;
tb3_post51.m_lb_N = table2array(IRF_tb3_post51(:,7))*adj;
tb3_post51.m_ub_N = table2array(IRF_tb3_post51(:,8))*adj;
%----------------------------------------------------------------------------
cd(cdir)

%% computations

tb3_allt.m_lbx_P   = tb3_allt.m_P   - 2.0*tb3_allt.m_se_P  ;  tb3_allt.m_ubx_P   = tb3_allt.m_P   + 2.0*tb3_allt.m_se_P;
tb3_allt.m_lbx_N   = tb3_allt.m_N   - 2.0*tb3_allt.m_se_N  ;  tb3_allt.m_ubx_N   = tb3_allt.m_N   + 2.0*tb3_allt.m_se_N;

tb3_post51.m_lbx_P = tb3_post51.m_P - 2.0*tb3_post51.m_se_P;  tb3_post51.m_ubx_P = tb3_post51.m_P + 2.0*tb3_post51.m_se_P;
tb3_post51.m_lbx_N = tb3_post51.m_N - 2.0*tb3_post51.m_se_N;  tb3_post51.m_ubx_N = tb3_post51.m_N + 2.0*tb3_post51.m_se_N;

%messi

%% Plot

colorLci  = [0.40 0.60 0.60];
colorPci  = [0.40 0.60 1.00];
colorNci  = [0.60 0.30 0.50];

colorL   = 'k';
colorP   = 'b';
colorN   = 'r';

czero   = [0.30 0.30 0.30];

xlb = 0; xub = 11;
ylb = -0.35; yub = 0.50;



%----------------------------------------------------------------------------
%---state
fig = figure(100101); clf;
subplot(1,2,1)
fill([0:H, H:-1:0],[tb3_allt.m_lb_P(1:H+1)', tb3_allt.m_ub_P(H+1:-1:1)'],colorPci,'edgecolor',colorPci); hold on
fill([0:H, H:-1:0],[tb3_allt.m_lb_N(1:H+1)', tb3_allt.m_ub_N(H+1:-1:1)'],colorNci,'edgecolor',colorNci,'FaceAlpha',.5,'EdgeAlpha',.5); % hold on
% f1 = fill([0:H, H:-1:0],[m_lb_N(1:H+1)', m_ub_N(H+1:-1:1)'],[0.60 0.30 0.5],'edgecolor',[0.60 0.30 0.5],'FaceAlpha',.5,'EdgeAlpha',.5);
yy = plot(0:H,tb3_allt.m_P(1:H+1),'b',0:H,tb3_allt.m_N(1:H+1),'r','LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('All sample','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('percentage points','Interpreter','LaTex','Fontsize',27)
leg = legend(yy,'Progressive','Non-Progressive');
set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northwest')
legend boxoff
xlim([xlb xub])
ylim([ylb yub])
hold off

subplot(1,2,2)
fill([0:H, H:-1:0],[tb3_post51.m_lb_P(1:H+1)', tb3_post51.m_ub_P(H+1:-1:1)'],colorPci,'edgecolor',colorPci); hold on
fill([0:H, H:-1:0],[tb3_post51.m_lb_N(1:H+1)', tb3_post51.m_ub_N(H+1:-1:1)'],colorNci,'edgecolor',colorNci,'FaceAlpha',.5,'EdgeAlpha',.5); % hold on
% f1 = fill([0:H, H:-1:0],[m_lb_N(1:H+1)', m_ub_N(H+1:-1:1)'],[0.60 0.30 0.5],'edgecolor',[0.60 0.30 0.5],'FaceAlpha',.5,'EdgeAlpha',.5);
yy = plot(0:H,tb3_post51.m_P(1:H+1),'b',0:H,tb3_post51.m_N(1:H+1),'r','LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Post Fed-Treasury Accord','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('percentage points','Interpreter','LaTex','Fontsize',27)
leg = legend(yy,'Progressive','Non-Progressive');
set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northwest')
legend boxoff
xlim([xlb xub])
ylim([ylb yub])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_24','-dpng','-r0')    
    cd(cdir)
end
%----------------------------------------------------------------------------