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
%---linear
IRF_tb3_allt_LINEAR  = readtable('LPM_FIGURE_23_1913_2006.csv','Delimiter',',');
tb3_allt.m_LINEAR    = table2array(IRF_tb3_allt_LINEAR(:,1))*adj;
tb3_allt.m_se_LINEAR = table2array(IRF_tb3_allt_LINEAR(:,2))*adj;
tb3_allt.m_lb_LINEAR = table2array(IRF_tb3_allt_LINEAR(:,3))*adj;
tb3_allt.m_ub_LINEAR = table2array(IRF_tb3_allt_LINEAR(:,4))*adj;
HORIZ                = table2array(IRF_tb3_allt_LINEAR(:,6));
H = max(HORIZ);

IRF_tb3_post51_LINEAR  = readtable('LPM_FIGURE_23_1951_2006.csv','Delimiter',',');
tb3_post51.m_LINEAR    = table2array(IRF_tb3_post51_LINEAR(:,1))*adj;
tb3_post51.m_se_LINEAR = table2array(IRF_tb3_post51_LINEAR(:,2))*adj;
tb3_post51.m_lb_LINEAR = table2array(IRF_tb3_post51_LINEAR(:,3))*adj;
tb3_post51.m_ub_LINEAR = table2array(IRF_tb3_post51_LINEAR(:,4))*adj;
% HORIZ           = table2array(IRF_tb3_allt_LINEAR(:,5));

IRF_tb3_post80_LINEAR  = readtable('LPM_FIGURE_23_1980_2006.csv','Delimiter',',');
tb3_post80.m_LINEAR    = table2array(IRF_tb3_post80_LINEAR(:,1))*adj;
tb3_post80.m_se_LINEAR = table2array(IRF_tb3_post80_LINEAR(:,2))*adj;
tb3_post80.m_lb_LINEAR = table2array(IRF_tb3_post80_LINEAR(:,3))*adj;
tb3_post80.m_ub_LINEAR = table2array(IRF_tb3_post80_LINEAR(:,4))*adj;

%----------------------------------------------------------------------------
cd(cdir)


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
%---linear
fig = figure(100100); clf
subplot(1,3,1)
fill([0:H, H:-1:0],[tb3_allt.m_lb_LINEAR(1:H+1)', tb3_allt.m_ub_LINEAR(H+1:-1:1)'],colorLci,'edgecolor',colorLci); hold on
yy = plot(0:H,tb3_allt.m_LINEAR(1:H+1),'k','LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('All sample','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('percentage points','Interpreter','LaTex','Fontsize',27)
% leg = legend(yy,'Progressive','Non-Progressive');
% set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northwest')
% legend boxoff
xlim([xlb xub])
ylim([ylb yub])
hold off

subplot(1,3,2)
fill([0:H, H:-1:0],[tb3_post51.m_lb_LINEAR(1:H+1)', tb3_post51.m_ub_LINEAR(H+1:-1:1)'],colorLci,'edgecolor',colorPci); hold on
yy = plot(0:H,tb3_post51.m_LINEAR(1:H+1),'k','LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Post Fed-Treasury Accord','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('percentage points','Interpreter','LaTex','Fontsize',27)
% leg = legend(yy,'Progressive','Non-Progressive');
% set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northwest')
% legend boxoff
xlim([xlb xub])
ylim([ylb yub])
hold off

subplot(1,3,3)
fill([0:H, H:-1:0],[tb3_post80.m_lb_LINEAR(1:H+1)', tb3_post80.m_ub_LINEAR(H+1:-1:1)'],colorLci,'edgecolor',colorPci); hold on
yy = plot(0:H,tb3_post80.m_LINEAR(1:H+1),'k','LineWidth',3);
plot(0:H,zeros(H+1,1),'-.','Color',czero,'LineWidth',3) 
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Post 1980','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('percentage points','Interpreter','LaTex','Fontsize',27)
% leg = legend(yy,'Progressive','Non-Progressive');
% set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northwest')
% legend boxoff
xlim([xlb xub])
ylim([ylb yub])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_23','-dpng','-r0')    
    cd(cdir)
end


% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 22 11];
% if (save_graph == 1 )
%     print('Robustness_with_levels_env','-dpng','-r0')
% end
%----------------------------------------------------------------------------