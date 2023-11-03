clear; close all; clc;

cdir = pwd;
odir = ('output');
fdir = ('figures');  


%% Load results

save_graph = 1;

cd(odir)
%--load ATR LPM result
IRF_ATR  = readtable('LPM_FIGURE_10_atr_psz.csv','Delimiter',',');
ATR.m    = table2array(IRF_ATR(:,1));
ATR.m_se = table2array(IRF_ATR(:,2));
ATR.m_lb = table2array(IRF_ATR(:,3));
ATR.m_ub = table2array(IRF_ATR(:,4));
HORIZ    = table2array(IRF_ATR(:,6));
H        = max(HORIZ);


%--load panel result
%IRF_ATRS  = readtable('LPM_FIGURE_25_bp_fullsample.csv','Delimiter',',');
IRF_ATRS  = readtable('LPM_FIGURE_25.csv','Delimiter',',');
ATRS.m    = table2array(IRF_ATRS(:,1));
ATRS.m_se = table2array(IRF_ATRS(:,2));
ATRS.m_lb = table2array(IRF_ATRS(:,3));
ATRS.m_ub = table2array(IRF_ATRS(:,4));
% HORIZ     = table2array(IRF_ATRS(:,5));
% H         = max(HORIZ);
cd(cdir)

%% Plot

xlb = 0;       xub = 15;
ylb = -0.010;  yub = 0.025;
cATR    = [0.50 0.70 0.60]; cATRa2    = [0.40, 0.60, 0.60]; 
cATRS   = [0.90 0.10 0.30]; cATRSa2   = [0.70, 0.30, 0.10]; 
cT50    = [0.80 0.50 0.50]; cT50a2    = [0.60, 0.40, 0.40]; 
cB50    = [0.50 0.50 0.80]; cB50a2    = [0.40, 0.40, 0.60]; 
cT50B50 = [0.60 0.60 0.90]; cT50B50a2 = [0.60, 0.90, 0.94]; 
czero   = [0.30 0.30 0.30];


fig = figure(1001); clf;
subplot(1,2,1)
fill([0:H, H:-1:0],[ATR.m_lb(1:H+1)', ATR.m_ub(H+1:-1:1)'],cATRa2,'edgecolor',cATRa2); hold on
plot(HORIZ,ATR.m,'color',cATR,'LineWidth',3.25); 
plot(HORIZ,zeros(1,H+1),'color',czero,'LineWidth',2.25)
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('Average Tax Rate - Federal','Interpreter','LaTex','Fontsize',30)
ylabel('percentage points difference','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',25)
xlim([xlb xub])
ylim([ylb yub])
hold off

subplot(1,2,2)
fill([0:H, H:-1:0],[ATRS.m_lb(1:H+1)', ATRS.m_ub(H+1:-1:1)'],cATRSa2,'edgecolor',cATRa2); hold on
plot(HORIZ,ATRS.m,'color',cATRS,'LineWidth',3.25); 
plot(HORIZ,zeros(1,H+1),'color',czero,'LineWidth',2.25)
set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
title('Average Tax Rate - State','Interpreter','LaTex','Fontsize',30)
% ylabel('percentage points difference','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',25)
xlim([xlb xub])
ylim([ylb yub])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_25','-dpng','-r0')
    cd(cdir)
end

