clear; close all; clc;

%% directories

cdir = pwd;
odir = ('output');
fdir = ('figures');

%% Load results

cd(odir)
IRF_ATR  = readtable('LPM_FIGURE_10_atr_psz.csv','Delimiter',',');
ATR.m    = table2array(IRF_ATR(:,1));
ATR.m_se = table2array(IRF_ATR(:,2));
ATR.m_lb = table2array(IRF_ATR(:,3));
ATR.m_ub = table2array(IRF_ATR(:,4));
HORIZ    = table2array(IRF_ATR(:,6));
H        = max(HORIZ);

IRF_T50  = readtable('LPM_FIGURE_10_atr_t50_psz.csv','Delimiter',',');
T50.m    = table2array(IRF_T50(:,1));
T50.m_se = table2array(IRF_T50(:,2));
T50.m_lb = table2array(IRF_T50(:,3));
T50.m_ub = table2array(IRF_T50(:,4));

IRF_B50  = readtable('LPM_FIGURE_10_atr_b50_psz.csv','Delimiter',',');
B50.m    = table2array(IRF_B50(:,1));
B50.m_se = table2array(IRF_B50(:,2));
B50.m_lb = table2array(IRF_B50(:,3));
B50.m_ub = table2array(IRF_B50(:,4));

IRF_T50B50  = readtable('LPM_FIGURE_10_datr_t50b50.csv','Delimiter',',');
T50B50.m    = table2array(IRF_T50B50(:,1));
T50B50.m_se = table2array(IRF_T50B50(:,2));
T50B50.m_lb = table2array(IRF_T50B50(:,3));
T50B50.m_ub = table2array(IRF_T50B50(:,4));

cd(cdir)


zz = [ATR.m, T50.m, B50.m, T50B50.m];

%% Plot

save_graph = 1;

xlb = 0; xub = 15;
ylb = -0.015;  yub = 0.055;
cATR    = [0.50 0.70 0.60]; cATRa2    = [0.40, 0.60, 0.60]; 
cT50    = [0.80 0.50 0.50]; cT50a2    = [0.60, 0.40, 0.40]; 
cB50    = [0.50 0.50 0.80]; cB50a2    = [0.40, 0.40, 0.60]; 
cT50B50 = [0.60 0.60 0.90]; cT50B50a2 = [0.60, 0.90, 0.94]; 
czero   = [0.30 0.30 0.30];

figure(1001)
subplot(1,3,1)
fill([0:H, H:-1:0],[ATR.m_lb(1:H+1)', ATR.m_ub(H+1:-1:1)'],cATRa2,'edgecolor',cATRa2); hold on
plot(HORIZ,ATR.m,'color',cATR,'LineWidth',3.25); 
% plot(HORIZ,MODEL.ATR(1:H+1),'color',cATR,'LineWidth',3.25); 
plot(HORIZ,zeros(1,H+1),'color',czero,'LineWidth',2.25)
title('Average Tax Rate','Interpreter','LaTex','Fontsize',30)
ylabel('percentage points difference','Interpreter','LaTex','Fontsize',27)
xlabel('Quarters','Interpreter','LaTex','Fontsize',25)
set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
xlim([xlb xub])
ylim([ylb yub])
hold off

subplot(1,3,2)
fill([0:H, H:-1:0],[T50.m_lb(1:H+1)', T50.m_ub(H+1:-1:1)'],cT50a2,'edgecolor',cT50a2); hold on
fill([0:H, H:-1:0],[B50.m_lb(1:H+1)', B50.m_ub(H+1:-1:1)'],cB50a2,'edgecolor',cB50a2); 
plot(HORIZ,T50.m,'color',cT50,'LineWidth',3.25); 
plot(HORIZ,B50.m,'color',cB50,'LineWidth',3.25); 
% plot(HORIZ,MODEL.T50(1:H+1),'color',cT50,'LineWidth',3.25); 
% plot(HORIZ,MODEL.B50(1:H+1),'color',cB50,'LineWidth',3.25); plot(HORIZ,zeros(1,H+1),'color',czero,'LineWidth',2.25)
plot(HORIZ,zeros(1,H+1),'color',czero,'LineWidth',2.25)
title('ATR by Group','Interpreter','LaTex','Fontsize',30)
set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
xlabel('Quarters','Interpreter','LaTex','Fontsize',25)
hold off
leg = legend('Top $50\%$','Bottom $50\%$');
set(leg,'Interpreter','LaTex','Fontsize',23)
xlim([xlb xub])
ylim([ylb yub])
legend boxoff

subplot(1,3,3)
fill([0:H, H:-1:0],[T50B50.m_lb(1:H+1)', T50B50.m_ub(H+1:-1:1)'],cT50B50a2,'edgecolor',cT50B50a2); hold on
plot(HORIZ,T50B50.m,'color',cT50B50,'LineWidth',3.25); 
% plot(HORIZ,MODEL.T50(1:H+1)-MODEL.B50(1:H+1),'color',cT50B50,'LineWidth',3.25); 
plot(HORIZ,zeros(1,H+1),'color',czero,'LineWidth',2.25)
title('ATR top $50\%$ - ATR bottom $50\%$','Interpreter','LaTex','Fontsize',30)
set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
xlabel('Quarters','Interpreter','LaTex','Fontsize',25)
xlim([xlb xub])
ylim([ylb yub])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_10','-dpng','-r0')
    cd(cdir)
end


%% Plot

xlb = 0; xub = 15;

cATR    = [0.50 0.70 0.60]; cATRa2    = [0.40, 0.60, 0.60]; 
cT50    = [0.80 0.50 0.50]; cT50a2    = [0.60, 0.40, 0.40]; 
cB50    = [0.50 0.50 0.80]; cB50a2    = [0.40, 0.40, 0.60]; 
cT50B50 = [0.60 0.60 0.90]; cT50B50a2 = [0.60, 0.90, 0.94]; 
czero   = [0.30 0.30 0.30];

figure(1002)
subplot(1,2,1)
fill([0:H, H:-1:0],[ATR.m_lb(1:H+1)', ATR.m_ub(H+1:-1:1)'],cATRa2,'edgecolor',cATRa2); hold on
plot(HORIZ,ATR.m,'color',cATR,'LineWidth',3.25); 
% plot(HORIZ,MODEL.ATR(1:H+1),'color',cATR,'LineWidth',3.25); 
plot(HORIZ,zeros(1,H+1),'color',czero,'LineWidth',2.25)
title('Average Tax Rate','Interpreter','LaTex','Fontsize',30)
ylabel('percentage points difference','Interpreter','LaTex','Fontsize',27)
xlabel('Quarters','Interpreter','LaTex','Fontsize',25)
set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
xlim([xlb xub])
hold off

subplot(1,2,2)
fill([0:H, H:-1:0],[T50.m_lb(1:H+1)', T50.m_ub(H+1:-1:1)'],cT50a2,'edgecolor',cT50a2); hold on
fill([0:H, H:-1:0],[B50.m_lb(1:H+1)', B50.m_ub(H+1:-1:1)'],cB50a2,'edgecolor',cB50a2); 
plot(HORIZ,T50.m,'color',cT50,'LineWidth',3.25); 
plot(HORIZ,B50.m,'color',cB50,'LineWidth',3.25); 
% plot(HORIZ,MODEL.T50(1:H+1),'color',cT50,'LineWidth',3.25); 
% plot(HORIZ,MODEL.B50(1:H+1),'color',cB50,'LineWidth',3.25); plot(HORIZ,zeros(1,H+1),'color',czero,'LineWidth',2.25)
plot(HORIZ,zeros(1,H+1),'color',czero,'LineWidth',2.25)
title('ATR by Group','Interpreter','LaTex','Fontsize',30)
set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
xlabel('Quarters','Interpreter','LaTex','Fontsize',25)
hold off
leg = legend('Top $50\%$','Bottom $50\%$');
set(leg,'Interpreter','LaTex','Fontsize',23)
xlim([xlb xub])
legend boxoff

