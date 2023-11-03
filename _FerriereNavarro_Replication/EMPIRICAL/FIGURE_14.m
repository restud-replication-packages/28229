clear; close all; clc

%% directories

cdir = pwd;
odir = ('output');
fdir = ('figures');

%% Load results

save_graph = 1;

cd(odir)
MULT_PROG  = readtable('LPM_FIGURE_14_pfed.csv','Delimiter',',');


m_P    = table2array(MULT_PROG(:,1));
m_se_P = table2array(MULT_PROG(:,2));
m_lb_P = table2array(MULT_PROG(:,3));
m_ub_P = table2array(MULT_PROG(:,4));
m_N    = table2array(MULT_PROG(:,5));
m_se_N = table2array(MULT_PROG(:,6));
m_lb_N = table2array(MULT_PROG(:,7));
m_ub_N = table2array(MULT_PROG(:,8));
pp_val = table2array(MULT_PROG(:,12));
HORIZ  = table2array(MULT_PROG(:,14));
cd(cdir)

H = max(HORIZ);

%%
xlb = 0; xub = 15;


fig = figure(13); clf;
fill([0:H, H:-1:0],[m_lb_P(1:H+1)', m_ub_P(H+1:-1:1)'],[0.40 0.60 1],'edgecolor',[0.40 0.60 1])
hold on
f1 = fill([0:H, H:-1:0],[m_lb_N(1:H+1)', m_ub_N(H+1:-1:1)'],[0.60 0.30 0.5],'edgecolor',[0.60 0.30 0.5],'FaceAlpha',.5,'EdgeAlpha',.5);
yy = plot(0:H,m_P(1:H+1),'b',0:H,m_N(1:H+1),'r','LineWidth',3);
hold on
    plot(0:H,zeros(H+1,1),'-.','Color',[.3 .3 .3],'LineWidth',3) 
hold off
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
ylabel('Cumulative Multiplier','Interpreter','LaTex','Fontsize',30)
leg = legend(yy,'Progressive','Non-Progressive');
set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northwest')
legend boxoff
xlim([xlb xub])
hold off
ylim([-0.7 1.1])

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )    
    cd(fdir)
    print('FIGURE_14','-dpng','-r0')
    % print('DATAP_Appendix_MULTxP_vAF','-dpng','-r0')
    cd(cdir)
end




