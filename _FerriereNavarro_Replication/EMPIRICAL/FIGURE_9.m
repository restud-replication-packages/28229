clear; close all; clc

%% directories

cdir = pwd;
odir = ('output');
fdir = ('figures');

%% Load data

cd(odir)
MULT_LINEAR  = readtable('LPM_FIGURE_9.csv','Delimiter',',');
L.m    = table2array(MULT_LINEAR(:,1));
L.m_se = table2array(MULT_LINEAR(:,2));
L.m_lb = table2array(MULT_LINEAR(:,3));
L.m_ub = table2array(MULT_LINEAR(:,4));
HORIZ  = table2array(MULT_LINEAR(:,6));
H      = max(HORIZ);

cd(cdir)

%% plot

save_graph = 1;

xlb = 0; xub = 15;
figure(13)

fill([0:H, H:-1:0],[L.m_lb(1:H+1)', L.m_ub(H+1:-1:1)'],[.4 .6 .6],'edgecolor',[.4 .6 .6])
hold on
yy = plot(0:H,L.m,'k','LineWidth',3.15);
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
ylabel('Cumulative Multiplier','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
% leg = legend(yy,'Linear');
% set(leg,'Interpreter','LaTex','Fontsize',17)
% legend boxoff
xlim([xlb xub])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_9','-dpng','-r0')
    cd(cdir)
end