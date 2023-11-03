clear; close all; clc

%% directories

cdir  = pwd;
odir  = ('output');
RZdir = ('RZ_codes');
fdir  = ('figures');


%% Load data

cd(odir)
MULT_LINEAR  = readtable('LPM_FIGURE_26_benchmark.csv','Delimiter',',');
L.m    = table2array(MULT_LINEAR(:,1));
L.m_se = table2array(MULT_LINEAR(:,2));
L.m_lb = table2array(MULT_LINEAR(:,3));
L.m_ub = table2array(MULT_LINEAR(:,4));
HORIZ  = table2array(MULT_LINEAR(:,6));
H      = max(HORIZ);

%---no mtr results
MULT_LINEAR_nomtr  = readtable('LPM_FIGURE_26_no_mtr.csv','Delimiter',',');
L_nomtr.m    = table2array(MULT_LINEAR_nomtr(:,1));
L_nomtr.m_se = table2array(MULT_LINEAR_nomtr(:,2));
L_nomtr.m_lb = table2array(MULT_LINEAR_nomtr(:,3));
L_nomtr.m_ub = table2array(MULT_LINEAR_nomtr(:,4));
% HORIZ  = table2array(MULT_LINEAR(:,10));
% H      = max(HORIZ);
cd(cdir)

%---Ramey results
cd(RZdir)

MULT_RZ  = readtable('RZresults_FNsample.csv','Delimiter',',');
RZ.H    = table2array(MULT_RZ(:,1));
RZ.m    = table2array(MULT_RZ(:,2));
RZ.se   = table2array(MULT_RZ(:,3));
RZ.m_lb = RZ.m - RZ.se;
RZ.m_ub = RZ.m + RZ.se;


MULT_RZfs = readtable('RZresults_fullsample.csv','Delimiter',',');
RZfs.H    = table2array(MULT_RZfs(:,1));
RZfs.m    = table2array(MULT_RZfs(:,2));
RZfs.se   = table2array(MULT_RZfs(:,3));
RZfs.m_lb = RZfs.m - RZfs.se;
RZfs.m_ub = RZfs.m + RZfs.se;

cd(cdir)

zz = [L.m, L.m_se,L_nomtr.m,L_nomtr.m_se];

zzRZ = [RZ.m, RZ.se, RZfs.m, RZfs.se];


%% Plot w/ RZ results

save_graph = 1;

xlb = 0; xub = 15;

fig=figure(1301); clf;
fill([0:H, H:-1:0],[L.m_lb(1:H+1)', L.m_ub(H+1:-1:1)'],[.4 .6 .6],'edgecolor',[.4 .6 .6]); hold on
fill([0:H, H:-1:0],[L_nomtr.m_lb(1:H+1)', L_nomtr.m_ub(H+1:-1:1)'],[.4 .4 .4],'edgecolor',[.4 .4 .4],'FaceAlpha',.75,'EdgeAlpha',.75)
fill([1:H, H:-1:1],[RZ.m_lb(2:H+1)', RZ.m_ub(H+1:-1:2)'],[.8 .4 .4],'edgecolor',[.8 .4 .4],'FaceAlpha',.40,'EdgeAlpha',.40 ) 
fill([1:H, H:-1:1],[RZfs.m_lb(2:H+1)', RZfs.m_ub(H+1:-1:2)'],[.8 .4 .7],'edgecolor',[.8 .4 .7],'FaceAlpha',.40,'EdgeAlpha',.40 ) 
hold on
yy(1) = plot(0:H,L.m,'k','LineWidth',3.15);
yy(2) = plot(0:H,L_nomtr.m,'color',[0.9 0.5 0.5],'LineWidth',3.15);
yy(3) = plot(RZ.H(1:H+1),RZ.m(1:H+1),'color',[0.15 0.40 0.85],'LineWidth',3.15,'marker','o','markersize',13);
yy(4) = plot(RZfs.H(1:H+1),RZfs.m(1:H+1),'color',[0.45 0.40 0.85],'LineWidth',3.15,'marker','x','markersize',13);
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
ylabel('Cumulative Multiplier','Interpreter','LaTex','Fontsize',30)
xlabel('Quarters','Interpreter','LaTex','Fontsize',27)
leg = legend('Benchmark','no AMTR','Ramey-Zubairy for 1913-2006','Ramey-Zubairy for 1889-2012',...
             'Benchmark','no AMTR','Ramey-Zubairy for 1913-2006','Ramey-Zubairy for 1889-2012');
set(leg,'Interpreter','LaTex','Fontsize',23,'Location','SouthEast','NumColumns',2)
% legyy = legend(yy,'Benchmark','no AMTR','Ramey-Zubairy for 1913-2006','Ramey-Zubairy for 1889-2012');
% set(legyy,'Interpreter','LaTex','Fontsize',23,'Location','SouthEast')
legend boxoff
xlim([xlb xub])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_26','-dpng','-r0')
    cd(cdir)
end



% title('State-dependent Output Multipliers')

