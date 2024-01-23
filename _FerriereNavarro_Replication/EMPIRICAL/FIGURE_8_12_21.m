clear; close all; clc

%% directories

cdir = pwd;
odir = ('output');
fdir = ('figures');

%% load data

save_graph = 1;

DATA  = readtable('DATA_MACRO_FN.csv','Delimiter',',');
QUARTER = table2array(DATA(:,1));
RGOVPC  = table2array(DATA(:,3));
NEWS    = table2array(DATA(:,4));
psoc    = table2array(DATA(:,18));
pfed    = table2array(DATA(:,19));

%% FIGURE 8: G plot


cc = [0.0 0.45 0.74];
ck = [0.0 0.00 0.00];

ylb1 = 0;    yub1 = 0.0105;
ylb2 = -25;  yub2 = 100;

qd1 = 1914.50; qd2 = 1939.50; qd3 = 1950.50;
qd4 = 1965.00; qd5 = 1980.00; qd6 = 2001.50;
qdv = [qd1, qd2, qd3, qd4, qd5, qd6];


fig = figure(10); clf;
subplot(2,1,1)
for qq = 1:6
    plot([qdv(qq) qdv(qq)],[ylb1 yub1],'LineWidth',2.25,'color',ck); hold on
end
plot(QUARTER,RGOVPC   ,'LineWidth',3.25,'color',cc); hold on
set(gca,'XGrid','off','YGrid','on','Fontsize',21) 
title('Government Spending','Interpreter','LaTex','Fontsize',30)
ylabel('Millions of 2005 \$ ','Interpreter','LaTex','Fontsize',27)
xlabel('Year','Interpreter','LaTex','Fontsize',27)
xlim([1913 2012])
ylim([ylb1 yub1])

hold off

subplot(2,1,2)
for qq = 1:6
    plot([qdv(qq) qdv(qq)],[ylb2 yub2],'LineWidth',2.25,'color',ck); hold on
end
plot(QUARTER,100*NEWS   ,'LineWidth',3.25,'color',cc); hold on
set(gca,'XGrid','off','YGrid','on','Fontsize',21) 
title('News','Interpreter','LaTex','Fontsize',30)
ylabel('\% of GDP ','Interpreter','LaTex','Fontsize',27)
xlabel('Year','Interpreter','LaTex','Fontsize',27)
xlim([1913 2012])
ylim([ylb2 yub2])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_8','-dpng','-r0')
    cd(cdir)
end

%% tax data - annual data 

inq1 = find(QUARTER == floor(QUARTER));
YEAR = floor(QUARTER(inq1));
pfedY = pfed(inq1);
psocY = psoc(inq1);

%% FIGURE 12: progressiviy benchmark



fig = figure(10101); clf;
plot(YEAR,pfedY   ,'LineWidth',3.25,'color',[0.9 0.4 0.4]); hold on
set(gca,'XGrid','off','YGrid','on','Fontsize',21) 
% title('$P$ measures','Interpreter','LaTex','Fontsize',30)
ylabel('Progressivity $\gamma$','Interpreter','LaTex','Fontsize',27)
xlabel('Year','Interpreter','LaTex','Fontsize',27)
xlim([1913 2012])
hold off
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_12','-dpng','-r0')
    cd(cdir)
end

%% FIGURE 21: progressiviy benchmark + psoc

fig = figure(1000); clf;
plot(YEAR,pfedY ,'LineWidth',3.25,'color',[0.9 0.4 0.4]); hold on
plot(YEAR,psocY,'LineWidth',3.25,'color',[0.4 0.4 0.9]);
set(gca,'XGrid','off','YGrid','on','Fontsize',21) 
% title('$P$ measures','Interpreter','LaTex','Fontsize',30)
xlabel('Year','Interpreter','LaTex','Fontsize',27)
xlim([1913 2012])
leg = legend('Benchmark','w/ Social Security');
set(leg,'Interpreter','LaTex','Fontsize',27)
legend boxoff
hold off
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_21','-dpng','-r0')
    cd(cdir)
end