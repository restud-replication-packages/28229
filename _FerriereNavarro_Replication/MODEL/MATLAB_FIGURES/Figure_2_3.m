clear; close all; clc;

%% set up

cdir = pwd;

savef = 1;  % =1 to save figures

CGdir = '../CODES/BENCHMARK/CG/OUTPUT'; 
MGdir = '../CODES/BENCHMARK/MG/OUTPUT';




%% Load CG results

cd(CGdir)
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]

lbd  = AGGREGATES(1); Gagg = AGGREGATES(2); Cagg = AGGREGATES(3);
Yagg = AGGREGATES(4); Lagg = AGGREGATES(5); EMP  = AGGREGATES(6);
Kagg = AGGREGATES(7); DB   = AGGREGATES(8); Iagg = AGGREGATES(9);


%%% Transition
load lbd_TR.txt;  load gma_TR.txt;  load wge_TR.txt;
load wgeH_TR.txt; load Lagg_TR.txt; load Yagg_TR.txt;
load Cagg_TR.txt; load G_TR.txt; 
cd(cdir)

T_TR = numel(lbd_TR);

%% Computations CG results

CG.mult = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-Gagg);

CG.IRFG  = 100*log(G_TR/Gagg);
CG.IRFL  = 100*log(Lagg_TR/Lagg);
CG.IRFC  = 100*log(Cagg_TR/Cagg);
CG.IRFYd = 100*log((G_TR+Cagg_TR)/(Gagg+Cagg));
CG.IRFY  = 100*log(Yagg_TR/Yagg);

%% Load MG results

cd(MGdir)
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]

lbd  = AGGREGATES(1); Gagg = AGGREGATES(2); Cagg = AGGREGATES(3);
Yagg = AGGREGATES(4); Lagg = AGGREGATES(5); EMP  = AGGREGATES(6);
Kagg = AGGREGATES(7); DB   = AGGREGATES(8); Iagg = AGGREGATES(9);



%%% Transition
load lbd_TR.txt; load gma_TR.txt; load wge_TR.txt;
load wgeH_TR.txt; load Lagg_TR.txt; load Yagg_TR.txt;
load Cagg_TR.txt; load Iagg_TR.txt; load Aagg_TR.txt;
load qk_TR.txt; load rk_TR.txt;
load G_TR.txt; load div_TR.txt; load Kagg_TR.txt;
load DB_TR.txt; load div_TR.txt;
cd(cdir)


%% Computations MG results

MG.mult = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-Gagg);

MG.IRFG  = 100*log(G_TR/Gagg);
MG.IRFL  = 100*log(Lagg_TR/Lagg);
MG.IRFC  = 100*log(Cagg_TR/Cagg);
MG.IRFYd = 100*log((G_TR+Cagg_TR)/(Gagg+Cagg));
MG.IRFY  = 100*log(Yagg_TR/Yagg);

%% Plot IRF

colorCG = [1.00 0.20 0.30];
colorMG = [0.30 0.40 0.90];

tlb = 0; tub = 15;
fig = figure(1001); clf;

subplot(1,3,1)
plot(0:T_TR-1,CG.IRFY,'color',colorCG,'LineWidth',3.25,'Marker','o','MarkerSize',4.25); hold on
plot(0:T_TR-1,MG.IRFY,'color',colorMG,'LineWidth',3.25);
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Output','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
ylim([-0.005 0.0375])
ylabel('$\%$ from steady-state','Interpreter','LaTex','Fontsize',27)
hold off
leg = legend('Constant Progressivity','Higher Progressivity');
set(leg,'Interpreter','LaTex','Fontsize',25,'Location','NorthEast')
legend boxoff

subplot(1,3,2)
plot(0:T_TR-1,CG.IRFL,'color',colorCG,'LineWidth',3.25,'Marker','o','MarkerSize',4.25); hold on
plot(0:T_TR-1,MG.IRFL,'color',colorMG,'LineWidth',3.25);
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Labor','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
% leg = legend('Constant Progressivity','Higher Progressivity');
% set(leg,'Interpreter','LaTex','Fontsize',25,'Location','SouthEast')
% legend boxoff
xlim([tlb tub])
ylim([-0.01 0.055])
hold off

subplot(1,3,3)
plot(0:T_TR-1,CG.IRFC,'color',colorCG,'LineWidth',3.25,'Marker','o','MarkerSize',4.25); hold on
plot(0:T_TR-1,MG.IRFC,'color',colorMG,'LineWidth',3.25);
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Consumption','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
% leg = legend('Constant Progressivity','Higher Progressivity');
% set(leg,'Interpreter','LaTex','Fontsize',25,'Location','SouthEast')
% legend boxoff
xlim([tlb tub])
ylim([-0.16 0.01])
hold off



if savef == 1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 22 11];
    cd ../FIGURES
    print('FIGURE_2','-dpng','-r0')
    cd(cdir)
end


%% Plot multiplier


tindex = [4; 8; 40];
fig = figure(102); clf;
% b = bar(tindex/4,[CG.mult(tindex) MG.mult(tindex)]);
b = bar([1 2 3],[CG.mult(tindex) MG.mult(tindex)]);
b(1).FaceColor = colorCG;
b(2).FaceColor = colorMG;
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
% title('Cumulative Output Multipliers','Interpreter','LaTex','Fontsize',27) % $\frac{\sum_{t}{Y_t - \bar{Y}}}{\sum_{t}{G_t - \bar{G}}}$
xlabel('Year','Interpreter','LaTex','Fontsize',27)
ylabel('Multiplier','Interpreter','LaTex','Fontsize',27)
set(gca, 'XTickLabel', {'1' '2','10'})
ylim([-0.1 0.4])
leg = legend('Constant Progressivity','Higher Progressivity','Location','SouthEast');
set(leg,'Interpreter','LaTex','Fontsize',25)
legend boxoff

if savef == 1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 22 11];
    cd ../FIGURES
    print('FIGURE_3','-dpng','-r0')
    cd(cdir)
end


