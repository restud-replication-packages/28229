clear;  close all; clc;

%% load results
save_graph = 1; 

cdir  = pwd;  % current directory

vec = [1; 4; 8; 12; 16; 40];

%% Debt 

mult_mat_div = ones(6,3,3); % 6 periods, CG-MG-DELTA, low-benchmark-high

%% Benchmark

mult_mat_debt = ones(6,3,3); % 6 periods, CG-MG-DELTA, **-benchmark-high

cd '../CODES/BENCHMARK/CG/OUTPUT'
load AGGREGATES.txt; G    = AGGREGATES(2); Yagg = AGGREGATES(4);
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)


mult_aux = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

mult_mat_div(:,1,1) = mult_aux(vec);

cd '../CODES/BENCHMARK/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)

mult_aux = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

mult_mat_div(:,2,1) = mult_aux(vec);
mult_mat_div(:,3,1) = mult_mat_div(:,2,1) - mult_mat_div(:,1,1);


%% Profits LESS CONCENTRATED

cd '../CODES/ROBUST_APPEND/PROFITLESS/CG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)

mult_aux = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

mult_mat_div(:,1,3) = mult_aux(vec);

cd '../CODES/ROBUST_APPEND/PROFITLESS/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)

mult_aux = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

mult_mat_div(:,2,3) = mult_aux(vec);

mult_mat_div(:,3,3) = mult_mat_div(:,2,3) - mult_mat_div(:,1,3);

%% Profits MORE CONCENTRATED

cd '../CODES/ROBUST_APPEND/PROFITMORE/CG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)


mult_aux = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

mult_mat_div(:,1,2) = mult_aux(vec);

cd '../CODES/ROBUST_APPEND/PROFITMORE/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)


mult_aux = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

mult_mat_div(:,2,2) = mult_aux(vec);

mult_mat_div(:,3,2) = mult_mat_div(:,2,2) - mult_mat_div(:,1,2);



%% Plot


colorCG = [1.00 0.20 0.30];
colorMG = [0.30 0.40 0.90];
% 
vec_plot = [2;3;6]; 

mat_graph_div = ones(6,2);
mat_graph_div(:,1) = squeeze(mult_mat_div(:,1,3));
mat_graph_div(:,2) = squeeze(mult_mat_div(:,2,3));

fig = figure(103); clf;

subplot(1,3,1)
b = bar([1 2 3],mat_graph_div(vec_plot,:));
b(1).FaceColor = colorCG;
b(2).FaceColor = colorMG;
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Less Concentrated','Interpreter','LaTex','Fontsize',27) % $\frac{\sum_{t}{Y_t - \bar{Y}}}{\sum_{t}{G_t - \bar{G}}}$
%xlabel('Year','Interpreter','LaTex','Fontsize',27)
set(gca, 'XTickLabel', {'1' '2','10'})
ylim([-0.02 0.4])
leg = legend('Constant Progressivity','Higher Progressivity');
set(leg,'Interpreter','LaTex','Fontsize',25,'Location','NorthEast')
legend boxoff
xlabel('Year','Interpreter','LaTex','Fontsize',27)





mat_graph_div = ones(6,2);
mat_graph_div(:,1) = squeeze(mult_mat_div(:,1,2));
mat_graph_div(:,2) = squeeze(mult_mat_div(:,2,2));




subplot(1,3,2)
b = bar([1 2 3],mat_graph_div(vec_plot,:));
b(1).FaceColor = colorCG;
b(2).FaceColor = colorMG;
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('More Concentrated','Interpreter','LaTex','Fontsize',27) % $\frac{\sum_{t}{Y_t - \bar{Y}}}{\sum_{t}{G_t - \bar{G}}}$
%xlabel('Year','Interpreter','LaTex','Fontsize',27)
set(gca, 'XTickLabel', {'1' '2','10'})
ylim([-0.02 0.4])
leg = legend('Constant Progressivity','Higher Progressivity');
set(leg,'Interpreter','LaTex','Fontsize',25,'Location','NorthEast')
legend boxoff
ylabel('Multipliers','Interpreter','LaTex','Fontsize',27)
xlabel('Year','Interpreter','LaTex','Fontsize',27)



mat_graph_diff= ones(6,3);
mat_graph_diff(:,1) = squeeze(mult_mat_div(:,3,1));
mat_graph_diff(:,2) = squeeze(mult_mat_div(:,3,3));
mat_graph_diff(:,3) = squeeze(mult_mat_div(:,3,2));

%fig = figure(102); clf;
color_rob2   = [0.95 0.75 0.1];
color_rob3   = [4 139 154]/255; %[0.2 0.6 0.4];
color_rob1 = 1.4*[0.6 0.6 0.6]; %[0.75 0.7 0.5];  %[161 149 121]/255; %[0.3 0.2 0.45];
color_rob4   = 1.2*[0.375 0.25 0.5625];

subplot(1,3,3)
b = bar([1 2 3],mat_graph_diff(vec_plot,:));
b(1).FaceColor = color_rob1;
b(2).FaceColor = color_rob2;
b(3).FaceColor = color_rob3;
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Difference in Multipliers','Interpreter','LaTex','Fontsize',27) % $\frac{\sum_{t}{Y_t - \bar{Y}}}{\sum_{t}{G_t - \bar{G}}}$
xlabel('Year','Interpreter','LaTex','Fontsize',27)
set(gca, 'XTickLabel', {'1' '2','10'})
ylim([-0.02 0.4])
leg = legend('Benchmark','Less Concentrated','More Concentrated');
set(leg,'Interpreter','LaTex','Fontsize',25,'Location','NorthWest')
legend boxoff





fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd ../FIGURES
    print('FIGURE_17','-dpng','-r0')
end




