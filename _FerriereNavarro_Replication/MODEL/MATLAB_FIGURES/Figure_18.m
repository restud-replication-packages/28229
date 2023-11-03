clear;  close all; clc;

%% load results
save_graph = 1; % 1 to save graphs 

cdir  = pwd;  % current directory

vec = [1; 4; 8; 40];

TT = 40;

% Delta multiplier
dmult_mat_wages = 0*ones(TT,3,3); % 4 periods, Benchmark/Flat LPE/Flat MPC, Benchmark/Flexible wages/Sticky wages


%% Benchmark 

cd '../CODES/BENCHMARK/CG/OUTPUT'
load AGGREGATES.txt; G    = AGGREGATES(2); Yagg = AGGREGATES(4);
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_CG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

cd '../CODES/BENCHMARK/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_MG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
dmult_mat_wages(:,1,1) = mult_mat_MG(1:TT) - mult_mat_CG(1:TT);

%% Benchmark FLAT LPE

cd '../CODES/ROBUST/FLATLPE/CG/OUTPUT'
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]
G    = AGGREGATES(2); Yagg = AGGREGATES(4);
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_CG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

cd '../CODES/ROBUST/FLATLPE/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_MG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
dmult_mat_wages(:,2,1) = mult_mat_MG(1:TT) - mult_mat_CG(1:TT);

%% Benchmark FLAT MPC

cd '../CODES/ROBUST/FLATMPC/CG/OUTPUT'
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]
G    = AGGREGATES(2); Yagg = AGGREGATES(4);
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_CG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

cd '../CODES/ROBUST/FLATMPC/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_MG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
dmult_mat_wages(:,3,1) = mult_mat_MG(1:TT) - mult_mat_CG(1:TT);


%% FLEXIBLE WAGES 

cd '../CODES/ROBUST/FLEXWAGES/CG/OUTPUT'
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]
G    = AGGREGATES(2); Yagg = AGGREGATES(4);
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_CG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

cd '../CODES/ROBUST/FLEXWAGES/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_MG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
dmult_mat_wages(:,1,2) = mult_mat_MG(1:TT) - mult_mat_CG(1:TT);


%% FLEX WAGES FLAT LPE


cd '../CODES/ROBUST_APPEND/FLATLPE_FLEXWAGES/CG/OUTPUT'
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]
G    = AGGREGATES(2); Yagg = AGGREGATES(4);
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_CG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

cd '../CODES/ROBUST_APPEND/FLATLPE_FLEXWAGES/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_MG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
dmult_mat_wages(:,2,2) = mult_mat_MG(1:TT) - mult_mat_CG(1:TT);

%% Benchmark FLAT MPC

cd '../CODES/ROBUST_APPEND/FLATMPC_FLEXWAGES/CG/OUTPUT'
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]
G    = AGGREGATES(2); Yagg = AGGREGATES(4);
%%% Transition
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_CG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

cd '../CODES/ROBUST_APPEND/FLATMPC_FLEXWAGES/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_MG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
dmult_mat_wages(:,3,2) = mult_mat_MG(1:TT) - mult_mat_CG(1:TT);


%% STICKY WAGES 

%fname = 'TRCG_SW6';
cd '../CODES/ROBUST_APPEND/STICKWAGES/CG/OUTPUT'
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]
G    = AGGREGATES(2); Yagg = AGGREGATES(4);
%%% Transition
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_CG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

cd '../CODES/ROBUST_APPEND/STICKWAGES/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_MG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
dmult_mat_wages(:,1,3) = mult_mat_MG(1:TT) - mult_mat_CG(1:TT);

%% Benchmark FLAT LPE

%fname = 'TRCG_SW6_FlatLPE';
cd '../CODES/ROBUST_APPEND/FLATLPE_STICKWAGES/CG/OUTPUT'
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]
G    = AGGREGATES(2); Yagg = AGGREGATES(4);
%%% Transition
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_CG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

cd '../CODES/ROBUST_APPEND/FLATLPE_STICKWAGES/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_MG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
dmult_mat_wages(:,2,3) = mult_mat_MG(1:TT) - mult_mat_CG(1:TT);

%% Benchmark FLAT MPC

cd '../CODES/ROBUST_APPEND/FLATMPC_STICKWAGES/CG/OUTPUT'
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]
G    = AGGREGATES(2); Yagg = AGGREGATES(4);
%%% Transition
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_CG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);

cd '../CODES/ROBUST_APPEND/FLATMPC_STICKWAGES/MG/OUTPUT'
load Yagg_TR.txt; load G_TR.txt;
cd(cdir)
mult_mat_MG = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
dmult_mat_wages(:,3,3) = mult_mat_MG(1:TT) - mult_mat_CG(1:TT);



%% Plot


colorCG = [1.00 0.20 0.30];
colorMG = [0.30 0.40 0.90];
% 
vec_plot = [1;2;3;4]; 
%vec_plot = [1;3];%;6]; 



color_rob2   = [0.95 0.75 0.1];
color_rob3   = [4 139 154]/255; %[0.2 0.6 0.4];
color_rob1 = 1.4*[0.6 0.6 0.6]; %[0.75 0.7 0.5];  %[161 149 121]/255; %[0.3 0.2 0.45];
color_rob4   = 1.2*[0.375 0.25 0.5625];




colorB   = [0.75 0.20 0.32];
colorFhG = [0.95 0.75 0.10];
colorb98 = [4 139 154]/255;
colorb96 = [0.375 0.25 0.5625];


cbench = colorB;   % [0.4 0.4 0.9];
cflpe  = colorFhG; % [0.4 0.4 0.1];
cfmpc  = colorb98; % [0.9 0.4 0.1];



tlb = 0; tub = 15;

ylb = 0.08; yub = 0.3;

figure(1301)


subplot(1,3,1)
plot(0:TT-1,dmult_mat_wages(:,1,1),'color',cbench,'LineWidth',3.25,'Marker','|','MarkerSize',7); hold on
plot(0:TT-1,dmult_mat_wages(:,2,1) ,'color',cflpe ,'LineWidth',3.25,'Marker','O','MarkerSize',7);
plot(0:TT-1,dmult_mat_wages(:,3,1) ,'color',cfmpc ,'LineWidth',3.25,'Marker','X','MarkerSize',7);
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Baseline: $\Theta_w=200$','Interpreter','LaTex','Fontsize',27)
ylabel('difference in output multipliers (p.p.)','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
leg = legend('Benchmark','Flatter \textit{lpe}','Lower \textit{mpc}');
set(leg,'Interpreter','LaTex','Fontsize',23)
legend boxoff
xlim([tlb tub])
ylim([ylb yub])
hold off



subplot(1,3,2)
plot(0:TT-1,dmult_mat_wages(:,1,2),'color',cbench,'LineWidth',3.25,'Marker','|','MarkerSize',7); hold on
plot(0:TT-1,dmult_mat_wages(:,2,2) ,'color',cflpe ,'LineWidth',3.25,'Marker','O','MarkerSize',7);
plot(0:TT-1,dmult_mat_wages(:,3,2) ,'color',cfmpc ,'LineWidth',3.25,'Marker','X','MarkerSize',7);
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Flexible Wages: $\Theta_w=0$','Interpreter','LaTex','Fontsize',27)
%ylabel('difference in output multipliers (p.p.)','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
% leg = legend('Benchmark','flat \textit{lpe}','lower \textit{mpc}');
% set(leg,'Interpreter','LaTex','Fontsize',23)
% legend boxoff
xlim([tlb tub])
ylim([ylb yub])

hold off



subplot(1,3,3)
plot(0:TT-1,dmult_mat_wages(:,1,3),'color',cbench,'LineWidth',3.25,'Marker','|','MarkerSize',7); hold on
plot(0:TT-1,dmult_mat_wages(:,2,3) ,'color',cflpe ,'LineWidth',3.25,'Marker','O','MarkerSize',7);
plot(0:TT-1,dmult_mat_wages(:,3,3) ,'color',cfmpc ,'LineWidth',3.25,'Marker','X','MarkerSize',7);
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('More Rigid Wages: $\Theta_w=600$','Interpreter','LaTex','Fontsize',27)
%ylabel('difference in output multipliers (p.p.)','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
% leg = legend('Benchmark','flat \textit{lpe}','lower \textit{mpc}');
% set(leg,'Interpreter','LaTex','Fontsize',23)
% legend boxoff
xlim([tlb tub])
ylim([ylb yub])

hold off


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd ../FIGURES
    print('FIGURE_18','-dpng','-r0')
end


