clear; close all; clc;

cdir = pwd;

fdirv = {'../CODES/BENCHMARK/CG','../CODES/BENCHMARK/MG',...
        '../CODES/ROBUST/FLATLPE/CG','../CODES/ROBUST/FLATLPE/MG',...
        '../CODES/ROBUST/FLATMPC/CG','../CODES/ROBUST/FLATMPC/MG'};

Nfdir = numel(fdirv);

save_graph = 1; % =1 to save graph


%% Load transitions

T_TR = 250;

GAGG_mat = 189*ones(T_TR,Nfdir); wge_mat  = 189*ones(T_TR,Nfdir); lbd_mat  = 189*ones(T_TR,Nfdir);
Pi_mat   = 189*ones(T_TR,Nfdir); rk_mat   = 189*ones(T_TR,Nfdir); YAGG_mat = 189*ones(T_TR,Nfdir);
CAGG_mat = 189*ones(T_TR,Nfdir); LAGG_mat = 189*ones(T_TR,Nfdir); IAGG_mat = 189*ones(T_TR,Nfdir);
MULT_mat = 189*ones(T_TR,Nfdir);

for ic = 1:Nfdir
    fname = fdirv{ic};
    cd(strcat(fname,'/OUTPUT'))
    load Lagg_TR.txt; load Yagg_TR.txt; load Cagg_TR.txt; load G_TR.txt;
    load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]

    G = AGGREGATES(2); Cagg = AGGREGATES(3);
    Yagg = AGGREGATES(4); Lagg = AGGREGATES(5); 

    
    mult = cumsum(Yagg_TR(1:T_TR)-Yagg)./cumsum(G_TR(1:T_TR)-G_TR(end));
    
    GAGG_mat(:,ic) = G_TR(1:T_TR);
    YAGG_mat(:,ic) = 100*log(Yagg_TR(1:T_TR)/Yagg);
    CAGG_mat(:,ic) = 100*log(Cagg_TR(1:T_TR)/Cagg);
    LAGG_mat(:,ic) = 100*log(Lagg_TR(1:T_TR)/Lagg);    

    MULT_mat(:,ic) = mult;
    
    cd(cdir)
end    

%% computations

dL_bench = LAGG_mat(:,2)-LAGG_mat(:,1);
dL_flpe  = LAGG_mat(:,4)-LAGG_mat(:,3);
dL_fmpc  = LAGG_mat(:,6)-LAGG_mat(:,5);

dC_bench = CAGG_mat(:,2)-CAGG_mat(:,1);
dC_flpe  = CAGG_mat(:,4)-CAGG_mat(:,3);
dC_fmpc  = CAGG_mat(:,6)-CAGG_mat(:,5);

dM_bench = MULT_mat(:,2)-MULT_mat(:,1);
dM_flpe  = MULT_mat(:,4)-MULT_mat(:,3);
dM_fmpc  = MULT_mat(:,6)-MULT_mat(:,5);

dI_bench = IAGG_mat(:,2)-IAGG_mat(:,1);
dI_flpe  = IAGG_mat(:,4)-IAGG_mat(:,3);
dI_fmpc  = IAGG_mat(:,6)-IAGG_mat(:,5);


%% plot

colorB   = [0.75 0.20 0.32];
colorFhG = [0.95 0.75 0.10];
colorb98 = [4 139 154]/255;
colorb96 = [0.375 0.25 0.5625];


cbench = colorB;   % [0.4 0.4 0.9];
cflpe  = colorFhG; % [0.4 0.4 0.1];
cfmpc  = colorb98; % [0.9 0.4 0.1];
cfall  = [0.1 0.8 0.4];

tlb = 0; tub = 15;
figure(1301)
subplot(1,3,1)
plot(0:T_TR-1,dM_bench,'color',cbench,'LineWidth',3.25,'Marker','|','MarkerSize',7); hold on
plot(0:T_TR-1,dM_flpe ,'color',cflpe ,'LineWidth',3.25,'Marker','O','MarkerSize',7);
plot(0:T_TR-1,dM_fmpc ,'color',cfmpc ,'LineWidth',3.25,'Marker','X','MarkerSize',7);
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Output Multiplier','Interpreter','LaTex','Fontsize',27)
ylabel('difference in multipliers (p.p.)','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
% leg = legend('Benchmark','flat \textit{lpe}','lower \textit{mpc}');
% set(leg,'Interpreter','LaTex','Fontsize',23)
% legend boxoff
xlim([tlb tub])
hold off
subplot(1,3,2)
plot(0:T_TR-1,dL_bench,'color',cbench,'LineWidth',3.25,'Marker','|','MarkerSize',7); hold on
plot(0:T_TR-1,dL_flpe ,'color',cflpe ,'LineWidth',3.25,'Marker','O','MarkerSize',7);
plot(0:T_TR-1,dL_fmpc ,'color',cfmpc ,'LineWidth',3.25,'Marker','X','MarkerSize',7);
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Labor','Interpreter','LaTex','Fontsize',27)
ylabel('difference in response (\%)','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
leg = legend('Benchmark','Flatter \textit{lpe}','Lower \textit{mpc}');
set(leg,'Interpreter','LaTex','Fontsize',25)
legend boxoff
xlim([tlb tub])
hold off
subplot(1,3,3)
plot(0:T_TR-1,dC_bench,'color',cbench,'LineWidth',3.25,'Marker','|','MarkerSize',7); hold on
plot(0:T_TR-1,dC_flpe ,'color',cflpe ,'LineWidth',3.25,'Marker','O','MarkerSize',7);
plot(0:T_TR-1,dC_fmpc ,'color',cfmpc ,'LineWidth',3.25,'Marker','X','MarkerSize',7);
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Consumption','Interpreter','LaTex','Fontsize',27)
ylabel('difference in response (\%)','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd ../FIGURES
    print('FIGURE_5','-dpng','-r0')
end
