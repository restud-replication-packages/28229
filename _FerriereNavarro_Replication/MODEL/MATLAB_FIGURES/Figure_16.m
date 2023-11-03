clear; close all; clc;

%% folders

save_graph = 1;
 
cdir  = pwd;

fdirv = {'../CODES/BENCHMARK/CG/OUTPUT','../CODES/BENCHMARK/MG/OUTPUT',...
        '../CODES/ROBUST_APPEND/HETB/CG/OUTPUT','../CODES/ROBUST_APPEND/HETB/MG/OUTPUT'};

Nfdir = numel(fdirv);

T_TR  = 100;

GAGG_mat = 189*ones(T_TR,Nfdir);
YAGG_mat = 189*ones(T_TR,Nfdir);
MULT_mat = 189*ones(T_TR,Nfdir);


for ic = 1:Nfdir
    fname = fdirv{ic};
    cd(fname)
    
    %---Steady-State
    load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, FC]
    
    lbd  = AGGREGATES(1); G = AGGREGATES(2); Cagg = AGGREGATES(3);
    Yagg = AGGREGATES(4); Lagg = AGGREGATES(5); EMP  = AGGREGATES(6);
    Kagg = AGGREGATES(7); DB   = AGGREGATES(8); Iagg = AGGREGATES(9);
    
    load prices.txt;
    wge = prices(1); rk = prices(2);  qk = prices(3); %     ra = prices(1); wgeH = prices(2); wge = prices(3); rk = prices(4); qk = prices(5);
    
    epsW = 7.0;
    wgeH = ((epsW-1.0D0)/epsW)*wge;
    
    %---Transition
    load G_TR.txt;      G_TR    = G_TR(1:T_TR);
    load Yagg_TR.txt;   Yagg_TR = Yagg_TR(1:T_TR);
    
    mult = cumsum(Yagg_TR-Yagg)./cumsum(G_TR-G);
    
    GAGG_mat(:,ic) = G_TR;
    YAGG_mat(:,ic) = 100*log(Yagg_TR/Yagg);    
    MULT_mat(:,ic) = mult;
    
    cd(cdir)
end

%% Plots

colorCG = [1.00 0.20 0.30];
colorMG = [0.30 0.40 0.90];

color_rob2   = [0.95 0.75 0.1];
color_rob3   = [4 139 154]/255; 
color_rob1 = 1.4*[0.6 0.6 0.6]; 
color_rob4   = 1.2*[0.375 0.25 0.5625];

vec = [4; 8; 40];  Nvec = numel(vec);

MATx      = 189*ones(Nvec,2);
MATx(:,1) = MULT_mat(vec,3);
MATx(:,2) = MULT_mat(vec,4);

MATd      = 189*ones(Nvec,2);
MATd(:,1) = MULT_mat(vec,2)-MULT_mat(vec,1);
MATd(:,2) = MULT_mat(vec,4)-MULT_mat(vec,3);

fig = figure(104); clf;

subplot(1,2,1)
b = bar([1 2 3],MATx);
b(1).FaceColor = colorCG;
b(2).FaceColor = colorMG;
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('$(\beta,B)$ Model','Interpreter','LaTex','Fontsize',27) % $\frac{\sum_{t}{Y_t - \bar{Y}}}{\sum_{t}{G_t - \bar{G}}}$
%xlabel('Year','Interpreter','LaTex','Fontsize',27)
set(gca, 'XTickLabel', {'1' '2','10'})
ylim([-0.01 0.55])
leg = legend('Constant Progressivity','Higher Progressivity');
set(leg,'Interpreter','LaTex','Fontsize',25,'Location','NorthWest')
legend boxoff
ylabel('Multipliers','Interpreter','LaTex','Fontsize',27)
xlabel('Year','Interpreter','LaTex','Fontsize',27)

subplot(1,2,2)
b = bar([1 2 3],MATd);
b(1).FaceColor = color_rob1;
b(2).FaceColor = color_rob3;
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Difference in Multipliers','Interpreter','LaTex','Fontsize',27) % $\frac{\sum_{t}{Y_t - \bar{Y}}}{\sum_{t}{G_t - \bar{G}}}$
xlabel('Year','Interpreter','LaTex','Fontsize',27)
set(gca, 'XTickLabel', {'1' '2','10'})
ylim([0 0.35])
leg = legend('Benchmark','$(\beta,B)$ Model');
set(leg,'Interpreter','LaTex','Fontsize',25,'Location','NorthWest')
legend boxoff


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd ../FIGURES
    print('FIGURE_16','-dpng','-r0')
end






