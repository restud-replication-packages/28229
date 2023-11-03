clear; close all; clc;

%% set up

cdir = pwd;

savef = 1;  % =1 to save figures

SSdir = '../CODES/BENCHMARK/CALIB/OUTPUT'; 
CGdir = '../CODES/BENCHMARK/CG_DECOMP/OUTPUT_1'; 
MGdir = '../CODES/BENCHMARK/MG_DECOMP/OUTPUT_1';

qind = [0.5;0.8;1]; Nqind = numel(qind);                                    % Size of groups


%% Load Steady State Results

cd(SSdir) %%%  Steady-State

load AGGREGATES.txt;

lbd  = AGGREGATES(1); Gagg = AGGREGATES(2); Cagg = AGGREGATES(3);
Yagg = AGGREGATES(4); Lagg = AGGREGATES(5); EMP  = AGGREGATES(6);
Kagg = AGGREGATES(7); DB   = AGGREGATES(8); Iagg = AGGREGATES(9);

load prices.txt;
wge = prices(1); rk = prices(2); qk = prices(3);

epsW = 7; MCW  = (epsW-1)/epsW; wgeH = MCW*wge;

load tax.txt;
gma = tax(2); TF = tax(3);
hbar = 0.33;

load avec.txt; load xvec.txt; load Svec.txt; load params.txt;

Na = params(1);  Nx = params(2);  Nbeta = params(3);  NS = Na*Nx;  Nh = 2;

%---steady-state policies and measure
load hpol.txt; load mu.txt;

cd(cdir)

% reshape steady-state

hpolx = hpol;
hpol = 189*ones(NS,Nbeta,Nh);
inn = 1;
for ib = 1:Nbeta
    ou = inn + NS-1;
    hpol(:,ib,1) = hpolx(inn:ou,1);
    hpol(:,ib,2) = hpolx(inn:ou,2);
    inn = ou+1;
end

% tax x group

yl_lng   = 189*ones(NS*Nbeta,1);
hp_lng   = 189*ones(NS*Nbeta,1);
mu_lng   = 189*ones(NS*Nbeta,1);
muw_lng  = 189*ones(NS*Nbeta,1);
taul_lng = 189*ones(NS*Nbeta,1);

ih=2;
inn = 1;
for ib = 1:Nbeta
    
    ylx   = wgeH*(Svec(:,2)*hbar); ylxx = max(ylx,1e-13);
    
    taulx = (ylx - lbd*(ylx.^(1-gma)))./ylxx;
    
    ou = inn + NS-1;
    yl_lng(inn:ou)    = ylx;    
    hp_lng(inn:ou)    = hpol(:,ib,ih);
    mu_lng(inn:ou)    = mu(:,ib);
    muw_lng(inn:ou)   = mu_lng(inn:ou).*hp_lng(inn:ou);
    taul_lng(inn:ou)  = taulx;
    
    inn = ou+1;
end
mu_lng  = mu_lng /sum(mu_lng(:));
muw_lng = muw_lng/sum(muw_lng(:));



%---sort by yl
[yls,iyls] = sort(yl_lng);
taul_xyl   = taul_lng(iyls);
hp_xyl     = hp_lng(iyls);
mu_xyl     = mu_lng(iyls);
muw_xyl    = muw_lng(iyls);
CDF_xyl    = cumsum(muw_xyl);

qtaul_xy = 189*ones(Nqind,1);  

for iq = 1:Nqind
    if iq == 1
        inn = (CDF_xyl<=qind(iq));
    else
        inn = (CDF_xyl>qind(iq-1)).*(CDF_xyl<=qind(iq));        
    end
    mm  = muw_xyl.*inn;
    
    xx  = taul_xyl.*inn; qtaul_xy(iq) = sum(xx.*mm)/sum(mm);
end

qtaul_xy_SS = qtaul_xy;


%% CG case 

xdir = CGdir; CGcase = 1;

run FIGURE_TAX_COMPUTATIONS

qtaul_xy_CG = qtaul_xy_t;
qyl_xy_CG   = qyl_xy_t;
mp_t_CG     = mp_t;

CG_IRF_G = 100*log(G_TR/Gagg);
CG_DB    = 100*(DB_TR/Yagg);

%% MG case 

xdir = MGdir; CGcase = 0;

run FIGURE_TAX_COMPUTATIONS

qtaul_xy_MG = qtaul_xy_t;
qyl_xy_MG   = qyl_xy_t;
mp_t_MG     = mp_t;

MG_IRF_G = 100*log(G_TR/Gagg);
MG_DB    = 100*(DB_TR/Yagg);

%% Make plot

colorCG = [1.00 0.20 0.30];
colorMG = [0.30 0.40 0.90];

tlb = 0; tub = 15;

fig = figure(10001); clf;
subplot(2,2,1)
plot(0:T_TR-1,CG_IRF_G,'LineWidth',3.25,'Color',colorCG,'Marker','o','MarkerSize',4.25); hold on
plot(0:T_TR-1,MG_IRF_G,'LineWidth',3.25,'Color',colorMG); hold on
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Government Spending','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
ylim([-0.05 1])
ylabel('$\%$ from steady-state','Interpreter','LaTex','Fontsize',27)
hold off

subplot(2,2,2)
plot(0:T_TR-1,CG_DB(1:T_TR),'LineWidth',3.25,'Color',colorCG,'Marker','o','MarkerSize',4.25); hold on
plot(0:T_TR-1,MG_DB(1:T_TR),'LineWidth',3.25,'Color',colorMG); hold on
plot(0:T_TR-1,zeros(T_TR,1)+CG_DB(end),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Public Debt','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
ylim([99.99 100.151])
ylabel('$\%$ of steady-state GDP','Interpreter','LaTex','Fontsize',27)
hold off
leg = legend('Constant Progressivity','Higher Progressivity','Steady-State','Location','SouthWest');
set(leg,'Interpreter','LaTex','Fontsize',25,'Location','NorthEast')
legend boxoff


%---tax botom-50%
subplot(2,3,4); iq = 1;
plot(0:T_TR-1,100*qtaul_xy_CG(iq,1:T_TR),'LineWidth',3.25,'Color',colorCG,'Marker','o','MarkerSize',4.25); hold on
plot(0:T_TR-1,100*qtaul_xy_MG(iq,1:T_TR),'LineWidth',3.25,'Color',colorMG); 
plot(0:T_TR-1,zeros(T_TR,1)+100*qtaul_xy_SS(iq),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Tax Rate Bottom $50\%$','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
ylabel('percentage points','Interpreter','LaTex','Fontsize',27)
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
ylim([100*qtaul_xy_SS(iq)-0.02 100*qtaul_xy_SS(iq)+0.25])
yticks([23.6 23.7 23.8])
hold off
subplot(2,3,5); iq = 2;
plot(0:T_TR-1,100*qtaul_xy_CG(iq,1:T_TR),'LineWidth',3.25,'Color',colorCG,'Marker','o','MarkerSize',4.25); hold on
plot(0:T_TR-1,100*qtaul_xy_MG(iq,1:T_TR),'LineWidth',3.25,'Color',colorMG); 
plot(0:T_TR-1,zeros(T_TR,1)+100*qtaul_xy_SS(iq),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Tax Rate Bottom $50\%$ to $80\%$','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
yticks([30.65 30.75 30.85])
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
ylim([100*qtaul_xy_SS(iq)-0.02 100*qtaul_xy_SS(iq)+0.25])
hold off
subplot(2,3,6); iq = 3;
plot(0:T_TR-1,100*qtaul_xy_CG(iq,1:T_TR),'LineWidth',3.25,'Color',colorCG,'Marker','o','MarkerSize',4.25); hold on
plot(0:T_TR-1,100*qtaul_xy_MG(iq,1:T_TR),'LineWidth',3.25,'Color',colorMG); 
plot(0:T_TR-1,zeros(T_TR,1)+100*qtaul_xy_SS(iq),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Tax Rate Top $20\%$','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
yticks([36.25 36.35 36.45])
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
ylim([100*qtaul_xy_SS(iq)-0.02 100*qtaul_xy_SS(iq)+0.25])
hold off



fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (savef == 1 )
    cd ../FIGURES
    print('FIGURE_1','-dpng','-r0')
end

