clear;  close all; clc;

cdir = pwd;

savef = 1;  % =1 to save figures

SSdir = '../CODES/BENCHMARK/CALIB/OUTPUT'; 

CGdir_1 = '../CODES/BENCHMARK/CG_DECOMP/OUTPUT_1'; % all sequences 
CGdir_2 = '../CODES/BENCHMARK/CG_DECOMP/OUTPUT_2'; % indirect
CGdir_3 = '../CODES/BENCHMARK/CG_DECOMP/OUTPUT_3'; % direct

MGdir_1 = '../CODES/BENCHMARK/MG_DECOMP/OUTPUT_1';
MGdir_2 = '../CODES/BENCHMARK/MG_DECOMP/OUTPUT_2';
MGdir_3 = '../CODES/BENCHMARK/MG_DECOMP/OUTPUT_3';


%% Load SS results

cd(SSdir)

%%%  Steady-State
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, Iagg]

lbd  = AGGREGATES(1); Gagg = AGGREGATES(2); Cagg = AGGREGATES(3);
Yagg = AGGREGATES(4); Lagg = AGGREGATES(5); EMP  = AGGREGATES(6);
Kagg = AGGREGATES(7); DB   = AGGREGATES(8); Iagg = AGGREGATES(9);

load prices.txt;
wge = prices(1); rk = prices(2); qk = prices(3);

epsW = 7; MCW  = (epsW-1)/epsW; wgeH = MCW*wge;

cd(cdir)

%% Load CG results

cd(CGdir_1)
load G_TR.txt;
load div_TR.txt; load lbd_TR.txt; load div_nobank_TR.txt;
load r_TR.txt; load wgeH_TR.txt;
load Yagg_TR.txt; load Avgetax.txt

div_TR = div_TR./Yagg_TR;
div_nobank_TR = div_nobank_TR./Yagg_TR;


CG.div   = div_TR;
CG.div_nobank   = div_nobank_TR;
CG.Tax   = Avgetax-Avgetax(end);
CG.wH    = wgeH_TR;
CG.r     = r_TR;

T_TR = numel(G_TR);

load Lagg_TR.txt;  CG.L_all = 100*(log(Lagg_TR)-log(Lagg));
load Cagg_TR.txt;  CG.C_all = 100*(log(Cagg_TR)-log(Cagg));

cd(cdir)

% Load CG indirect
cd(CGdir_2)
load Lagg_TR.txt;  CG.L_ind = 100*(log(Lagg_TR)-log(Lagg));
load Cagg_TR.txt;  CG.C_ind = 100*(log(Cagg_TR)-log(Cagg));
cd(cdir)

% Load CG taxes
cd(CGdir_3)
load Lagg_TR.txt;  CG.L_tax = 100*(log(Lagg_TR)-log(Lagg));
load Cagg_TR.txt;  CG.C_tax = 100*(log(Cagg_TR)-log(Cagg));
cd(cdir)
 

%%
tlb = 0; tub = 15;

index = [4;12;40];
indexy = [1 2 3];

colorW   = [0.95 0.75 0.1];
colorr   = [4 139 154]/255; 
colorind = [0.75 0.55 0];
colordiv = 1.4*[0.6 0.6 0.6];

colorT   = 1.2*[0.375 0.25 0.5625];

%
figure(101)
subplot(2,3,1)
plot(0:T_TR-1,100*(CG.Tax),'LineWidth',3.2,'Color',colorT,'LineStyle','--'); hold on
plot(0:T_TR-1,100*(CG.r-rk),'LineWidth',3.25,'Color',colordiv,'Marker','square','MarkerSize',7,'MarkerFaceColor',colorr,'MarkerEdgeColor',colorr);
plot(0:T_TR-1,100*CG.div_nobank,'LineWidth',3.25,'Color',colordiv,'LineStyle',':','Marker','o','MarkerSize',7,'MarkerFaceColor',colorind,'MarkerEdgeColor',colorind);
plot(0:T_TR-1,100*log(CG.wH/wgeH),'LineWidth',3.25,'Color',colordiv,'LineStyle','-.','Marker','diamond','MarkerSize',7,'MarkerFaceColor',colorW,'MarkerEdgeColor',colorW); hold on
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); hold off
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Taxes and Prices','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
leg = legend('$\Delta \tau_t^\ell$','$\Delta r_t$','$\Delta d_t/Y_t$','$\Delta \ln w_t^h$');
set(leg,'Interpreter','LaTex','Fontsize',25);
legend boxoff
ylabel('percentage point','Interpreter','LaTex','Fontsize',27)
ylim([-0.07 0.12])
yticks([-0.15 -0.10 -0.05 0 0.05 0.10 0.15])% -2*pi -pi 0 pi 2*pi 3*pi])
%xlabel('Quarter','Interpreter','LaTex','Fontsize',27)

hold off

subplot(2,3,2)
plot(0:T_TR-1,CG.L_tax,'LineWidth',3.25,'Color',colorT,'LineStyle','--');hold on
plot(0:T_TR-1,CG.L_ind,'LineWidth',3.25,'Color',colordiv,'LineStyle',':','MarkerSize',7,'Marker','+');
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
leg = legend('Direct effects of taxes','Indirect effects');
set(leg,'Interpreter','LaTex','Fontsize',25);
legend boxoff
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Labor','Interpreter','LaTex','Fontsize',27)
ylabel('$\%$ from steady-state','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
%xlabel('Quarter','Interpreter','LaTex','Fontsize',27)

ylim([-0.05 0.1])

text(1,-0.075,0.1,'Constant Progressivity','Interpreter','LaTex','Fontsize',30)

subplot(2,3,3)
plot(0:T_TR-1,CG.C_ind,'LineWidth',3.25,'Color',colordiv,'LineStyle',':','MarkerSize',7,'Marker','+');hold on
plot(0:T_TR-1,CG.C_tax,'LineWidth',3.25,'Color',colorT,'LineStyle','--','MarkerSize',7);
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Consumption','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
%xlabel('Quarter','Interpreter','LaTex','Fontsize',27)


%% Load MG results

cd(MGdir_1)
load G_TR.txt;
load div_TR.txt; load lbd_TR.txt; load div_nobank_TR.txt;
load r_TR.txt; load wgeH_TR.txt;
load Yagg_TR.txt; load Avgetax.txt

div_TR = div_TR./Yagg_TR;
div_nobank_TR = div_nobank_TR./Yagg_TR;

MG.div   = div_TR;
MG.div_nobank   = div_nobank_TR;
MG.Tax   = Avgetax-Avgetax(end);
MG.wH    = wgeH_TR;
MG.r     = r_TR;

T_TR = numel(G_TR);

load Lagg_TR.txt;  MG.L_all = 100*(log(Lagg_TR)-log(Lagg));
load Cagg_TR.txt;  MG.C_all = 100*(log(Cagg_TR)-log(Cagg));

cd(cdir)

% Load MG indirect
cd(MGdir_2)
load Lagg_TR.txt;  MG.L_ind = 100*(log(Lagg_TR)-log(Lagg));
load Cagg_TR.txt;  MG.C_ind = 100*(log(Cagg_TR)-log(Cagg));
cd(cdir)

% Load MG taxes
cd(MGdir_3)
load Lagg_TR.txt;  MG.L_tax = 100*(log(Lagg_TR)-log(Lagg));
load Cagg_TR.txt;  MG.C_tax = 100*(log(Cagg_TR)-log(Cagg));
cd(cdir)


%%

%figure(12)
subplot(2,3,4)
plot(0:T_TR-1,100*(MG.Tax),'LineWidth',3.2,'Color',colorT,'LineStyle','--'); hold on
plot(0:T_TR-1,100*log(MG.wH/wgeH),'LineWidth',3.25,'Color',colordiv,'LineStyle','-.','Marker','diamond','MarkerSize',9,'MarkerFaceColor',colorW,'MarkerEdgeColor',colorW); hold on
plot(0:T_TR-1,100*(MG.r-rk),'LineWidth',3.25,'Color',colordiv,'Marker','square','MarkerSize',9,'MarkerFaceColor',colorr,'MarkerEdgeColor',colorr);
plot(0:T_TR-1,100*MG.div_nobank,'LineWidth',3.25,'Color',colordiv,'LineStyle',':','Marker','o','MarkerSize',9,'MarkerFaceColor',colorind,'MarkerEdgeColor',colorind);
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); hold off
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Taxes and Prices','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
% leg = legend('$\Delta \tau_t^\ell$','$\Delta \ln w_t$ $(\%)$','$\Delta r_t$','$\Delta d_t/Y_t$');
% set(leg,'Interpreter','LaTex','Fontsize',25);
% legend boxoff
ylabel('percentage point','Interpreter','LaTex','Fontsize',27)
ylim([-0.07 0.12])
yticks([-0.15 -0.10 -0.05 0 0.05 0.10 0.15])% -2*pi -pi 0 pi 2*pi 3*pi])
hold off

xlabel('Quarter','Interpreter','LaTex','Fontsize',27)

subplot(2,3,5)
plot(0:T_TR-1,MG.L_tax,'LineWidth',3.25,'Color',colorT,'LineStyle','--');hold on
plot(0:T_TR-1,MG.L_ind,'LineWidth',3.25,'Color',colordiv,'LineStyle',':','MarkerSize',7,'Marker','+');
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
% leg = legend('Direct effects of taxes','Indirect effects');
% set(leg,'Interpreter','LaTex','Fontsize',25);
% legend boxoff
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Labor','Interpreter','LaTex','Fontsize',27)
ylabel('$\%$ from steady-state','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])

ylim([-0.05 0.1])
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)
text(-0,-0.15,0,'Constant Progressivity','Interpreter','LaTex','Fontsize',27)



subplot(2,3,6)
plot(0:T_TR-1,MG.C_tax,'LineWidth',3.25,'Color',colorT,'LineStyle','--');hold on
plot(0:T_TR-1,MG.C_ind,'LineWidth',3.25,'Color',colordiv,'LineStyle',':','MarkerSize',7,'Marker','+');
plot(0:T_TR-1,zeros(T_TR,1),'color',[0.6 0.6 0.6],'LineWidth',2.25,'LineStyle','--'); 
set(gca,'XGrid','on','YGrid','on','Fontsize',27,'TickLabelInterpreter', 'LaTex')  
title('Consumption','Interpreter','LaTex','Fontsize',27)
xlim([tlb tub])
xlabel('Quarter','Interpreter','LaTex','Fontsize',27)


fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (savef == 1 )
    cd ../FIGURES
    print('FIGURE_4','-dpng','-r0')
end

