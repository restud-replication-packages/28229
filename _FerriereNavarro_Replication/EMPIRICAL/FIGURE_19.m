clear; close all; clc

%% directories

cdir = pwd;
odir = ('output');
fdir = ('figures');

%% load data

save_graph = 1;

cd(odir)
shocks_DATA  = readtable('shocks.csv','Delimiter',',');
quarter = table2array(shocks_DATA(:,1));
news    = table2array(shocks_DATA(:,2)); % news = 100*news;
ind     = table2array(shocks_DATA(:,3)); ind(1) = ind(2);
bp      = table2array(shocks_DATA(:,4));
cd(cdir)

%% computations - bar

colorPci  = [0.40 0.60 1.00];
colorNci  = [0.80 0.30 0.20];
%---bp progressive
bpx = bp(9:end); indx = ind(9:end); data = bpx(indx==1);
Nbp = 20; bpmin = min(data); bpmax = max(data); bp_xprg = linspace(bpmin,bpmax,Nbp)';
bp_bprg = histbins_fixvec(data,bp_xprg);

%---bp non-progressive
data = bpx(indx==0); 
Nbp = 20; bpmin = min(data); bpmax = max(data); bp_xnprg = linspace(bpmin,bpmax,Nbp)';
bp_bnprg = histbins_fixvec(data,bp_xnprg);

%---news progressive
data = news(indx==1);
Nnews = 20; newsmin = min(data); newsmax = max(data); news_xprg = linspace(newsmin,newsmax,Nnews)';
news_bprg = histbins_fixvec(data,news_xprg);

%---news non-progressive
data = news(indx==0);
Nnews = 20; newsmin = min(data); newsmax = max(data); news_xnprg = linspace(newsmin,newsmax,Nnews)';
news_bnprg = histbins_fixvec(data,news_xnprg);

fig = figure(1013);  clf;
subplot(1,2,1)
bar(bp_xprg ,bp_bprg ,'FaceColor',colorPci,'EdgeColor',colorPci); hold on
bar(bp_xnprg,bp_bnprg,'FaceColor',colorNci,'EdgeColor',colorNci,'FaceAlpha',.75,'EdgeAlpha',.5)
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('BP distribution','Interpreter','LaTex','Fontsize',30)
ylabel('density ($\%$)','Interpreter','LaTex','Fontsize',27)
xlabel('shock','Interpreter','LaTex','Fontsize',27)
ylim([0.0 0.50])
leg = legend('Progressive','Non-Progressive');
set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northeast')
legend boxoff
hold off
subplot(1,2,2)
bar(news_xprg ,news_bprg ,'FaceColor',colorPci,'EdgeColor',colorPci); hold on
bar(news_xnprg,news_bnprg,'FaceColor',colorNci,'EdgeColor',colorNci,'FaceAlpha',.75,'EdgeAlpha',.5)
set(gca,'XGrid','off','YGrid','on','Fontsize',21) %
title('RZ distribution','Interpreter','LaTex','Fontsize',30)
ylabel('density ($\%$)','Interpreter','LaTex','Fontsize',27)
xlabel('shock','Interpreter','LaTex','Fontsize',27)
ylim([0.0 1.0])
leg = legend('Progressive','Non-Progressive');
set(leg,'Interpreter','LaTex','Fontsize',30,'Location','northeast')
legend boxoff
hold off

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    cd(fdir)
    print('FIGURE_19','-dpng','-r0')    
    cd(cdir)
end

