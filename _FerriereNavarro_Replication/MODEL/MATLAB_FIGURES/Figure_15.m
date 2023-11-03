clear; close all; clc;

%% load results

save_graph = 1;
cdir = pwd;

odirv = {'../CODES/BENCHMARK/CALIB/OUTPUT','../CODES/ROBUST_APPEND/HETB/CALIB/OUTPUT'};
     
NC = numel(odirv);

qvec = [0.2 0.4 0.6 0.8 1.0];  Nq = numel(qvec);

qas_xas_xC  = 189*ones(Nq,NC);  qyl_xyl_xC  = 189*ones(Nq,NC);
qmpc_xas_xC = 189*ones(Nq,NC);  qmpc_xyl_xC = 189*ones(Nq,NC);
qlpe_xas_xC = 189*ones(Nq,NC);  qlpe_xyl_xC = 189*ones(Nq,NC);
qwrk_xas_xC = 189*ones(Nq,NC);  qwrk_xyl_xC = 189*ones(Nq,NC);
qash_xas_xC = 189*ones(Nq,NC);

%% loop in calibs

for iC = 1:NC
    odir = odirv{iC};
    
    cd(odir)
    
    load avec.txt
    load xvec.txt
    load Svec.txt
    load mu.txt
    load mpc.txt
    load lpe.txt
    load yls.txt
    load hpol.txt;
    load params.txt
    cd(cdir)
    
    Na = params(1);  Nx = params(2);  Nbeta = params(3);  NS = Na*Nx;
    
    
    %---------------------------------------------------------------------------------------------
    %---reshape
    mu_rs   = 189*ones(Na,Nx,Nbeta);
    mpc_rs  = 189*ones(Na,Nx,Nbeta);
    lpe_rs  = 189*ones(Na,Nx,Nbeta);
    is = 1; 
    for ix = 1:Nx
    for ia = 1:Na
        mu_rs(ia,ix,:)  = mu(is,:);
        mpc_rs(ia,ix,:) = mpc(is,:);
        lpe_rs(ia,ix,:) = lpe(is,:);   
        is = is+1;
    end
    end

    mua = 189*ones(Na,1);  mua_xb = 189*ones(Na,Nbeta);
    for ia = 1:Na
        mm      = mu_rs(ia,:,:);
        mua(ia) = sum(mm(:));
        for ib = 1:Nbeta
            mm            = mu_rs(ia,:,ib);
            mua_xb(ia,ib) = sum(mm(:));
        end
    end

    hwrk_rs = 189*ones(Na,Nx,Nbeta);
    ih = 2;
    inn = 1;
    for ib = 1:Nbeta    
        for ix = 1:Nx
        for ia = 1:Na
            hwrk_rs(ia,ix,ib) = hpol(inn,ih);

            inn = inn+1;
        end
        end
    end

    mu_lng   = 189*ones(NS*Nbeta,1);
    mpc_lng  = 189*ones(NS*Nbeta,1);
    lpe_lng  = 189*ones(NS*Nbeta,1);
    yl_lng   = 189*ones(NS*Nbeta,1);
    a_lng    = 189*ones(NS*Nbeta,1);
    x_lng    = 189*ones(NS*Nbeta,1);
    hwrk_lng = 189*ones(NS*Nbeta,1);
    inn = 1;
    for ib = 1:Nbeta
        ou = inn + NS-1;
        mu_lng(inn:ou)   = mu(:,ib);
        mpc_lng(inn:ou)  = mpc(:,ib);
        lpe_lng(inn:ou)  = lpe(:,ib);
        yl_lng(inn:ou)   = yls(:,ib);
        a_lng(inn:ou)    = Svec(:,1);
        x_lng(inn:ou)    = Svec(:,2);    
        inn = ou +1;

    end
    inn = 1;
    for ib = 1:Nbeta
        for ix = 1:Nx
        for ia = 1:Na
            hwrk_lng(inn) = hwrk_rs(ia,ix,ib);
            inn = inn + 1;
        end
        end
    end
    %---------------------------------------------------------------------------------------------
    %---------------------------------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------------------------------
    %---computation    

    %---sort by assets
    [ass, ias] = sort(a_lng);
    mu_xas   = mu_lng(ias);
    mpc_xas  = mpc_lng(ias);
    lpe_xas  = lpe_lng(ias);
    hwrk_xax = hwrk_lng(ias);

    CDF_xas = cumsum(mu_xas);
    
    qas_xas  = 189*ones(Nq,1);
    qmpc_xas = 189*ones(Nq,1);
    qlpe_xas = 189*ones(Nq,1);
    qwrk_xas = 189*ones(Nq,1);
    for iq = 1:Nq
        if iq == 1
            inn = (CDF_xas<=qvec(iq));
        else
            inn = (CDF_xas>qvec(iq-1)).*(CDF_xas<=qvec(iq));        
        end
        mm  = mu_xas.*inn;
        xx  = ass.*inn;      qas_xas(iq)  = sum(xx.*mm)/sum(mm);
        xx  = mpc_xas.*inn;  qmpc_xas(iq) = sum(xx.*mm)/sum(mm);
        xx  = lpe_xas.*inn;  qlpe_xas(iq) = sum(xx.*mm)/sum(mm);
        xx  = hwrk_xax.*inn; qwrk_xas(iq) = sum(xx.*mm)/sum(mm);
    end
    
    qash_xas = qas_xas./sum(qas_xas);
    
    qas_xas_xC(:,iC)  = qas_xas;
    qash_xas_xC(:,iC) = qash_xas;
    qmpc_xas_xC(:,iC) = qmpc_xas;
    qlpe_xas_xC(:,iC) = qlpe_xas;
    qwrk_xas_xC(:,iC) = qwrk_xas;

    

    imlb = find(CDF_xas<0.5,1,'last');
    imub = find(CDF_xas>0.5,1,'first'); % mm = mu_xas(imlb:imub);
    a_mdn = sum(ass(imlb:imub).*mu_xas(imlb:imub))/sum(mu_xas(imlb:imub));

    qamnd_xas = qas_xas./a_mdn;

    %---top 10 and top 1
    top10and1 = 189*ones(1,2);
    inn = (CDF_xas>0.90); mm = mu_xas.*inn; xx  = ass.*inn; top10and1(1) = sum(xx.*mm)/sum(ass.*mu_xas);
    inn = (CDF_xas>0.99); mm = mu_xas.*inn; xx  = ass.*inn; top10and1(2) = sum(xx.*mm)/sum(ass.*mu_xas);
    

    %---sort by yl
    % iylp = find(yl_lng)
    [yl_xyl, iyl] = sort(yl_lng);
    mu_xyl  = mu_lng(iyl);
    mpc_xyl = mpc_lng(iyl);
    lpe_xyl = lpe_lng(iyl);

    CDF_xyl = cumsum(mu_xyl);
    
    qyl_xyl  = 189*ones(Nq,1);
    qmpc_xyl = 189*ones(Nq,1);
    qlpe_xyl = 189*ones(Nq,1);
    for iq = 1:Nq
        if iq == 1
            inn = (CDF_xyl<=qvec(iq));
        else
            inn = (CDF_xyl>qvec(iq-1)).*(CDF_xyl<=qvec(iq));        
        end
        mm  = mu_xyl.*inn;
        xx  = yl_xyl.*inn;   qyl_xyl(iq)  = sum(xx.*mm)/sum(mm);
        xx  = mpc_xyl.*inn;  qmpc_xyl(iq) = sum(xx.*mm)/sum(mm);
        xx  = lpe_xyl.*inn;  qlpe_xyl(iq) = sum(xx.*mm)/sum(mm);
    end
    
    qyl_xyl_xC(:,iC)  = qyl_xyl;
    qmpc_xyl_xC(:,iC) = qmpc_xyl;
    qlpe_xyl_xC(:,iC) = qlpe_xyl;
    % qwrk_xyl_xC(:,iC) = qwrk_xyl;
    %---------------------------------------------------------------------------------------------
    %---------------------------------------------------------------------------------------------
    
    
end

Mlpe = -100*qlpe_xyl_xC;

%% data

EMPxFV = [(0.78+0.61)/2, (0.78+0.81)/2, (0.82+0.81)/2, (0.85+0.81)/2, (0.81+0.80)/2];

%% Plot

 ccD = [0.10 0.70 0.70];  

color_rob2   = [0.95 0.75 0.1];
color_rob3   = [4 139 154]/255; %[0.2 0.6 0.4];
color_rob1 = 1.4*[0.6 0.6 0.6]; %[0.75 0.7 0.5];  %[161 149 121]/255; %[0.3 0.2 0.45];
color_rob4   = 1.2*[0.375 0.25 0.5625];

% cc = [0.40 0.40 0.90;...
%       0.85 0.40 0.40;...
%       0.70 0.70 0.70];
  
cc = [color_rob1;...
      color_rob3];
  
fig = figure(11011); clf;
subplot(2,2,1)
iC=1;plot(1:Nq,qwrk_xas_xC(:,iC),'color',cc(iC,:),'LineWidth',3.25,'Marker','o');  hold on
iC=2;plot(1:Nq,qwrk_xas_xC(:,iC),'color',cc(iC,:),'LineWidth',3.25,'Marker','d');

set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
xlabel('Assets quintile','Interpreter','LaTex','Fontsize',24)
ylim([0.30 1.0])
title('Hours worked','Interpreter','LaTex','Fontsize',31)
leg = legend('Benchmark','$(\beta,B)$ model');
legend boxoff
set(leg,'Interpreter','LaTex','Fontsize',23,'Location','SouthWest')
hold off

subplot(2,2,2)
iC=1;plot(1:Nq,qash_xas_xC(:,iC),'color',cc(iC,:),'LineWidth',3.25,'Marker','o');  hold on
iC=2;plot(1:Nq,qash_xas_xC(:,iC),'color',cc(iC,:),'LineWidth',3.25,'Marker','d');
set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
xlabel('Assets quintile','Interpreter','LaTex','Fontsize',24)
title('Assets as $\%$ of total','Interpreter','LaTex','Fontsize',31)
hold off

subplot(2,2,3)
iC=1;plot(qas_xas_xC(:,iC),qmpc_xas_xC(:,iC),'color',cc(iC,:),'LineWidth',3.25,'Marker','o'); hold on
iC=2;plot(qas_xas_xC(:,iC),qmpc_xas_xC(:,iC),'color',cc(iC,:),'LineWidth',3.25,'Marker','d');
set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
xlabel('Assets by quintile','Interpreter','LaTex','Fontsize',24)
title('$mpc$','Interpreter','LaTex','Fontsize',31)
hold off
subplot(2,2,4)
iC=1;plot(qyl_xyl_xC(:,iC),-100*qlpe_xyl_xC(:,iC),'color',cc(iC,:),'LineWidth',3.25,'Marker','o'); hold on
iC=2;plot(qyl_xyl_xC(:,iC),-100*qlpe_xyl_xC(:,iC),'color',cc(iC,:),'LineWidth',3.25,'Marker','d');
set(gca,'XGrid','on','YGrid','on','Fontsize',21)  
xlabel('Labor income by quintile','Interpreter','LaTex','Fontsize',24)
xlim([0 1.85])
title('$lpe^{\tau}$','Interpreter','LaTex','Fontsize',31)
hold off




fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 22 11];
if (save_graph == 1 )
    print('../FIGURES/FIGURE_15','-dpng','-r0')
end