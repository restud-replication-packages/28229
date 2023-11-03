clear; close all; clc;

%% load results

cdir = pwd;
odir = '../OUTPUT';

cd(odir)

load avec.txt
load xvec.txt
load Svec.txt
load mu.txt
load mpc.txt
load lpe.txt
load yls.txt
load params.txt
load hpolw.txt
load tax.txt

cd(cdir)

% Na = numel(avec); Nx = numel(xvec); Nbeta = 4;  NS = Na*Nx;

Na = params(1);  Nx = params(2);  Nbeta = params(3);  NS = Na*Nx;
lbd = tax(1); gma = tax(2); TF = tax(3);


%% reshape

mu_rs  = 189*ones(Na,Nx,Nbeta);
mpc_rs = 189*ones(Na,Nx,Nbeta);
lpe_rs = 189*ones(Na,Nx,Nbeta);
hrs_rs = 189*ones(Na,Nx,Nbeta);
is = 1;
for ix = 1:Nx
for ia = 1:Na
    mu_rs(ia,ix,:)  = mu(is,:);
    mpc_rs(ia,ix,:) = mpc(is,:);
    lpe_rs(ia,ix,:) = lpe(is,:);
    hrs_rs(ia,ix,:) = hpolw(is,:);
    is = is+1;
end
end


tax   = (yls - lbd*yls.^(1-gma)) - TF; 
trate = (1.0 - lbd*yls.^(-gma));


mu_lng  = mu(:,1);
mpc_lng = mpc(:,1);
lpe_lng = lpe(:,1);
hrs_lng = hpolw(:,1);
yl_lng  = yls(:,1);
a_lng   = Svec(:,1);
x_lng   = Svec(:,2);
tax_lng = tax(:,1);
trate_lng = trate(:,1);


for ib = 2:Nbeta 
    mu_lng    = [mu_lng;  mu(:,ib)  ];
    mpc_lng   = [mpc_lng; mpc(:,ib) ];
    lpe_lng   = [lpe_lng; lpe(:,ib) ];
    hrs_lng   = [hrs_lng; hpolw(:,ib) ];
    yl_lng    = [yl_lng;  yls(:,ib) ];
    a_lng     = [a_lng;   Svec(:,1)];
    x_lng     = [x_lng; Svec(:,2)];
    tax_lng   = [tax_lng; tax(:,ib)];
    trate_lng = [trate_lng; trate(:,ib)];
end

%% computations

% Graph measure
% figure
% plot(avec,sum(mu_rs(:,:,1),2)+sum(mu_rs(:,:,2),2)+sum(mu_rs(:,:,3),2),'linewidth',2)

%---sort by mpc
[mpc_xmpc, impc] = sort(mpc_lng);
mu_xmpc  = mu_lng(impc);
lpe_xmpc = lpe_lng(impc);

CDF_xmpc = cumsum(mu_xmpc);

% qvec = [0.2 0.4 0.6 0.8 1.0];  Nq = numel(qvec);
qvec = linspace(0.1,1,10);  Nq = numel(qvec);
qmpc_xmpc = 189*ones(Nq,1);
qlpe_xmpc = 189*ones(Nq,1);
for iq = 1:Nq
    if iq == 1
        inn = (CDF_xmpc<=qvec(iq));
    else
        inn = (CDF_xmpc>qvec(iq-1)).*(CDF_xmpc<=qvec(iq));        
    end
    mm  = mu_xmpc.*inn;
    xx  = mpc_xmpc.*inn;  qmpc_xmpc(iq) = sum(xx.*mm)/sum(mm);
    xx  = lpe_xmpc.*inn;  qlpe_xmpc(iq) = sum(xx.*mm)/sum(mm);
end

%---sort by assets
[ass, ias] = sort(a_lng);
mu_xas  = mu_lng(ias);
mpc_xas = mpc_lng(ias);
lpe_xas = lpe_lng(ias);
hrs_xas = hrs_lng(ias);

CDF_xas = cumsum(mu_xas);

a_med = ass(find(CDF_xas>=0.5,1,'first'));

qvec = [0.2 0.4 0.6 0.8 1.0];  Nq = numel(qvec);
% qvec = linspace(0.1,1,10);  Nq = numel(qvec);
qas_xas  = 189*ones(Nq,1);
qmpc_xas = 189*ones(Nq,1);
qlpe_xas = 189*ones(Nq,1);
qhrs_xas = 189*ones(Nq,1);
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
    xx  = hrs_xas.*inn;  qhrs_xas(iq) = sum(xx.*mm)/sum(mm);
end

qash_xas = qas_xas./sum(qas_xas);

display(['Asset share by asset quintile = ',num2str(qash_xas')])
%qhrs_xas'

%---sort by x
[xss, ixs] = sort(x_lng);
mu_xxs  = mu_lng(ixs);
mpc_xxs = mpc_lng(ixs);
lpe_xxs = lpe_lng(ixs);

CDF_xxs = cumsum(mu_xxs);

qvec = [0.2 0.4 0.6 0.8 1.0];  Nq = numel(qvec);
% qvec = linspace(0.1,1,10);  Nq = numel(qvec);
qxs_xxs  = 189*ones(Nq,1);
qmpc_xxs = 189*ones(Nq,1);
qlpe_xxs = 189*ones(Nq,1);
for iq = 1:Nq
    if iq == 1
        inn = (CDF_xxs<=qvec(iq));
    else
        inn = (CDF_xxs>qvec(iq-1)).*(CDF_xxs<=qvec(iq));        
    end
    mm  = mu_xxs.*inn;
    xx  = xss.*inn;      qxs_xxs(iq)  = sum(xx.*mm)/sum(mm);
    xx  = mpc_xxs.*inn;  qmpc_xxs(iq) = sum(xx.*mm)/sum(mm);
    xx  = lpe_xxs.*inn;  qlpe_xxs(iq) = sum(xx.*mm)/sum(mm);
end

%---sort by yl
% iylp = find(yl_lng)
[yl_xyl, iyl] = sort(yl_lng);
mu_xyl    = mu_lng(iyl);
mpc_xyl   = mpc_lng(iyl);
lpe_xyl   = lpe_lng(iyl);
trate_xyl = trate_lng(iyl);

CDF_xyl = cumsum(mu_xyl);

qvec = [0.2 0.4 0.6 0.8 1.0];  Nq = numel(qvec);
% qvec = linspace(0.1,1,10);  Nq = numel(qvec);
qyl_xyl   = 189*ones(Nq,1);
qmpc_xyl  = 189*ones(Nq,1);
qlpe_xyl  = 189*ones(Nq,1);
qrate_xyl = 189*ones(Nq,1);
for iq = 1:Nq
    if iq == 1
        inn = (CDF_xyl<=qvec(iq));
    else
        inn = (CDF_xyl>qvec(iq-1)).*(CDF_xyl<=qvec(iq));        
    end
    mm  = mu_xyl.*inn;
    xx  = yl_xyl.*inn;    qyl_xyl(iq)   = sum(xx.*mm)/sum(mm);
    xx  = mpc_xyl.*inn;   qmpc_xyl(iq)  = sum(xx.*mm)/sum(mm);
    xx  = lpe_xyl.*inn;   qlpe_xyl(iq)  = sum(xx.*mm)/sum(mm);
    xx  = trate_xyl.*inn; qrate_xyl(iq) = sum(xx.*mm)/sum(mm);    
end


display(['MPC   by asset quint = ', num2str(qmpc_xas')])
display(['LPE^tau by inc quint = ', num2str(-100*qlpe_xyl')])

%% Correlation MPC INCOME

ylat_lng = lbd*(yl_lng).^(1-gma) + TF;

inn  = find(hrs_lng>=0.01);

mux  = mu_lng(inn); mux = mux./sum(mux);
mpcx = mpc_lng(inn);
yyx  = ylat_lng(inn);  % yt_lng  yl_lng  ylat_lng ytat_lng

xx = mpcx.*mux; mpc_mn = sum(xx);  dmpc = mpcx-mpc_mn;  
xx = yyx .*mux; yyx_mn = sum(xx);  dyyx = yyx -yyx_mn;

xx = (dmpc.^2).*mux; mpc_sd = sqrt(sum(xx));
xx = (dyyx.^2).*mux; yyx_sd = sqrt(sum(xx));

xx = dmpc.*dyyx.*mux; cov_x = sum(xx(:));

corr_x = cov_x/(mpc_sd*yyx_sd);

display(['Corr mpc with after-tax labor income = ',num2str(corr_x)])   



yat_lng = lbd*(yl_lng).^(1-gma) + (-1+1.035^0.25)*a_lng + TF; 

mux  = mu_lng(inn); mux = mux./sum(mux);
mpcx = mpc_lng(inn);
yyx  = yat_lng(inn);

xx = mpcx.*mux; mpc_mn = sum(xx);  dmpc = mpcx-mpc_mn;  
xx = yyx .*mux; yyx_mn = sum(xx);  dyyx = yyx -yyx_mn;

xx = (dmpc.^2).*mux; mpc_sd = sqrt(sum(xx));
xx = (dyyx.^2).*mux; yyx_sd = sqrt(sum(xx));

xx = dmpc.*dyyx.*mux; cov_x = sum(xx(:));

corr_x = cov_x/(mpc_sd*yyx_sd);

display(['Corr mpc with after-tax total income = ',num2str(corr_x)])