clear; close all; clc;

%% load results

cdir = pwd;
odir = '../OUTPUT';

cd(odir)
load Svec.txt
load mu.txt
load mpc_persist.txt
load params.txt
load yls.txt
load tax.txt
load hpolw.txt

cd(cdir)

mpc = mpc_persist;

Na = params(1);  Nx = params(2);  Nbeta = params(3);  NS = Na*Nx;
lbd = tax(1); gma = tax(2); TF = tax(3);

%% reshape

mu_rs  = 189*ones(Na,Nx,Nbeta);
mpc_rs = 189*ones(Na,Nx,Nbeta);
hrs_rs = 189*ones(Na,Nx,Nbeta);

is = 1;
for ix = 1:Nx
for ia = 1:Na
    mu_rs(ia,ix,:)  = mu(is,:);
    mpc_rs(ia,ix,:) = mpc(is,:);
    hrs_rs(ia,ix,:) = hpolw(is,:);
    is = is+1;
end
end


mu_lng  = mu(:,1);
mpc_lng = mpc(:,1);
a_lng   = Svec(:,1);
yl_lng  = yls(:,1);
hrs_lng = hpolw(:,1);


for ib = 2:Nbeta 
    mu_lng    = [mu_lng;  mu(:,ib)  ];
    mpc_lng   = [mpc_lng; mpc(:,ib) ];
    a_lng     = [a_lng;   Svec(:,1)];
    yl_lng    = [yl_lng;  yls(:,ib) ];
    hrs_lng   = [hrs_lng; hpolw(:,ib) ];

end

%% computations

% Graph measure
% figure
% plot(avec,sum(mu_rs(:,:,1),2)+sum(mu_rs(:,:,2),2)+sum(mu_rs(:,:,3),2),'linewidth',2)


%---sort by assets
[ass, ias] = sort(a_lng);
mu_xas  = mu_lng(ias);
mpc_xas = mpc_lng(ias);

CDF_xas = cumsum(mu_xas);

qvec = [0.2 0.4 0.6 0.8 1.0];  Nq = numel(qvec);
% qvec = linspace(0.1,1,10);  Nq = numel(qvec);
qas_xas  = 189*ones(Nq,1);
qmpc_xas = 189*ones(Nq,1);
for iq = 1:Nq
    if iq == 1
        inn = (CDF_xas<=qvec(iq));
    else
        inn = (CDF_xas>qvec(iq-1)).*(CDF_xas<=qvec(iq));        
    end
    mm  = mu_xas.*inn;
    xx  = ass.*inn;      qas_xas(iq)  = sum(xx.*mm)/sum(mm);
    xx  = mpc_xas.*inn;  qmpc_xas(iq) = sum(xx.*mm)/sum(mm);
end

display(['Persistent MPC by asset quintile = ',num2str(qmpc_xas')])

%% Correlations


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

