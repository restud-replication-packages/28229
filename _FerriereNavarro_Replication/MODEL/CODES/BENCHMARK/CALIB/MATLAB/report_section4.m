clear; close all; clc;

%% load results

cdir = pwd;
odir = '../OUTPUT';

cd(odir)

load avec.txt
load xvec.txt
load Svec.txt
load mu.txt
load mpc_persist.txt
load lpe.txt
load yls.txt
load params.txt
load hpolw.txt
load tax.txt


%---Steady-State
load AGGREGATES.txt;  %AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, FC]                       

lbd  = AGGREGATES(1);
G    = AGGREGATES(2);
Cagg = AGGREGATES(3);
Yagg = AGGREGATES(4);
Lagg = AGGREGATES(5);
EMP  = AGGREGATES(6);
Kagg = AGGREGATES(7);
DB   = AGGREGATES(8);
FC   = AGGREGATES(9);
cd(cdir)

mpc = mpc_persist; 

% Na = numel(avec); Nx = numel(xvec); Nbeta = 4;  NS = Na*Nx;

Na = params(1);  Nx = params(2);  Nbeta = params(3);  NS = Na*Nx; NN = NS*Nbeta;
lbd = tax(1); gma = tax(2); TF = tax(3);

dlta  = 0.0235D0; Iagg  = dlta*Kagg;

%% reshape

mu_rs  = 189*ones(Na,Nx,Nbeta); mpc_rs = 189*ones(Na,Nx,Nbeta);
lpe_rs = 189*ones(Na,Nx,Nbeta); hrs_rs = 189*ones(Na,Nx,Nbeta);
is = 1;
for ix = 1:Nx
for ia = 1:Na
    mu_rs(ia,ix,:)  = mu(is,:); mpc_rs(ia,ix,:) = mpc(is,:);
    lpe_rs(ia,ix,:) = lpe(is,:); hrs_rs(ia,ix,:) = hpolw(is,:);
    is = is+1;
end
end

tax   = (yls - lbd*yls.^(1-gma)) - TF; 
trate = (1.0 - lbd*yls.^(-gma));

mu_lng  = mu(:,1); mpc_lng = mpc(:,1);
lpe_lng = lpe(:,1); hrs_lng = hpolw(:,1);
yl_lng  = yls(:,1); a_lng   = Svec(:,1);
x_lng   = Svec(:,2); tax_lng = tax(:,1);
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


%%---sort by yl + use measure mu*yl
wy_lng = yl_lng./sum(yl_lng(:));
py_lng = mu_lng.*yl_lng;  
py_lng = py_lng/sum(py_lng(:));

[yl_xyl, iyl] = sort(yl_lng);
wy_xyl    = wy_lng(iyl);
py_xyl    = py_lng(iyl);
mu_xyl    = mu_lng(iyl);
mpc_xyl   = mpc_lng(iyl);
lpe_xyl   = lpe_lng(iyl);
trate_xyl = trate_lng(iyl);

CDF_xyl = cumsum(mu_xyl);
qvec = [0.2 0.4 0.6 0.8 1.0];  Nq = numel(qvec);   % qvec = linspace(0.1,1,10);  Nq = numel(qvec);


%% Labor 

lpe_mn = sum((-100*lpe_lng).*py_lng);
dlpe   = (-100*lpe_lng)-lpe_mn;

mpc_mn = sum(mpc_lng.*py_lng);
dmpc   = mpc_lng-mpc_mn;

mpc_mn_mu  = sum(mpc_lng.*mu_lng);
dmpc_mu    = mpc_lng - mpc_mn_mu;

% Compute equivalent taxes
Yl = sum(yl_lng.*mu_lng);
trate_mn = sum(trate_lng.*py_lng);

% Normalize tax rate in CG
dtaxCG_boe = 0.01*ones(NN,1);

% Revenues CG
DREVCG = sum(   ( (1+dtaxCG_boe).*(1-dtaxCG_boe.*(-100*lpe_lng)) - 1).* trate_lng .* yl_lng .* mu_lng);
DREVCG_app = dtaxCG_boe(1) * ( trate_mn*Yl- sum((-100*lpe_lng).*trate_lng.*yl_lng.*mu_lng)    )  ;

inx       = find(CDF_xyl>qvec(4),1,'first');
yl_thresh = yl_xyl(inx);
inn       = (yl_lng>=yl_thresh);

% Compute tax rate in MG to generate the same tax revenues 
den = ( sum(trate_lng.*yl_lng.*mu_lng.*inn)- sum((-100*lpe_lng).*trate_lng.*yl_lng.*mu_lng.*inn)    ); 
dtax_g5 =  DREVCG_app/den;

dtaxMG_boe    = (1-inn).*0 + inn.*dtax_g5;
DREVMG_app = dtax_g5 * ( sum(trate_lng.*yl_lng.*mu_lng.*inn)- sum((-100*lpe_lng).*trate_lng.*yl_lng.*mu_lng.*inn)    )  ;

DTCG = dtaxCG_boe.*(1-(-100*lpe_lng)).*trate_lng.*yl_lng;
DTMG = dtaxMG_boe.*(1-(-100*lpe_lng)).*trate_lng.*yl_lng;


% Compute dL in each case

%---dL CG all
dtaxCG_boe_mn = sum(dtaxCG_boe.*py_lng);

xx = dlpe.*(dtaxCG_boe-dtaxCG_boe_mn); 
Cov_lpexdt_CG = sum(xx.*py_lng);

dLCG = - lpe_mn*dtaxCG_boe_mn - Cov_lpexdt_CG;


%---dL MG all
dtaxMG_boe_mn = sum(dtaxMG_boe.*py_lng);

xx = dlpe.*(dtaxMG_boe-dtaxMG_boe_mn); 
Cov_lpexdt_MG = sum(xx.*py_lng);

dLMG = - lpe_mn*dtaxMG_boe_mn - Cov_lpexdt_MG;


display(['CG: change in taxes = ',num2str(100*dtaxCG_boe_mn),'%'])
display(['CG: change in labor total = ',num2str(100*dLCG),'%'])

display(['MG: change in taxes = ',num2str(100*dtaxMG_boe_mn),'%'])
display(['MG: change in labor total = ',num2str(100*dLMG),'%'])
display(['MG: change in labor term1 = ',num2str(-100*lpe_mn*dtaxMG_boe_mn),'%'])
display(['MG: change in labor term2 = +',num2str(-100*Cov_lpexdt_MG),'%'])

display('***')

%% Consumption

% Compute the Tax Burden term in each case

xx = dmpc_mu.*(DTCG-DREVCG_app);
Cov_mpcxDT_mu_CG = sum(xx.*mu_lng);

dCCG_t1 = - mpc_mn_mu.*DREVCG_app - Cov_mpcxDT_mu_CG; 


xx = dmpc_mu.*(DTMG-DREVMG_app);
Cov_mpcxDT_mu_MG = sum(xx.*mu_lng);

dCMG_t1 = - mpc_mn_mu.*DREVMG_app - Cov_mpcxDT_mu_MG; 

% Compute the labor supply channel in each case

xxa = (-100*lpe_lng).*dtaxCG_boe; xxa_mn = sum(xxa.*py_lng);
xx = dmpc.*(xxa-xxa_mn);
Cov_mpcxdtlpe_CG = sum(xx.*py_lng);

dCCG_t2 = - (mpc_mn*(-dLCG) +  Cov_mpcxdtlpe_CG    )*Yl;


xxa = (-100*lpe_lng).*dtaxMG_boe; xxa_mn = sum(xxa.*py_lng);
xx = dmpc.*(xxa-xxa_mn);
Cov_mpcxdtlpe_MG = sum(xx.*py_lng);

dCMG_t2 = - (mpc_mn*(-dLMG) +  Cov_mpcxdtlpe_MG    )*Yl;


dCCG = dCCG_t1 + dCCG_t2; 

dCMG = dCMG_t1 + dCMG_t2; 

% Expressin percentage terms

dCCG_decomp_t1 = ([dCCG;dCCG_t1;-mpc_mn_mu.*DREVCG_app;-Cov_mpcxDT_mu_CG]./Cagg)*100;
dCMG_decomp_t1 = ([dCMG;dCMG_t1;-mpc_mn_mu.*DREVMG_app;-Cov_mpcxDT_mu_MG]./Cagg)*100;


dCCG_decomp_t2 = ([dCCG;dCCG_t2; -(mpc_mn*(-dLCG))*Yl;-Cov_mpcxdtlpe_CG*Yl]./Cagg)*100;
dCMG_decomp_t2 = ([dCMG;dCMG_t2; -(mpc_mn*(-dLMG))*Yl;-Cov_mpcxdtlpe_MG*Yl]./Cagg)*100;


% Report results

%[dCCG_decomp_t1 dCMG_decomp_t1 dCCG_decomp_t2 dCMG_decomp_t2]

display(['CG: change in consumption, tax burden channel = ',num2str(dCCG_decomp_t1(2)),'%'])
display(['MG: change in consumption, tax burden channel = ',num2str(dCMG_decomp_t1(2)),'%'])

display(['Tax burden channel, MG fall relative to CG = ',num2str(dCMG_decomp_t1(2)/dCCG_decomp_t1(2))])

display('***')

display(['CG: change in consumption, total = ',num2str(dCCG_decomp_t1(1)),'%'])
display(['CG: change in consumption, total relative to tax burden = ',num2str(dCCG_decomp_t1(1)/dCCG_decomp_t1(2)),'%'])

display(['MG: change in consumption, total = ',num2str(dCMG_decomp_t1(1)),'%'])
display(['MG: change in consumption, total relative to tax burden = ',num2str(dCMG_decomp_t1(1)/dCMG_decomp_t1(2)),'%'])




