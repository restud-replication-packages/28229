clear; close all; clc;

%% load results

cdir = pwd;
odir = '../OUTPUT';

cd(odir)

load avec.txt
load xvec.txt
load Svec.txt
load mu.txt
load params.txt
load tax.txt

load mpc_annual.txt

mpc = mpc_annual;

load apol.txt

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

% Na = numel(avec); Nx = numel(xvec); Nbeta = 4;  NS = Na*Nx;

Na = params(1);  Nx = params(2);  Nbeta = params(3);  NS = Na*Nx;
lbd = tax(1); gma = tax(2); TF = tax(3);


%% reshape

mu_rs  = 189*ones(Na,Nx,Nbeta);
mpc_rs = 189*ones(Na,Nx,Nbeta);
is = 1;
for ix = 1:Nx
for ia = 1:Na
    mu_rs(ia,ix,:)  = mu(is,:);
    mpc_rs(ia,ix,:) = mpc(is,:);
    is = is+1;
end
end

a_pol_rs = 189*ones(Na,Nx,Nbeta,2);
is = 1;
for ib = 1:Nbeta
    for ix = 1:Nx
        for ia=1:Na
            a_pol_rs(ia,ix,ib,1) = apol(is,1);
            a_pol_rs(ia,ix,ib,2) = apol(is,2);
            is = is+1;
        end
    end
end



mu_lng  = mu(:,1);
mpc_lng = mpc(:,1);
a_lng   = Svec(:,1);

for ib = 2:Nbeta 
    mu_lng    = [mu_lng;  mu(:,ib)  ];
    mpc_lng   = [mpc_lng; mpc(:,ib) ];
    a_lng     = [a_lng;   Svec(:,1)];
end

%% computations


%---sort by assets
[ass, ias] = sort(a_lng);
mu_xas  = mu_lng(ias);
mpc_xas = mpc_lng(ias);

CDF_xas = cumsum(mu_xas);

qvec = [0.2 0.4 0.6 0.8 1.0];  Nq = numel(qvec);
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


display(['MPC by asset quint = ', num2str(qmpc_xas')])



      
      