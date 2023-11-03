clear; close all; clc;

%% load results

cdir = pwd;
odir = '../OUTPUT';

cd(odir)

load mu.txt
load params.txt

load mpe_smallsize.txt
load ysim_smallsize.txt

mpe = mpe_smallsize; ysim = ysim_smallsize;

cd(cdir)


Na = params(1);  Nx = params(2);  Nbeta = params(3);  NS = Na*Nx;


%% reshape

mu_rs  = 189*ones(Na,Nx,Nbeta);
mpe_rs = 189*ones(Na,Nx,Nbeta);
ysim_rs = 189*ones(Na,Nx,Nbeta);
is = 1;
for ix = 1:Nx
for ia = 1:Na
    mu_rs(ia,ix,:)   = mu(is,:);
    mpe_rs(ia,ix,:)  = mpe(is,:);
    ysim_rs(ia,ix,:) = ysim(is,:);
    is = is+1;
end
end


mu_lng    = mu(:,1);
mpe_lng   = mpe(:,1);
ysim_lng  = ysim(:,1);

for ib = 2:Nbeta 
    mu_lng    = [mu_lng;   mu(:,ib)  ];
    mpe_lng   = [mpe_lng;  mpe(:,ib) ];
    ysim_lng  = [ysim_lng; ysim(:,ib)];
end

mpe_lng = -mpe_lng/5;  % Average over 5 years

%% computations


%---sort by assets
[~, iy] = sort(ysim_lng);
mu_xy  = mu_lng(iy);
mpe_xy = mpe_lng(iy);

CDF_xy = cumsum(mu_xy);


qvec = [0.25 0.5 0.75 1.0];  Nq = numel(qvec);
qmpe_xy = 189*ones(Nq,1);
for iq = 1:Nq
    if iq == 1
        inn = (CDF_xy<=qvec(iq));
    else
        inn = (CDF_xy>qvec(iq-1)).*(CDF_xy<=qvec(iq));        
    end
    mm  = mu_xy.*inn;
    xx  = mpe_xy.*inn;  qmpe_xy(iq) = sum(xx.*mm)/sum(mm);
end

display(['MPE small size by inc quart = ', num2str(qmpe_xy')])

%%

cd(odir)

load mpe_mediumsize.txt
load ysim_mediumsize.txt

mpe = mpe_mediumsize; ysim = ysim_mediumsize;

cd(cdir)


%% reshape

mpe_rs = 189*ones(Na,Nx,Nbeta);
ysim_rs = 189*ones(Na,Nx,Nbeta);
is = 1;
for ix = 1:Nx
for ia = 1:Na
    mpe_rs(ia,ix,:)  = mpe(is,:);
    ysim_rs(ia,ix,:) = ysim(is,:);
    is = is+1;
end
end

mpe_lng   = mpe(:,1);
ysim_lng  = ysim(:,1);

for ib = 2:Nbeta 
    mpe_lng   = [mpe_lng;  mpe(:,ib) ];
    ysim_lng  = [ysim_lng; ysim(:,ib)];
end

mpe_lng = -mpe_lng/5;  % Average over 5 years

%% computations


%---sort by assets
[~, iy] = sort(ysim_lng);
mu_xy  = mu_lng(iy);
mpe_xy = mpe_lng(iy);

CDF_xy = cumsum(mu_xy);


qvec = [0.25 0.5 0.75 1.0];  Nq = numel(qvec);
qmpe_xy = 189*ones(Nq,1);
for iq = 1:Nq
    if iq == 1
        inn = (CDF_xy<=qvec(iq));
    else
        inn = (CDF_xy>qvec(iq-1)).*(CDF_xy<=qvec(iq));        
    end
    mm  = mu_xy.*inn;
    xx  = mpe_xy.*inn;  qmpe_xy(iq) = sum(xx.*mm)/sum(mm);
end

display(['MPE medium size by inc quart = ', num2str(qmpe_xy')])
      

%%

cd(odir)

load mpe_largesize.txt
load ysim_largesize.txt

mpe = mpe_largesize; ysim = ysim_largesize;

cd(cdir)


%% reshape

mpe_rs = 189*ones(Na,Nx,Nbeta);
ysim_rs = 189*ones(Na,Nx,Nbeta);
is = 1;
for ix = 1:Nx
for ia = 1:Na
    mpe_rs(ia,ix,:)  = mpe(is,:);
    ysim_rs(ia,ix,:) = ysim(is,:);
    is = is+1;
end
end

mpe_lng   = mpe(:,1);
ysim_lng  = ysim(:,1);

for ib = 2:Nbeta 
    mpe_lng   = [mpe_lng;  mpe(:,ib) ];
    ysim_lng  = [ysim_lng; ysim(:,ib)];
end

mpe_lng = -mpe_lng/5;  % Average over 5 years

%% computations


%---sort by assets
[~, iy] = sort(ysim_lng);
mu_xy  = mu_lng(iy);
mpe_xy = mpe_lng(iy);

CDF_xy = cumsum(mu_xy);


qvec = [0.25 0.5 0.75 1.0];  Nq = numel(qvec);
qmpe_xy = 189*ones(Nq,1);
for iq = 1:Nq
    if iq == 1
        inn = (CDF_xy<=qvec(iq));
    else
        inn = (CDF_xy>qvec(iq-1)).*(CDF_xy<=qvec(iq));        
    end
    mm  = mu_xy.*inn;
    xx  = mpe_xy.*inn;  qmpe_xy(iq) = sum(xx.*mm)/sum(mm);
end

display(['MPE large size by inc quart = ', num2str(qmpe_xy')])
      
      
      