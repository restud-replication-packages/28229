cd(xdir)

%%% Transition
load lbd_TR.txt; load gma_TR.txt;
load wgeH_TR.txt;
load G_TR.txt; load DB_TR.txt;

%---transition policies and measure
load hpol_TR.txt; load mu_TR.txt;

cd(cdir)

T_TR = numel(lbd_TR);

%% reshape

mu_t = 189*ones(NS,Nbeta,T_TR);
inn = 1;
for tt = 1:T_TR
    for ib = 1:Nbeta
        ou = inn+NS-1;
        mu_t(:,ib,tt) = mu_TR(inn:ou);                
        inn = ou+1;
    end
end

hpol_t = 189*ones(NS,Nbeta,T_TR);
inn = 1;
for tt = 1:T_TR
    for ib = 1:Nbeta
        ou = inn+NS-1;
        hpol_t(:,ib,tt) = 1 - hpol_TR(inn:ou);       % this si prob of working         
        inn = ou+1;
    end
end

%% computation x t

qtaul_xy_t = 189*ones(Nqind,T_TR);  qyl_xy_t = 189*ones(Nqind,T_TR); mp_t = 189*ones(T_TR,1);

for tt = 1:T_TR
    yl_lng   = 189*ones(NS*Nbeta,1);    
    hp_lng   = 189*ones(NS*Nbeta,1);
    mu_lng   = 189*ones(NS*Nbeta,1);
    muw_lng  = 189*ones(NS*Nbeta,1);
    taul_lng = 189*ones(NS*Nbeta,1);
    
    mp_t(tt) = 0.0;
    inn = 1;
    for ib = 1:Nbeta
        ylx   = wgeH_TR(tt)*(Svec(:,2)*hbar); ylxx = max(ylx,1e-13);
        
        taulSS = (ylx - lbd*(ylx.^(1-gma)))./ylxx;
        taulxt = (ylx - lbd_TR(tt)*(ylx.^(1-gma_TR(tt))))./ylxx;
        
        if CGcase == 1
            taulx = taulxt; 
        else  % MG case
            taulx = (taulxt<taulSS).*taulSS + (taulxt>=taulSS).*taulxt;           
        end
        
        mmx      = mu(:,ib); hhx = hpol_t(:,ib,tt);
        mmp      = (taulxt>taulSS).*mmx.*hhx;  
        mp_t(tt) = mp_t(tt) + sum(mmp);
        
        
        ou = inn + NS-1;
        yl_lng(inn:ou)   = ylx;        
        hp_lng(inn:ou)   = hpol_t(:,ib);
        mu_lng(inn:ou)   = mu_t(:,ib);
        muw_lng(inn:ou)  = mu_lng(inn:ou).*hp_lng(inn:ou);
        taul_lng(inn:ou) = taulx;
        inn = ou+1;
    end    
    mu_lng  = mu_lng  /sum(mu_lng(:));
    muw_lng = muw_lng/sum(muw_lng(:));
    
    %---sort by yl
    [yls,iyls] = sort(yl_lng);
    taul_xyl   = taul_lng(iyls);
    hp_xyl     = hp_lng(iyls);
    mu_xyl     = mu_lng(iyls);
    muw_xyl    = muw_lng(iyls);
    CDF_xyl    = cumsum(muw_xyl);
    
    for iq = 1:Nqind
        if iq == 1
            inn = (CDF_xyl<=qind(iq));
        else
            inn = (CDF_xyl>qind(iq-1)).*(CDF_xyl<=qind(iq));
        end
        mm  = muw_xyl.*inn; 
        xx  = taul_xyl.*inn; qtaul_xy_t(iq,tt) = sum(xx.*mm)/sum(mm);
        xx  = yls.*inn;      qyl_xy_t(iq,tt)   = sum(xx.*mm)/sum(mm);
    end
    
end




