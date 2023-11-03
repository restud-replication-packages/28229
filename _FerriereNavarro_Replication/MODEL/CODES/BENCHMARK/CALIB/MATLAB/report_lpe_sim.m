clear; clc; close all; 


%% Load Results

cd ../OUTPUT

load 'avec.txt'; load 'xvec.txt'; load 'mu.txt'; load Px.txt;

Na = numel(avec); Nx = numel(xvec); Nbeta = 3; Nh = 2;

load 'apol.txt'; load 'cpol.txt'; load 'hpol.txt';

load AGGREGATES.txt; % lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, FC
lbd  = AGGREGATES(1); Gagg = AGGREGATES(2);  Cagg = AGGREGATES(3);  
Yagg = AGGREGATES(4); Lagg = AGGREGATES(5);  EMP = AGGREGATES(6);  
Hagg = AGGREGATES(6)/3; Kagg = AGGREGATES(7);  TF =  0.081*Yagg;

load prices.txt
wge = prices(1); wgeH = 6/7*wge; rk  = prices(2);

hbar = 0.33; hvec = [0,hbar];

Pb = 0*ones(3,3);
Pb(1,2) = 1/200; Pb(2,1) = 0.5/200;
Pb(2,3) = Pb(2,1); Pb(3,2) = Pb(1,2);
Pb(1,1) = 1-1/200; Pb(2,2) = Pb(1,1); Pb(3,3) = Pb(1,1);


%% Reshape steady-state

mu_ss    = 189*ones(Na,Nx,Nbeta);

is = 1;
for ix = 1:Nx
for ia = 1:Na
    mu_ss(ia,ix,:)  = mu(is,:); is = is+1;
end
end

apol_ss  = 189*ones(Na,Nx,Nbeta,Nh);
cpol_ss  = 189*ones(Na,Nx,Nbeta,Nh);
hprob_ss = 189*ones(Na,Nx,Nbeta,Nh);   % a probability
hpol_ss  = 189*ones(Na,Nx,Nbeta,Nh);   % a probability
muw_ss   = 189*ones(Na,Nx,Nbeta,Nh);

aw_ss    = 189*ones(Na,Nx,Nbeta,Nh);
x_ss     = 189*ones(Na,Nx,Nbeta,Nh);

incm_ss  = 189*ones(Na,Nx,Nbeta,Nh);
yl_ss    = 189*ones(Na,Nx,Nbeta,Nh); 
lpol_ss  = 189*ones(Na,Nx,Nbeta,Nh);

b_ss    = 189*ones(Na,Nx,Nbeta,Nh);
a_ss    = 189*ones(Na,Nx,Nbeta,Nh);
h_ss    = 189*ones(Na,Nx,Nbeta,Nh);

is = 1;
for ib = 1:Nbeta
    b_ss(:,:,ib,:)  = ib;
    for ix = 1:Nx
        x_ss(:,ix,ib,:) = (ix);
        for ia=1:Na
            aw_ss(ia,ix,ib,:) = avec(ia);
            a_ss(ia,ix,ib,:)  = ia;
            for ih = 1:Nh
                apol_ss(ia,ix,ib,ih)  = apol(is,ih);
                cpol_ss(ia,ix,ib,ih)  = cpol(is,ih);
                hprob_ss(ia,ix,ib,ih) = hpol(is,ih);
                muw_ss(ia,ix,ib,ih)   = mu_ss(ia,ix,ib)*hprob_ss(ia,ix,ib,ih);
                
                hpol_ss(ia,ix,ib,ih)  = hvec(ih);
                h_ss(ia,ix,ib,ih)  = ih;
                
                lpol_ss(ia,ix,ib,ih)   = hpol_ss(ia,ix,ib,ih).*xvec(ix);
                yl_ss(ia,ix,ib,ih)     = wgeH*lpol_ss(ia,ix,ib,ih);
                incm_ss(ia,ix,ib,ih)   = yl_ss(ia,ix,ib,ih) + rk*aw_ss(ia,ix,ib,ih);
                
            end
            is = is+1;
        end
    end
end


% Reshape Px

Px = reshape(Px(:),Nx,Nx);


%% LPE SIMULATION 

% Simulate N households by (a,x,beta), and YN by (a,x,beta,h)


tauk = 0.35; gamma = 0.1;

Naxh = 25*1000;     % number of households
Ty = 30;            % years of simulation
Tsim = 4*Ty;        % number of quarters


% Make a vector of mu, non-sorted
mu_v = muw_ss(:); mu_cdf = cumsum(mu_v);

% Keep track of (a,x,beta,h) for each point of this vector
a_v  = a_ss(:); x_v  = x_ss(:); b_v  = b_ss(:); h_v  = h_ss(:);

% Random draw
draw = rand([Naxh 1]); asim = 189*ones(Naxh,Tsim);
xsim = asim; hsim = asim; bsim = asim; ysim = asim; ylsim = ysim; yasim = ysim; csim = ysim; csim2 = csim;

for i = 1:Naxh
   index = find(draw(i)<=mu_cdf,1,'first');
   if ((mu_cdf(index+1)-draw(i))<(draw(i)-mu_cdf(index))); index = index + 1; end
   if draw(i)>mu_cdf(end); index = numel(mu_cdf); end
   asim(i,1) = a_v(index);
   xsim(i,1) = x_v(index);
   hsim(i,1) = h_v(index);
   bsim(i,1) = b_v(index);
   
   ylsim(i,1) = wgeH*hvec(hsim(i,1)).*xvec(xsim(i,1));
   ysim(i,1)  = ylsim(i,1) + rk*avec(asim(i,1));

   yasim(i,1) = lbd*( ylsim(i,1))^(1-gamma) + (1-tauk)*rk*avec(asim(i,1)) + TF;% after-tax income
   
   
end


% Now simulate over the time sample.

for t = 2:Tsim

    % Draw two uniform on [0 1], one for x' and one for beta' and one for h
    draw_np = rand([Naxh 3]);

    for i = 1:Naxh
       % Compute a_np
       aaux      = apol_ss(asim(i,t-1),xsim(i,t-1),bsim(i,t-1),hsim(i,t-1));
       asim(i,t) = find(aaux>=avec,1,'last');
       csim(i,t-1) = yasim(i,t-1) + avec(asim(i,t-1)) - aaux ;
       csim2(i,t-1) = cpol_ss(asim(i,t-1),xsim(i,t-1),bsim(i,t-1),hsim(i,t-1));
       
       % Compute x_np
       xaux = draw_np(i,1);
       Pxx = cumsum(Px(xsim(i,t-1),:)); Pxx(end) = Pxx(end)+1e-5;
       xsim(i,t) = find(xaux<=Pxx,1,'first');
       
       % Compute b_np
       baux = draw_np(i,2);
       bsim(i,t) = find(xaux<=cumsum(Pb(bsim(i,t-1),:)),1,'first'); %bsim(i,t-1);
       
       % Compute h_np 
       haux = draw_np(i,3); 
       hpolaux = squeeze(cumsum(hprob_ss(asim(i,t),xsim(i,t),bsim(i,t),:)));
       hsim(i,t) = find(haux<=hpolaux,1,'first');
     
       ylsim(i,t) = wgeH*hvec(hsim(i,t)).*xvec(xsim(i,t));% + r*a_vec(asim(i,1));       

       ysim(i,t)  = (ylsim(i,t)) + rk*avec(asim(i,t));
       yasim(i,t) = lbd*(ylsim(i,t))^(1-gamma) + (1-tauk)*rk*avec(asim(i,t)) + TF;


    end
    
end

 
% Aggregate to annualize

csima = ones(Naxh,Ty); 
hsima  = csima; wsima = hsima; 
yasima = hsima; ylsima = hsima; ysima = yasima;
hwage  = ones(4,1);

for y = 1:Ty
   a = 4*(y-1)+1; b = a+3;
   for i = 1:Naxh
       csima(i,y) = mean(csim(i,a:b));
       for j = 1:4
           if (hvec(hsim(i,a+j-1))>0) 
               hwage(j,1) = (lbd*(wgeH*hbar.*xvec(xsim(i,a+j-1))).^(1-gamma))./hbar;
           else
               %hwage(j,1)=0;
               hwage(j,1) = (lbd*(wgeH*hbar.*xvec(xsim(i,a+j-1))).^(1-gamma))./hbar;

           end
       end
       wsima(i,y) = mean(hwage);
       hsima(i,y) = mean(hvec(hsim(i,a:b)));
       yasima(i,y) = mean(yasim(i,a:b));
       ylsima(i,y) = mean(ylsim(i,a:b));
       ysima(i,y)  = mean(ysim(i,a:b));
   end
    
end


% Vectorize everything
wsima  = wsima(:); hsima  = hsima(:); yasima = yasima(:);
ylsima = ylsima(:); ysima  = ylsima(:); csima  = csima(:);


% Drop non working, keep employed
hindex = (hsima>0);

wsime = wsima(hindex); csime = csima(hindex); hsime = hsima(hindex);

yasime = yasima(hindex); ylsime = ylsima(hindex); ysime  = ysima(hindex);


%% Micro elasticities
% 
% Total

Y = log(hsime); X = [log(wsime) log(csime)]; XX = [ones(length(X),1) X]; beta00 = XX\Y;
disp(['LPE average = ',num2str(beta00(2))])

 
%% Micro elasticities by quintile

qind = [0.2 0.4 0.6 0.8 1];
Nqind = numel(qind);
Nsime = numel(hsime);

y_sorting = wsime;

[xx_sort, ~] = sort(y_sorting);
q_thresh = xx_sort(floor(qind*Nsime));

beta1 = ones(2,Nqind);
beta2 = ones(3,Nqind);

YF = []; YY=[];

a=0;
for iq = 1:Nqind
    b=q_thresh(iq);
    
    hindex = (y_sorting<=b)&(y_sorting>=a);
    
    wsimx = wsime(hindex); csimx = csime(hindex); hsimx = hsime(hindex);    
    
    Y = log(hsimx(:)); X2 = [log(wsimx(:)) log(csimx(:))];
    XX2 = [ones(length(Y),1) X2];
    beta2(:,iq) = XX2\Y;

    a=b;
    
end

% Display results

disp(['LPE per quintile = ',num2str(beta2(2,:))])


