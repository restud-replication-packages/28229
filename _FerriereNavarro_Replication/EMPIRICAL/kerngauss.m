function  [fout, Fout] = kerngauss(ydata,xvec,hb)

Ny = numel(ydata);
Nx = numel(xvec);

fout = 189*ones(Nx,1);
Fout = 189*ones(Nx,1);

for ix = 1:Nx    
    x = xvec(ix);
    
    %---PDF
    u        = (ydata-x)/hb;    
    phi      = normpdf(u,0,1);    
    fout(ix) = sum(phi(:))/(Ny*hb);
    
    %---CDF
    u        = (x-ydata)/hb;    
    PHI      = normcdf(u,0,1);    
    Fout(ix) = sum(PHI(:))/Ny;
end



return