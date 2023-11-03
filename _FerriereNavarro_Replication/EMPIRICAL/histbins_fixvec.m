function pout = histbins_fixvec(xin,xvec)

Nx   = numel(xvec);
pout = zeros(Nx,1);

for ix = 1:Nx
    if ix == 1
        xr   = 0.5*(xvec(1)+xvec(2));        
        xind = (xin<=xr);  
    elseif ix == Nx
        xl   = 0.5*(xvec(Nx-1)+xvec(Nx));        
        xind = (xin>xl);
    else  % ix = 2,...,Nxout-1
        xl = 0.5*(xvec(ix-1)+xvec(ix)  );
        xr = 0.5*(xvec(ix)  +xvec(ix+1));
        xind = (xl<xin).*(xin<=xr);  
    end
    pout(ix) = sum( xind );
end
pout = pout/sum(pout(:));

return