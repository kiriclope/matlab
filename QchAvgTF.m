function out = QchAvgTF(u,b)

    function out = phi(x)
        out = exp(-x.^2./2)./sqrt(2.*pi);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function out = Phi(x)
        out = .5.*(1+erf( x./sqrt(2) ) ) ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(b>0)
        out = u.*Phi( u./sqrt(b) ) + sqrt(b).*phi( u./sqrt(b) ) ;
    else
        out = abs(u) ;
    end
end
    