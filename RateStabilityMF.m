function [lbd r] = RateStabilityMF(nbpop,K)
    
    model = 'LIF';
    dir = 'L5ChC2' ;
    
    Iext = ExternalInput(model,nbpop,dir) ;
    J = ImportJab(model,nbpop,dir) ;
    
    Tr = zeros(1,nbpop) ;
    Tr(1) = 10 ;
    Tr(2:nbpop) = 10 ;
    Tr(4) = 40 ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MF limit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    options = optimset('Display','off') ;
    r0 = BalRatesMF(model,nbpop,dir,Iext,J,0) ;
    r = fsolve(@SteadyState,r0,options) ;

    fprintf('Rates ') ;
    fprintf('%.3f ', r) ;
    fprintf('\n') ;

    M = linEq(r) ;
    lbd = eig(M) ;

    fprintf('lbd ')
    for i=1:length(M) 
        fprintf('%.3f +i %.3f | ', real(lbd(i)), imag(lbd(i)) ) ; 
    end
    fprintf('\n') ;

    function Eq = SteadyState(r)
        h = zeros(1,nbpop) ;
        Eq = [] ;
        for i=1:nbpop
            h(i) = Iext(i) ; 
            for j=1:nbpop
                h(i) = h(i) + J(i,j) * r(j) ; 
            end
            Eq = [Eq ; r(i) - F( sqrt(K) * h(i) ) ] ; 
        end
    end
    
    function M = linEq(r) 
        for i=1:nbpop
            h(i) = Iext(i) ;
            for j=1:nbpop
                h(i) = h(i) + J(i,j) * r(j) ;
            end
        end
        
        M = zeros(nbpop,nbpop) ;
        
        for i=1:nbpop
            M(i,i) = -1 / Tr(i) / sqrt(K) ;
            for j=1:nbpop
                M(i,j) =  J(i,j) / Tr(i) ;
            end
        end
    end
    
    function out = F(x) 
        if(x>0)
            out = x ;
        else
            out = 0 ;
        end
    end
    
    function out = dF(x)
        if(x>0)
            out = 1 ;
        else
            out = 0 ;
        end
    end

end    