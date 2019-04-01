function [lbd r] = RateStabilityK(nbpop,K)
    
    model = 'LIF';
    dir = 'L5ChC2' ;
    
    Iext = ExternalInput(model,nbpop,dir) ;
    J = ImportJab(model,nbpop,dir) ;
    
    Tr = zeros(1,nbpop) ;
    Tr(1) = 10 ;
    Tr(2:nbpop) = 10 ;
    Tr(4) = 20 ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Finite N, finite K with quench
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nbN = nbNeuron(nbpop,4,0,[]) ;
    Cpt = CptNeuron(nbpop,nbN) ;
    N = Cpt(nbpop+1) ;

    count = 0 ;
    while count<1
        count = count+1 ;
        try
            [u b] = RateInputDist(model,nbpop,dir,Iext,K,1,J,false) ;
        catch
            fprintf('Error computing inputs\n') 
            u = zeros(1,nbpop) ;
            b = zeros(1,nbpop) ;
        end
    end

    r = QchAvgTF(u,b).' ;
    fprintf('Finite K: Rates ') 
    fprintf('%.3f | ', r)
    fprintf('\n')

    z = normrnd(0,1,1,N) ;
    sympref('HeavisideAtOrigin',0) ;

    for i=1:nbpop
        for j=1:nbpop
            Gain(i,Cpt(j)+1:Cpt(j+1)) = heaviside( u(i) + sqrt(b(i)) .* z(Cpt(j)+1:Cpt(j+1) ) ) ;
        end
        AvgGain(i) = mean( Gain(i,:) ) ;
    end

    Id = eye(nbpop) ;

    for i = 1:nbpop
        for j = 1:nbpop
            G(i,j) = AvgGain(i) * J(i,j) / Tr(i) ;
            Id(i,j) = Id(i,j) / Tr(i) ;
        end
    end
    
    M = ( -Id + sqrt(K) .* G ) ;
    lbd = eig(M) ;

    fprintf('lbd ')
    for i=1:length(M) 
        fprintf('%.3f +i %.3f | ', real(lbd(i)), imag(lbd(i)) ) ; 
    end
    fprintf('\n') ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = phi(x)
        out = exp(-x.^2./2)./sqrt(2.*pi);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = Phi(x)
        out = .5.*(1+erf( x./sqrt(2) ) ) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = QchAvgTF(u,b)
        if(b>0)
            out = u.*Phi( u./sqrt(b) ) + sqrt(b).*phi( u./sqrt(b) ) ;
        else
            out = u ;
        end
    end

end    