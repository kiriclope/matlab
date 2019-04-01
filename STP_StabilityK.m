function [lbd r x] = STP_StabilityK(u,r0)

    nbpop = 1 ; 
    K = 1 ; 

    dir = 'Test' ;
    J = ImportJab('STP',nbpop,dir,0) ;
    Iext = ExternalInput('STP',nbpop,dir) ;

    switch nbpop

      case 1
        Tau = 20 ;
        Trec = 200 ;
        Tfac = 1500 ;
        % u = 0.4 ;
        Iext = r0(1) ;
        J = 1 ;

      case 2
        Tau = [20 10] ;
        Trec = [200 0 ; 0 0] ;
        Tfac = [1500 0 ; 0 0] ;
        u = [.4 1 ; 1 1] ;
        Iext(1) = r0(1) ;

      case 3

        Tau = [20 10 10] ;
        Trec = [ 0 0 0 ; 0 0 0 ; .1 .8 0 ]  ;
        Tfac = [ 0 0 0 ; 0 0 0 ; 1. .3 0] ;
        u = [ 1 1 1 ; 1 1 1 ; .03 .35 1 ] ;

    end
         
    options = optimset('Display','off') ;
    r = fsolve(@SteadyState,rand(1,nbpop),options) ;
    for i=1:nbpop
        for j=1:nbpop
            x(i,j) = 1 / ( 1 + Trec(i,j) * u(i,j) * r(j) ) ;
        end
    end

    fprintf('Rates ') ;
    fprintf('%.3f ', r) ;
    fprintf('\n') ;
    fprintf('x ') ;
    fprintf('%.3f ', x) ;
    fprintf('\n') ;

    M = linEqFast(u,x,r) ;
    M = CheckSTP(M) ;
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
                h(i) = h(i) + J(i,j) * u(i,j) * r(j) / ( 1 + Trec(i,j) * u(i,j) * r(j) ) ; 
            end
            Eq = [Eq ; r(i) - F( sqrt(K) * h(i) ) ] ; 
        end
    end
    
    function M = CheckSTP(M) 
        M( ~any(M,2), : ) = [];  % delete rows of zeros
        M( :, ~any(M,1) ) = [];  %columns
    end
    
    function M = linEqFast(u,x,r)

        for i=1:nbpop
            h(i) = Iext(i) ;
            for j=1:nbpop
                h(i) = h(i) + J(i,j) * u(i,j) * r(j) / ( 1 + Trec(i,j) * u(i,j) * r(j) ) ;
            end
        end

        M = zeros(nbpop*(nbpop+1),nbpop*(nbpop+1)) ; 

        for i=1:nbpop
            M(i,i) = -1 / Tau(i) ; 
            for j=1:nbpop
                if(Trec(i,j)~=0)

                    M(i,j) = M(i,j) + dF(sqrt(K) * h(i) ) * sqrt(K) * J(i,j) * u(i,j) * x(i,j) / Tau(i) ; 
                    M(i,i+j*nbpop) = M(i,i+j*nbpop) + dF( sqrt(K) * h(i)) * sqrt(K) * J(i,j) * u(i,j) * r(j) ;

                    M(i+j*nbpop,i*j+nbpop) = -1 / Trec(i,j) ;
                    M(i+j*nbpop,i) = - u(i,j) * x(i,j) ;
                    M(i+j*nbpop,i*j+nbpop) =  M(i+j*nbpop,i+j*nbpop) - u(i,j) * r(j) ;
                    
                else 
                    M(i,j) = M(i,j) + dF(sqrt(K) * h(i)) * sqrt(K) * J(i,j) / Tau(i) ;
                end

            end
        end
    end
    
    % function out = F(x) 
    %     out = .5*( 1 + erf( x/sqrt(2) ) ) ;
    % end
    
    % function out = dF(x) 
    %     out = exp(-x.^2) / sqrt(2 * pi) ;
    % end

    function out = F(x) 
        a = 1.5 ;
        out = a * log( 1 + exp(x/a) ) ;
    end

    function out = dF(x)
        a = 1.5 ;
        out = exp(x/a) / ( 1 + exp(x/a) ) ;
    end

end