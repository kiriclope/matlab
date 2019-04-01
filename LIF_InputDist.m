function [u b] = LIF_InputDist(model,nbpop,dir,Iext,K,g,J,IF_DISP,u,b)

    rng('shuffle') ;
    warning off ;

    if(isempty(Iext)) 
        Iext = ExternalInput(model,nbpop,dir) ;
    end
    if( nargin<7 | isempty(J) )
        J = ImportJab(model,nbpop,dir) ; 
    end
    if nargin<8
        IF_DISP = 1 ; 
    end
    
    J = g.*J ; 
    J2 = J.*J ; 
    detJ = det(J) ; 
    MFRates = linsolve(J,-Iext.') ; 

    if IF_DISP
        fprintf('MF Rates : ')
        fprintf('%.3f | ',MFRates)
        fprintf('\n')
    end
    
    flag = 0 ; 
    ntrial = 0 ; 
    ntrialmax = 100 ; 
    if(K==Inf)
        xmax = g ; 
    else
        xmax = sqrt(K)*g ; 
    end

    x0 = zeros(1,3*nbpop) ; 

    for i=1:nbpop
        XMAX = xmax*rand ;
        rd = - XMAX + 2*XMAX*rand ; 
        sig = XMAX*rand ; 
        x0(i) = normrnd(rd,sig) ;
        x0(i+nbpop) = XMAX*rand ;
        x0(i+2*nbpop) = XMAX*rand ;
    end

    if IF_DISP
        fprintf('CI : ')
        fprintf(' u :')
        fprintf(' %.3f |', x0(1:nbpop) )
        fprintf(' a :')
        fprintf(' %.3f |', x0(nbpop+1:2*nbpop) )
        fprintf(' b :')
        fprintf(' %.3f |', x0(2*nbpop+1:3*nbpop) )
        fprintf(' Error ')
        fprintf('%.3f | ', SelfCstEq(x0).')
        fprintf('\n')
    end

    options = optimset('Display','off') ;
    TOLERANCE = 1E-5 ;
    fval = 0 ; 

    while flag<=0 || any(abs(SelfCstEq(x0))>TOLERANCE )
        [x0,fval,flag] = fsolve(@SelfCstEq,x0,options) ;

        if IF_DISP 
            fprintf('ntrial %d',ntrial)
            fprintf(' flag %d',flag)
            fprintf(' u :')
            fprintf(' %.3f |', x0(1:nbpop) )
            fprintf(' b :')
            fprintf(' %.3f |', x0(nbpop+1:end) )
            fprintf(' Error : ')
            fprintf('%.3f | ', SelfCstEq(x0))
            fprintf('\r')
        end

        if flag<=0 || any(abs(SelfCstEq(x0))>TOLERANCE) 

            for i=1:nbpop
                XMAX = xmax*rand ; 
                rd = - XMAX + 2*XMAX*rand ; 
                sig = XMAX*rand ; 
                x0(i) = normrnd(rd,sig) ; 
                x0(i+nbpop) = XMAX*rand ; 
            end
            % x0(1:nbpop) = normrnd(rd,sig,1,nbpop) ; 
            % x0(nbpop+1:end) = XMAX*rand(1,nbpop) ;             
        end

        if( any( QchAvgTF(x0(1:nbpop),x0(nbpop+1:end) )<=0 ) )
            flag = -1 ; 
        end

        if(ntrial>ntrialmax)
            break ;
        end
        ntrial = ntrial + 1 ;

    end

    if(flag>0 && all( QchAvgTF(x0(1:nbpop),x0(nbpop+1:end) )>=0 ) )
        u = x0(1:nbpop) ;
        b = x0(nbpop+1:end) ;
        
        if IF_DISP
            fprintf('\nu : ')
            fprintf('%.3f | ',u)
            fprintf('\nb : ')
            fprintf('%.3f | ',b)
            fprintf('\nRates : ')
            fprintf('%.3f ', QchAvgTF(u,b) )
            fprintf('\n')
        end
    elseif(IF_DISP)
        fprintf('\nNo solution found \n')
        u = NaN(1,nbpop) ;
        b = NaN(1,nbpop) ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Utils                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Eq = SelfCstEq(x)
        u = x(1:nbpop) ;
        a = x(nbpop+1:2*nbpop) ; 
        b = x(2*nbpop+1:3*nbpop) ; 
        
        R = J * QchAvgTF(u,a,b).' ; 

        Equ = u.'./sqrt(K) - ( Iext.' + R ) ; 
        Eqa = a.' - J2 * QchAvgTF2(u,a,b).' ; 
        Eqb = b.' - R + a.' ; 
        
        Eq = [Equ ; Eqa ; Eqb] ; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    function out = phi(x) 
        out = exp(-x.^2./2) ./ sqrt(2.*pi) ; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = QchAvgTF(u,a,b)
        if(all(a)>0 && all(b)>0)
            out = integral( @(x) phi(x) .* Ricciardi(u, sqrt(a) .* x, sqrt(b) ) , -100, 100, 'ArrayValued', true) ; 
        else
            out = 0 ;
        end 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = QchAvgTF2(u,a,b) 
        if(all(a)>0 && all(b)>0)
            out = integral( @(x) phi(x) .* Ricciardi(u, sqrt(a) * x, sqrt(b) ).^2 , -100, 100, 'ArrayValued', true) ; 
        else
            out = 0 ;
        end 
    end
    
    function out = Ricciardi(u,a,b) 
        out = zeros(1,nbpop) ; 
        for i=1:nbpop 
            if(a(i)>0 && b(i)>0) 
                xUp = ( 1 - u(i) - a(i) ) ./ sqrt(b(i)) ; 
                xDn = ( 0 - u(i) - a(i) ) ./ sqrt(b(i)) ; 
                dum = integral( @(x) sqrt(pi) .* exp(x^2) .* ( 1 + erf(x) ) , xDn, xUp , 'ArrayValued', true) ; 
                if(dum~=0 && ~isnan(dum)) 
                    out = 1/dum ; 
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
