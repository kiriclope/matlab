function [u b] = RateInputDist(model,nbpop,dir,Iext,K,g,J,IF_DISP)

    rng('shuffle') ;
    warning off ;

    if(isempty(Iext)) 
        Iext = ExternalInput(model,nbpop,dir) ;
    end
    if( nargin<7 | isempty(J) )
        J = ImportJab(model,nbpop,dir) ;
    end
    if nargin<8
        IF_DISP = false ;
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
    xmax = 1 ;
    x0 = zeros(1,2*nbpop) ;

    % x0(1:nbpop) = -x0max + x0max*rand(1,nbpop) ; 
    % x0(nbpop+1:end) = x0max*rand(1,nbpop) ;

    rd = - xmax + 2*xmax*rand ;
    sig = 2*xmax*rand ;
    x0(1:nbpop) = normrnd(rd,sig,1,nbpop) ;
    x0(nbpop+1:end) = xmax*rand(1,nbpop) ;
     
    if IF_DISP
        fprintf('CI : ')
        fprintf(' u :')
        fprintf(' %.3f |', x0(1:nbpop) )
        fprintf(' b :')
        fprintf(' %.3f |', x0(nbpop+1:end) )
        fprintf(' Error ')
        fprintf('%.3f | ', SelfCstEq(x0).')
        fprintf('\n')
    end

    options = optimset('Display','off') ;
    TOLERANCE = 1E-5 ;
    fval = 0 ;

    while flag<=0 | any(abs(SelfCstEq(x0))>TOLERANCE)
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

        if flag<=0 | any(abs(SelfCstEq(x0))>TOLERANCE) 
            rd = - xmax + 2*xmax*rand ;
            sig = 2*xmax*rand ;
            x0(1:nbpop) = normrnd(rd,sig,1,nbpop) ;
            x0(nbpop+1:end) = xmax*rand(1,nbpop) ;
            % x0(1:nbpop) = -x0max + x0max*rand(1,nbpop) ;
            % x0(nbpop+1:end) = x0max*rand(1,nbpop) ;
        end

        if(ntrial>ntrialmax)
            break ;
        end
        ntrial = ntrial + 1 ;

    end

    if(flag>0)
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
        fprintf('\nNo solution found for u \n')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Utils                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Eq = SelfCstEq(x)
        u = x(1:nbpop) ;
        b = x(nbpop+1:end) ;

        Equ = u.'./sqrt(K) - (Iext.'+ J*QchAvgTF(u,b).') ; 
        Eqb = b.' - J2*QchAvgTF2(u,b).' ;

        Eq = [Equ;Eqb] ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function out = TF(x)
        out = x.*heaviside(x) ;
    end

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = QchAvgTF2(u,b)
        if(b>0)
            out = (u.^2+b).*Phi( u./sqrt(b) ) + u.*sqrt(b).*phi( u./sqrt(b) ) ;
        else
            out = u.^2 ;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

