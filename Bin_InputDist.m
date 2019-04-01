function [u a b] = Bin_InputDist(nbpop,dir,Iext,K,theta,IF_DISP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [u a b] = Bin_InputDist(nbpop,dir,Iext,K,theta)
% Computes input distribution for binary network in 
% the balance state , u mean input, a total variance,
% b quench variance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rng('shuffle') ;
    warning off ;

    if(isempty(Iext))
        Iext = ExternalInput('Binary',nbpop,dir) ;
    end

    if nargin<5
        theta = 1 ;
    end
    if nargin<6 
        IF_DISP = true ;
    end
    
    J = ImportJab('Binary',nbpop,dir) ;
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
    ntrialmax = 500 ;
    xmax = 1 ;
    x = zeros(1,2.*nbpop) ;
    %% x(1:nbpop) = -xmax + 2.*xmax*rand(1,nbpop) ;
    
    rd = - xmax + 2*xmax*rand ;
    sig = 2*xmax*rand ;
    x(1:nbpop) = normrnd(rd,sig,1,nbpop) ;
    x(nbpop+1:end) = xmax*rand(1,nbpop) ;
    
    if IF_DISP
        fprintf('CI : ')
        fprintf('%.3f | ', x)
        fprintf('Error ')
        fprintf('%.3f | ', SelfCstEq(x).')
        fprintf('\n')
    end
    
    options = optimset('Display','off') ;
    TOLERANCE = 1E-5 ;
    fval = 0 ;
    
    Rates = 0 ;
    while flag<=0 | any(abs(SelfCstEq(x))>TOLERANCE | all(Rates)<=0 )
        try
            [x,fval,flag] = fsolve(@SelfCstEq,x,options) ;
            % [u,res,fval,flag] = lsqnonlin(@SelfCstEq,u,[0 0],Inf,options) ;
        end

        u = x(1:nbpop) ;
        a = x(nbpop+1:end) ;
        Rates = QchAvgTF(u,a) ;

        if IF_DISP
            fprintf('ntrial %d',ntrial)
            fprintf(' flag %d',flag)
            fprintf(' u :')
            fprintf(' %.3f |', x(1:nbpop))
            fprintf(' a :')
            fprintf(' %.3f |', x(nbpop+1:end))
            fprintf(' Error : ')
            fprintf('%.3f | ', SelfCstEq(x))
            fprintf('\r')
        end

        if flag<=0 | any(abs(SelfCstEq(x))>TOLERANCE) | all(Rates)<=0
            rd = - xmax + 2*xmax*rand ;
            sig = 2*xmax*rand ;
            x(1:nbpop) = normrnd(rd,sig,1,nbpop) ;
            x(nbpop+1:end) = xmax*rand(1,nbpop) ;
            % x(1:nbpop) = -xmax + 2.*xmax*rand(1,nbpop) ;
            % x(nbpop+1:end) = xmax*rand(1,nbpop) ;
        end

        if(ntrial>ntrialmax)
            break ;
        end
        ntrial = ntrial + 1 ;

    end

    if(flag>0)
        u = x(1:nbpop) ;
        a = x(nbpop+1:end) ;

        if IF_DISP            
            fprintf('\nu : ')
            fprintf('%.3f | ',u)
            fprintf('\na : ')
            fprintf('%.3f | ',a)
            fprintf('\nRates : ')
            fprintf('%.3f ',QchAvgTF(u,a))
            fprintf('\n')
        end

        %%%%%%%%%%%%%%%%%%%
        % Quench variance %
        %%%%%%%%%%%%%%%%%%%

        ntrial = 0 ;
        bflag = 0 ;
        b = a.*rand(1,nbpop)-eps ;
        [nodes,weights] = GaussHermite_2(50); 

        if IF_DISP
            fprintf('CI : ')
            fprintf('%.3f | ', b)
            fprintf('Error ')
            fprintf('%.3f | ', QchSelfCstEq(b).')
            fprintf('\n')
        end

        options = optimset('Display','off') ;
        % options = optimset('Display','off','Algorithm','levenberg-marquardt') ;

        while bflag<=0 | any(abs(QchSelfCstEq(b))>TOLERANCE)
            try
                [b,fval,bflag] = fsolve(@QchSelfCstEq,b,options) ;
                % [b,res,fval,flag] = lsqnonlin(@QchSelfCstEq,b,[0 0],a,options) ;
            end
            
            if IF_DISP
                fprintf('ntrial %d',ntrial)
                fprintf(' flag %d',bflag)
                fprintf(' out :')
                fprintf(' %.3f |', b)
                fprintf(' Error : ')
                fprintf('%.3f | ', QchSelfCstEq(b))
                fprintf('\r')
            end

            if bflag<=0 | any(abs(QchSelfCstEq(b))>TOLERANCE)
                b = a.*rand(1,nbpop)-eps ;
            end

            if(ntrial>ntrialmax)
                break ;
            end
            ntrial = ntrial + 1 ;
        end

        if(bflag>0 & IF_DISP)
            fprintf('\nb : ')
            fprintf('%.3f | ',b)
            fprintf('\n')
        elseif(IF_DISP)
            fprintf('\nNo solution found for b\n')
        end

    elseif(IF_DISP)
        fprintf('\nNo solution found for u \n')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Utils                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Eq = SelfCstEq(x)
        u = x(1:nbpop) ;
        a = x(nbpop+1:end) ;

        Equ = u.'./sqrt(K) - (Iext.'+ J*QchAvgTF(u,a).') ;
        Eqa = a.' - J2*QchAvgTF(u,a).' ;
        Eq = [Equ;Eqa] ;
    end

    function Eqb = QchSelfCstEq(b)
        Eqb = b.' - J2*QchAvgTF2(u,a,b).' ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = TF(x)
        out = heaviside(x-theta) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = phi(x)
        out = exp(-x.^2./2)./sqrt(2.*pi);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % function out = QchAvgTF(u,a)
    %     if(a>0)
    %         out = .5*erfc( ( theta-u )./sqrt(2.*a) ) ;
    %     else
    %         out = zeros(1,length(a)) ;
    %     end
    % end
    
    function out = QchAvgTF(u,a)
        out = zeros(1,nbpop) ;
        for i=1:nbpop
            if(a(i)>0)
                out(i) = .5*erfc( ( theta-u(i) )./sqrt(2.*a(i) ) ) ;
            else
                out(i) = 0 ;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = QchAvgTF2(u,a,b)
        out = zeros(1,nbpop) ;
        for i =1:nbpop
            try
                integrand=@(x) .25.*erfc( ( theta-u(i)-sqrt(b(i)).*x )./sqrt( 2.*( a(i)-b(i) ) ) ).^2.*phi(x) ;
                out(i) = integral(integrand,-Inf,Inf) ;
            catch
                out(i) = 0 ;
            end
            % integrand=@(x) .25.*erfc( ( theta-u(i)-sqrt(2.*b(i)).*x )./sqrt(2.*( a(i)-b(i) ) ) ).^2 ;
            % out(i) = sum(integrand(nodes).*weights) ;
        end
    end

end

