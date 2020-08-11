function [u a b] = Bin_InputDist(model,nbpop,dir,Iext,K,theta,IF_DISP,u,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [u a b] = Bin_InputDist(nbpop,dir,Iext,K,theta)
% Computes input distribution for binary network in 
% the balance state , u mean input, a total variance,
% b quench variance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rng('shuffle') ;
    warning off ;

    if(isempty(Iext))
        Iext = ExternalInput(model,nbpop,dir) ;
    end

    if nargin<5
        theta = 1 ;
    end
    if nargin<6 
        IF_DISP = true ;
    end

    b=zeros(1,nbpop) ;
     
    J = ImportJab(model,nbpop,dir) ;
    J2 = J.*J ;

    detJ = det(J) ;
    MFRates = linsolve(J,-Iext.') ;

    D = [[1 1 1]; [1 1 1] ;[1 0 0]] ;
    %D = [[1 1 1]; [1 1 1] ;[1/sqrt(K) 1/sqrt(K) 0]] ;
    
    for i=1:nbpop
        for j=1:nbpop
            J(i,j) = J(i,j) * D(i,j) ; 
            J2(i,j) = J2(i,j) * D(i,j) * D(i,j) ; 
        end
    end
    
    if IF_DISP
        fprintf('MF Rates : ')
        fprintf('%.3f | ',MFRates)
        fprintf('\n')
    end

    flag = 0 ;
    ntrial = 0 ;
    ntrialmax = 500 ;
    xmax = 1 ; 
    x0 = zeros(1,2.*nbpop) ;
    %% x(1:nbpop) = -xmax + 2.*xmax*rand(1,nbpop) ;
    
    % XMAX = xmax*rand ;
    % rd = - XMAX + 2*XMAX*rand ;
    % sig = 2*XMAX*rand ;
    % x(1:nbpop) = normrnd(rd,sig,1,nbpop) ;
    % x(nbpop+1:end) = XMAX*rand(1,nbpop) ;
    if(nargin<9) 
            for i=1:nbpop
                XMAX = xmax*rand ;
                rd = - XMAX + 2*XMAX*rand ; 
                sig = XMAX*rand ; 
                x0(i) = normrnd(rd,sig) ;
                x0(i+nbpop) = XMAX*rand ;
            end
    else 
        if( ~isempty(u) )
            x0(1:nbpop) = u ;
            x0(nbpop+1:end) = a ; 
        else
            for i=1:nbpop
                XMAX = xmax*rand ;
                rd = - XMAX + 2*XMAX*rand ; 
                sig = XMAX*rand ; 
                x0(i) = normrnd(rd,sig) ;
                x0(i+nbpop) = XMAX*rand ;
            end
        end
    end
         
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

    Sol = 0 ;
    while( flag<=0 || any(abs(SelfCstEq(x))>TOLERANCE) || all(Sol)<=0 )
        try
            [x,fval,flag] = fsolve(@SelfCstEq,x0,options) ;
            % [u,res,fval,flag] = lsqnonlin(@SelfCstEq,u,[0 0],Inf,options) ;
        end

        u = x(1:nbpop) ;
        a = x(nbpop+1:end) ;
        Sol = QchAvgTF(u,a) ;

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

        if flag<=0 | any(abs(SelfCstEq(x))>TOLERANCE) || all(Sol)<=0
            XMAX = xmax*rand ;
            rd = - XMAX + 2*XMAX*rand ;
            sig = 2*XMAX*rand ;
            x0(1:nbpop) = normrnd(rd,sig,1,nbpop) ;
            x0(nbpop+1:end) = XMAX*rand(1,nbpop) ;
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

        Equ = u.' - sqrt(K) .* (Iext.'+ J*QchAvgTF(u,a).') ;
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
        idx = find(a>0) ;
        nidx= find(a<=0) ;
        out(idx) = .5*erfc( ( theta-u(idx) )./sqrt(2.*a(idx) ) ) ;
        out(nidx) = 0 ;
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

