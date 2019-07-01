function m = Space_InputDist(model,nbpop,dir,Iext,K,Crec,Cff,PROFILE,IF_DISP)
        
    warning off ;
    if nargin<9
        IF_DISP=1 ;
        if nargin<8
            PROFILE='RING';
            if nargin<7
                Cff = [] ;
                if nargin<6
                    Crec=[] ;
                    if nargin<5
                        K=1 ;
                        if nargin<4
                            Iext = [] ;
                        end
                    end
                end
            end
        end
    end
    
    J = Import_Jab(model,nbpop,dir) ;
    J2 = J.*J ;
    detJ = det(J) ;

    MFRates = linsolve(J,-Iext.') ;
    MFmodu = (Cff./Crec).*MFRates.' ;

    if IF_DISP
        fprintf('MF Rates :')
        fprintf(' m0 ')
        fprintf('%.3f | ',MFRates)
        fprintf(' m1 ')
        fprintf('%.3f | ',MFmodu)
        fprintf('\n')
    end                                 

    flag = 0 ;
    ntrial = 0 ;
    ntrialmax = 500 ;
    xmax = 2 ;

    x = zeros(1,4.*nbpop) ;
    
    % u0,u1
    % rd = - xmax + 2*xmax*rand ;
    % sig = 2*xmax*rand ;
    % x(1:nbpop) = normrnd(rd,sig,1,nbpop) ;
    x(1:nbpop) = -xmax + 2.*xmax*rand(1,nbpop) ;
    x(nbpop+1:2*nbpop) = xmax*rand(1,nbpop) ;
    % a0,a1
    x(2*nbpop+1:3*nbpop) = xmax*rand(1,nbpop) ;
    x(3*nbpop+1:end) = xmax*rand(1,nbpop) ;
    
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

    while flag<=0 | any(abs(SelfCstEq(x))>TOLERANCE)
        try
            [x,fval,flag] = fsolve(@SelfCstEq,x,options) ;
            % [u,res,fval,flag] = lsqnonlin(@SelfCstEq,u,[0 0],Inf,options) ;
        end

        if IF_DISP
            fprintf('ntrial %d',ntrial)
            fprintf(' flag %d',flag)
            % fprintf(' u0 :')
            % fprintf(' %.3f |', x(1:nbpop))
            % fprintf(' u1 :')
            % fprintf(' %.3f |', x(nbpop+1:2*nbpop))
            % fprintf(' a0 :')
            % fprintf(' %.3f |', x(2*nbpop+1:3*nbpop))
            % fprintf(' a1 :')
            % fprintf(' %.3f |', x(3*nbpop+1:end))
            fprintf(' Error : ')
            fprintf('%.3f | ', SelfCstEq(x))
            fprintf('\r')
        end

        if flag<=0 | any(abs(SelfCstEq(x))>TOLERANCE)
            % rd = - xmax + 2*xmax*rand ;
            % sig = 2*xmax*rand ;
            % x(1:nbpop) = normrnd(rd,sig,1,nbpop) ;
            % u0,u1
            x(1:nbpop) = -xmax + 2.*xmax*rand(1,nbpop) ;
            x(nbpop+1:2*nbpop) = xmax*rand(1,nbpop) ;
            % a0,a1
            x(2*nbpop+1:3*nbpop) = xmax*rand(1,nbpop) ;
            x(3*nbpop+1:end) = xmax*rand(1,nbpop) ;
        end

        if(ntrial>ntrialmax)
            x = NaN(1,2*nbpop) ;
            break ;
        end
        ntrial = ntrial + 1 ;

    end

    if(flag>0)

        u0 = x(1:nbpop) ;
        u1 = x(nbpop+1:2*nbpop) ;

        a0 = x(2*nbpop+1:3*nbpop) ;
        a1 = x(3*nbpop+1:end) ;

        CircVar = 1 - m1(u0,u1,a0,a1) ./ m0(u0,u1,a0,a1) ./ 2 ;

        M0 = m0(u0,u1,a0,a1) ;
        M1 = m1(u0,u1,a0,a1) ;
        
        if IF_DISP            
            fprintf('\n m0 : ')
            fprintf('%.3f ',M0)
            fprintf(' m1 : ')
            fprintf('%.3f ',M1)
            fprintf(' CircVar : ')
            fprintf('%.3f ',CircVar)
            fprintf('\n')
        end

        if IF_DISP
            figname=sprintf('Spatial_Profile') ;
            fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
            cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
            
            phi = -pi./2:.01:pi./2 ;
            for i=1:nbpop
                plot(phi, M0(i)+M1(i).*cos(2.*phi),'-','color',cl{i})
                plot(phi, Iext(i)*(1+Cff(i).*cos(2.*phi)),'--','color',cl{i})
            end
            
            xlabel('\phi')
            ylabel('m(\phi), I^{ext}_i')
            xlim([-pi./2 pi./2])
            ylim([0 1])
            set(gca,'xtick',[-pi./2:pi/4:pi./2],'xticklabel',{'-\pi/2','-\pi/4', '0', '\pi/4', '\pi/2'})
        end
    end
    
    function Eq = SelfCstEq(x)
        u0 = x(1:nbpop) ;
        u1 = x(nbpop+1:2.*nbpop) ;

        a0 = x(2.*nbpop+1:3.*nbpop) ;
        a1 = x(3.*nbpop+1:end) ;

        Equ0 = u0.'./sqrt(K) - ( Iext.' + J*m0(u0,u1,a0,a1).' ) ; 
        Equ1 = u1.'./sqrt(K) - ( Cff.*Iext ).' + J*(Crec.*m1(u0,u1,a0,a1)).' ; 
        
        Eqa0 = a0.' - J2*m0(u0,u1,a0,a1).' ; 
        Eqa1 = a1.' - J2*(Crec.*m1(u0,u1,a0,a1)).' ; 
        
        Eq = [Equ0;Equ1;Eqa0;Eqa1] ; 
    end

    function out = m0(u0,u1,a0,a1)
        u = @(phi) u0 + cos(2.*phi) .*u1 ; 
        a = @(phi) a0 + cos(2.*phi) .*a1 ; 
        out = integral( @(phi) TpsAvgTF(u(phi),a(phi)),0,2*pi,'ArrayValued',true)./pi ;
    end

    function out = m1(u0,u1,a0,a1)
        u = @(phi) u0 + cos(2*phi) .*u1 ; 
        a = @(phi) a0 + cos(2*phi) .*a1 ; 
        out = 2.*integral(@(phi) TpsAvgTF(u(phi),a(phi)).*cos(2.*phi),0,2*pi,'ArrayValued',true)./pi ;
    end

    function out = TpsAvgTF(u,a)
        out = zeros(1,nbpop) ;
        for i=1:nbpop
            if(a(i)>0)
                out(i) = .5.*erfc( ( theta-u(i) )./sqrt(2.*a(i) ) ) ;
            else
                out(i) = 0 ;
            end
        end
    end

end
