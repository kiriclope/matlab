function [M0 M1] = BinRing_InputDist(model,nbpop,dir,Iext,K,Crec,Cff,theta,IF_DISP)
%% Balance => Crec>Cff (Vreeswijk les Houches) 
    warning off ;
    if nargin<8
        theta = 1 ;
    end
    if nargin<9 
        IF_DISP = true ;
    end
    
    IF_NORM = 0 ;

    if(isempty(Iext)) 
        Iext = ExternalInput(model,nbpop,dir)*.1 ; 
    end
       
    J = ImportJab(model,nbpop,dir) ; 
    J2 = J.*J ;

    detJ = det(J) ;
    MFRates = linsolve(J,-Iext.') ;
    
    D = [[1 .025 1]; [1 .025 1] ;[0 .025 0]] ; 

    % D = [[.5 0 1 0]; [0 1 0 0] ; [1 1 0 1]; [0 0 0 .25]] ; 
    
    for i=1:nbpop
        for j=1:nbpop
            G(i,j) = J(i,j) * D(i,j) * Crec(j) ; 
            G2(i,j) = J2(i,j) * D(i,j) * Crec(j) ; 
        end
    end
    
    G0 = G 
    I0 = Cff.*Iext ; 
    Idx = [] ;
    
    IntoSOM = D(3,:) ; 
    IntoSOM(3) = [] ; 

    zeroSOM = find(IntoSOM==0) ; 
    o1toSOM = find(IntoSOM>1/sqrt(K)) ; 
    o0toSOM = intersect(find(IntoSOM>0),find(IntoSOM<=1/sqrt(K)) ) ; 
    
    if( any(zeroSOM) || any(o0toSOM) )
        if( length(o0toSOM)~=length(IntoSOM) )
            G0(:, o0toSOM ) = 0 ; 
        end        
        if( any(o1toSOM) ) 
            G0(:, o1toSOM ) = 0 ; 
        end         
        G0(3,:) = 0 ; 
        Idx = [Idx 3] ; 
    end 
    
    for i=nbpop:-1:1
        if( all(D(i,:)<=1./sqrt(K)) || all( all(G0(i,:)==0)))
            G0(i,:) = [] ;
            I0(i) = [] ; 
            Idx = [Idx i] ; 
        end
    end 

    for i=nbpop:-1:1 
        if(all(G0(:,i)==0))
            G0(:,i) = [] ; 
        end
    end
    
    G0 
    
    modu = linsolve(G0,-I0.') ; 
    MFmodu = modu ;
    if(length(modu)~=nbpop) 
        if(nbpop>3)
            MFmodu(4) = MFmodu(3) ;
        end
        MFmodu(3) = MFmodu(2) ; 
        MFmodu(2) = 0 ; 
    end

    if IF_DISP
        fprintf('MF Rates :')
        fprintf(' m0 ')
        fprintf('%.3f | ',MFRates)
        fprintf(' m1 ')
        fprintf('%.3f | ',MFmodu)
        fprintf('\n')
    end                                 

    if(any(abs(MFmodu./MFRates)>1))
        fprintf('Error: m1>m0 in the MF\n')
    end

    flag = 0 ;
    ntrial = 0 ;
    ntrialmax = 100 ; 
    xmax = 1 ;

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
    
    % if IF_DISP
    %     fprintf('CI : ')
    %     fprintf('%.3f | ', x)
    %     fprintf('Error ')
    %     fprintf('%.3f | ', SelfCstEq(x).')
    %     fprintf('\n')
    % end
    
    options = optimset('Display','off') ;
    TOLERANCE = 1E-5 ; 
    fval = 1 ; 
    flag = -1 ; 

    % [x,fval,flag] = fsolve(@SelfCstEq,x,options) ;

    while flag<=0 % | any(abs(SelfCstEq(x))>TOLERANCE)
        try
            [x,fval,flag] = fsolve(@SelfCstEq,x,options) ;
        end
        % [u,res,fval,flag] = lsqnonlin(@SelfCstEq,u,[0 0],Inf,options) ;

        if IF_DISP
            fprintf('ntrial %d',ntrial)
            fprintf(' flag %d',flag)
            fprintf(' u0 :')
            fprintf(' %.3f |', x(1:nbpop))
            fprintf(' u1 :')
            fprintf(' %.3f |', x(nbpop+1:2*nbpop))
            fprintf(' a0 :')
            fprintf(' %.3f |', x(2*nbpop+1:3*nbpop))
            fprintf(' a1 :')
            fprintf(' %.3f |', x(3*nbpop+1:end))
            fprintf(' Error : ')
            fprintf('%.3f | ', SelfCstEq(x))
            fprintf('\r')
        end

        if flag<=0 || any(abs(SelfCstEq(x))>TOLERANCE)
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
            x = NaN(1,3*nbpop) ;
            flag = -2 ;
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
            if( ishandle( findobj('type','figure','name',figname) ) )
                fig = findobj('type','figure','name',figname) ; 
                figure(fig); hold on ; 
            else
                fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
                xlabel('\phi')
                ylabel('m(\phi)/m0') 
                xlim([-pi./2 pi./2])
                %ylim([0 2])
                set(gca,'xtick',[-pi./2:pi/4:pi./2],'xticklabel',{'-\pi/2','-\pi/4', '0', '\pi/4', '\pi/2'})
            end
            cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
            
            phi = -pi./2:.01:pi./2 ;

            if(~IF_NORM)
                for i=1:nbpop
                    plot(phi,(M0(i)+M1(i).*cos(2.*phi+pi)),'-','color',cl{i})
                    if(i==1 || i==3)
                        plot(phi, (MFRates(i)+MFmodu(i).*cos(2.*phi+pi)),'--','color',cl{i}) 
                    else
                        plot(phi, (MFRates(i)+MFmodu(i).*cos(2.*phi)),'--','color',cl{i}) 
                    end
                end
                ylabel('m(\phi)') 
            else
                for i=1:nbpop
                    plot(phi,(M0(i)+M1(i).*cos(2.*phi+pi))./M0(i),'-','color',cl{i})
                    if(i==1)
                        plot(phi, (MFRates(i)+MFmodu(i).*cos(2.*phi+pi))./MFRates(i),'--','color',cl{i})                 
                    else
                        plot(phi, (MFRates(i)+MFmodu(i).*cos(2.*phi))./MFRates(i),'--','color',cl{i}) 
                    end
                end
            end
        end
    end
    
    function Eq = SelfCstEq(x)
        u0 = x(1:nbpop) ;
        u1 = x(nbpop+1:2.*nbpop) ;

        a0 = x(2.*nbpop+1:3.*nbpop) ;
        a1 = x(3.*nbpop+1:end) ;
        
        Equ0 = u0.'./sqrt(K) - ( Iext.' + J*m0(u0,u1,a0,a1).' ) ;
        Equ1 = u1.'./sqrt(K) - ( (Cff.*Iext).' + G2*m1(u0,u1,a0,a1).' ) ;

        % Equ0 = - ( Iext.' + J*m0(u0,u1,a0,a1).' ) ;
        % Equ1 = - (( Cff.*Iext ).' + G2*(Crec.*m1(u0,u1,a0,a1)).') ;
        %Equ1 = u1.'./sqrt(K) - ( ( Cff.*Iext ).' + G*m1(u0,u1,a0,a1).') ;
        
        Eqa0 = a0.' - J2*m0(u0,u1,a0,a1).' ;
        Eqa1 = a1.' - G2*m1(u0,u1,a0,a1).' ;
        %Eqa1 = a1.' - G2*m1(u0,u1,a0,a1).' ;
        
        Eq = [Equ0;Equ1;Eqa0;Eqa1] ;
    end

    function out = m0(u0,u1,a0,a1)
        u = @(phi) u0 + cos(2.*phi) .*u1 ;
        a = @(phi) a0 + cos(2.*phi) .*a1 ;
        out = integral( @(phi) TpsAvgTF(u(phi),a(phi)),0,pi,'ArrayValued',true)./pi ;
    end

    function out = m1(u0,u1,a0,a1)
        u = @(phi) u0 + cos(2*phi) .*u1 ;
        a = @(phi) a0 + cos(2*phi) .*a1 ;
        out = 2.*integral(@(phi) TpsAvgTF(u(phi),a(phi)).*cos(2.*phi),0,pi,'ArrayValued',true)./pi ;
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
