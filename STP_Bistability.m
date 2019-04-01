function out = STP_Bistability(nbpop) 

    warning off ;
    IF_DISP = 0 ;
    IF_SAVE = 0 ;
    options = optimset('Display','off') ; 
 
    function out = InvGauss(t,lbd,mu) %% ISI probability distribution 
        out =  sqrt( lbd / 2 * pi * t.^3 ) * exp( -0.5 * lbd * (t-mu).^2 / t / mu.^2 ) ; 
    end
    
    function out = LaplaceInvGauss(s,lbd,mu)
        out = exp( lbd ./ mu - sqrt( lbd ./ mu.^2 + 2 .* s ) ./ sqrt( 1 ./ lbd ) ) ; 
    end
    
    function out = Pr(U, Tf, Td, lbd, mu)  %% probability of release Pr 
        F = LaplaceInvGauss( 1./Tf , lbd, mu) ; 
        H = LaplaceInvGauss( 1./Td + 1./Tf , lbd, mu) ; 
        D = LaplaceInvGauss( 1./Td , lbd, mu) ; 
        
        y = U .* F ./ ( 1 - ( 1 - U ) .* F ) ; 
        xy = ( 1 - H./F ) .* y ; 
        x = ( 1 - ( 1 + ( 1 - U ) .* xy ) .* D ) ./ ( 1 - ( 1 - U ) .* D ) ; 
        out = U .* x + ( 1 - U ) .* xy ; 
    end
    
    % function out = Pr(U,Tf,Td,lbd,mu)
    %     u_ij = U .* ( 1 + Tf .* mu ) ./ ( 1 + U .* Tf .* mu ) ; 
    %     x_ij = 1 ./ ( 1 + Td .* u_ij .* mu ) ;                 
    %     out = u_ij .* x_ij  ;
    % end
    
    function out = rSofrE(mu)
        out = Iext(2) + ( J(2,1)  - J(2,2) .* J(3,1) ./ J(3,2) .* Pr(U(1),Tf(1),Td(1),lbd,mu) ) .* mu ./ J(2,3) ;
    end

    function out = PCPVSOM(mu,T)
        out = Iext(1) - J(1,3) ./ J(2,3) .* Pr(U(3),Tf(3),Td(3),lbd,rSofrE(mu)) .* Iext(2) + ( J(1,1) - J(1,3) ./ J(2,3) .* J(2,1) .* Pr(U(3),Tf(3),Td(3),lbd,rSofrE(mu)) ) .* mu ...
              +  J(1,3)  ./ J(2,3) .* J(3,1)  ./ J(3,2) .* ( J(2,2) .* Pr(U(3),Tf(3),Td(3),lbd,rSofrE(mu)) - J(1,2) .* J(2,3) ./ J(1,3) ) .* Pr(U(1),T,Td(1),lbd,mu) .* mu ;        
    end

    function out = rIofrE(mu) %% for PC-PV-VIP networks 
        out = ( J(2,3) * Iext(1) - J(1,3) * Iext(2) + ( J(2,3) * J(1,1) - J(1,3) * J(2,1) ) * mu ) / ( J(2,3) * J(1,2) - J(1,3) * J(2,2) ) ; 
    end
    
    switch nbpop
      case 2
        Td = .2 ;
        Tf = .6 ;
        U = .055 ;
        
        dir = 'STP' ; 
        J = [1.92 0.18; 0.25 0.20] ;
        Iext = [1 1] ;
        
        % J = abs( ImportJab('STP',nbpop,dir,0) ) 
        % Iext = ExternalInput('STP',nbpop,dir) ; 
        % A = J(1,2) * J(2,1) / J(1,1) / J(2,2) ; 
        % r0 = 1 / J(2,1) * ( J(2,2) / J(1,2) * Iext(1) - Iext(2) ) ; 
        
        % dir = 'EtoI' ; 
        % J = abs( ImportJab('STP',nbpop,dir,0) ) ; 
        % Iext = ExternalInput('STP',nbpop,dir) ; 
        % A = ( J(2,2) ./ J(1,2) .* Iext(1) - Iext(2) ) ./ J(2,1) ; 
        % r0 = J(1,1) .* J(2,2) / J(1,2) ./ J(2,1) ; 
        
      case 3 
        Td = .2  ; 
        Tf = .5 ; 
        U = .1 ; 
        
        % Td = [.2 0 .2] ; 
        % Tf = [.6 0 .6] ; 
        % U = [.055 0 .055] ; 
        
        dir = 'STP' ; 
        J = abs( ImportJab('STP',nbpop,dir,0) ) ; 
        Iext = ExternalInput('STP',nbpop,dir) ; 
        
        A = J(3,2) / J(3,1) * ( J(2,3) * J(1,1) - J(1,3) * J(2,1) ) / ( J(2,3) * J(1,2) - J(1,3) * J(2,2) ) ; 
        r0 = - ( J(2,3) * Iext(1) - J(1,3) * Iext(2) ) / ( J(2,3) * J(1,2) - J(1,3) * J(2,2) ) ; 
        
        % A = ( J(1,2) * J(2,3) * J(3,1) + J(1,3) * J(3,2) * J(2,1) - J(2,2) * J(1,3) * J(1,1) ) / J(1,1) / J(2,3) / J(3,2) ; 
        % r0 = J(3,2) * ( J(2,3) * Iext(1) - J(1,3) * Iext(2) ) / ( J(1,2) * J(2,3)* J(3,1) + J(1,3) * J(3,2) * J(2,1) - J(2,2) * J(1,3) *J(3,1) ) ; 

    end

    
    lbd = 1 ;
    BISTABLE = 1 ;
    while(BISTABLE<3)
        BALANCE = 0 ;
        while(~BALANCE)
            
            G = rand(nbpop,nbpop) + .1 ;
            G(:,2:nbpop) = - G(:,2:nbpop) ;
            Iext = rand(nbpop,1) + .1 ;

            if(nbpop==3)
                Iext(3) = 0 ;
                G(nbpop,nbpop) = 0 ;
            end

            J = abs(G) ;

            switch nbpop
              case 2
                A = J(1,2) * J(2,1) / J(1,1) / J(2,2) ; 
                r0 = 1 / J(2,1) * ( J(2,2) / J(1,2) * Iext(1) - Iext(2) ) ; 
              case 3
                A = J(3,2) / J(3,1) * ( J(2,3) * J(1,1) - J(1,3) * J(2,1) ) / ( J(2,3) * J(1,2) - J(1,3) * J(2,2) ) ; 
                r0 = - ( J(2,3) * Iext(1) - J(1,3) * Iext(2) ) / ( J(2,3) * J(1,2) - J(1,3) * J(2,2) ) ; 
            end

            BALANCE = A>0 && r0>0 ;

        end
        
        Z = AllZeros(@(r) PrBal(r)-Pr(U(1),Tf(1),Td(1),lbd,1./r),.001,100) ;
        fprintf('Zeros ') 
        fprintf('%.3f ', Z) 
        fprintf('\r') 
        BISTABLE = length(unique(Z)) ;
    end
          
    fprintf('\n')
          
    WriteParam('STP',nbpop,dir,Iext,G,Z) ;
    
    fprintf('Zeros ') 
    fprintf('%.3f ', Z) 
    fprintf('\n') 
    
    fprintf('Iext ')       
    fprintf('%.3f ', Iext )
    fprintf('\n')
    fprintf('Rates ')       
    fprintf('%.3f ', Z)
    fprintf('\r')
    fprintf('\nJij \n')           
    for i=1:nbpop
        for j=1:nbpop
            fprintf('%.3f ', G(i,j) )
        end
        fprintf('\n')
    end
          
    fprintf('A %.3f r0 %.3f \n',A,r0) ; 

    r = 1:50 ; 
    lbd = 1. ; 
    if(IF_DISP)
        plot(r, Pr(U(1),Tf(1),Td(1),lbd,1./r), 'k-' ) ; hold on ; 
    end
    %% from the balanced equation => new expression of Pr 
    
    function out = PrBal(r)        
        out = A .* ( 1 - r0./r ) ; 
    end     

    if(IF_DISP)
            
        if(nbpop==2 || Td(2)==0 ) 
            fprintf('2pop \n')
            plot(r,PrBal(r),'r-') 
        else 
            plot(r,PrBal(r) .* Pr(U(2),Tf(2),Td(2),lbd,1./rIofrE(r) ) ,'r-') 
        end 
        hold off ; 

        xlabel('\nu_E')
        ylabel('P_R')
        ylim([0 .5])

    end

    %%%%%%%%%%%%%%%%%%%%%%%
    %% bifurcation diagram 
    %%%%%%%%%%%%%%%%%%%%%%%

    Sol = [] ; 

    UpState=0 ; 
    X = rand(1) .* ones(3,1) ; 
    
    lbd = 1. ; 
    % T=.01:.01:1 ; 
    T=.1:.05:1 ; 

    for i=1:length(T)
        fprintf('Tf %.3f ', T(i)) 

        pause(.2)
        if(IF_DISP)
            figure(1) ;
            r = .0001:.1:50 ;
            plot(r, Pr(U(1),T(i),Td(1),lbd,1./r), 'k-' ) ; hold on ; 
            plot(r,PrBal(r),'r--') 
            % if(nbpop==2 || Td(2)==0 )
            %     plot(r,PrBal(r),'r--') 
            % else
            %     plot(r, PrBal(r) .* Pr(U(2),Tf(2),Td(2),lbd,1./rIofrE(r) ) ,'r-') 
            % end
            ylim([0 1])
            hold off ; 
        end
 
        Z = AllZeros(@(r) PrBal(r)-Pr(U(1),T(i),Td(1),lbd,1./r),.0001,100) ; 
        % if(nbpop==2)
        %     Z = AllZeros(@(r) PrBal(r)-Pr(U(1),T(i),Td(1),lbd,1./r),.0001,100) ;
        % else
        %     % Z = AllZeros(@(r) PrBal(r) .* Pr(U(2),Tf(2),Td(2),lbd,1./rIofrE(r)) - Pr(U(1),T(i),Td(1),lbd,1./r), .001, 100) ;
        %     Z = AllZeros(@(r) PCPVSOM(r,T(i)), .001, 100) ;
        % end

        if(length(Z)<3)
            if(length(Z)==2)
                X = [Z(1);nan;Z(2)] ;
            else 
                if(UpState) 
                    X = [nan;nan;Z] ;
                else 
                    X = [Z;nan;nan] ;
                end 
            end 
        else 
            X = Z.' ;
            UpState=1 ;
        end
        
        fprintf(' X :') 
        fprintf(' %.3f |', X ) 
        fprintf('\n') ;
        
        Sol = [Sol X] ;

    end
    
    fprintf('\n') ; 

    figname=sprintf('RateEvsTf') ;
    fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

    plot(T,Sol(1,:),'b-') 
    plot(T,Sol(2,:),'b--') 
    plot(T,Sol(3,:),'b-') 
    
    xlabel('T_f') 
    ylabel('Rates') 

    if(IF_SAVE)
        
        figdir = sprintf(['../STP/Figures/%dpop/%s/'],nbpop,dir) ; 
        try
            mkdir(figdir) ; 
        end
        fprintf('Writing to %s%s.svg \n',figdir,figname) ;
        
        file = fullfile(figdir, 'Param.txt') ;
        fileID = fopen(file,'w') ; 

        fprintf(fileID,'Iext ') ;
        fprintf(fileID,'%.3f ',Iext) ;
        fprintf(fileID,'\n') ;

        fprintf(fileID,'Jij \n') ;

        for i=1:nbpop
            for j=1:nbpop
                fprintf(fileID,'%f ', J(i,j) ) ;
            end
            fprintf(fileID,'\n') ; 
        end
        fprintf(fileID,'\n') ;

        fprintf(fileID,'Td ') ;
        fprintf(fileID,'%.3f ',Td) ;
        fprintf(fileID,'\n') ;

        fprintf(fileID,'Tf ') ;
        fprintf(fileID,'%.3f ',Tf) ;
        fprintf(fileID,'\n') ;

        fprintf(fileID,'U ') ;
        fprintf(fileID,'%.3f ',U) ;
        fprintf(fileID,'\n') ;

        fclose(fileID) ;
        
        ProcessFigure(fig, fullfile(figdir,figname)) ;
    end

    function z=AllZeros(f,xmin,xmax,N)
    % Inputs :
    % f : function of one variable
    % [xmin - xmax] : range where f is continuous containing zeros
    % N : control of the minimum distance (xmax-xmin)/N between two zeros
        if (nargin<4)
            N=100;
        end
        dx=(xmax-xmin)/N;
        x2=xmin;
        y2=f(x2);
        z=[];
        for i=1:N
            x1=x2;
            y1=y2;
            x2=xmin+i*dx;
            y2=f(x2);
            if (y1.*y2<=0)                              % Rolle's theorem : one zeros (or more) present
                z=[z,fsolve(f,(x2*y1-x1*y2)/(y1-y2),options)]; % Linear approximation to guess the initial value in the [x1,x2] range.
            end
        end
        if(isempty(z))
            z = nan ;
        end
    end
    
end