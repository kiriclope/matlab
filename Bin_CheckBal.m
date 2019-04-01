function Bin_CheckBal(nbpop,dir,Iext,C0,K,file,n,g,IF_Nk,IF_RING,Crec,Cff,IF_IEXT,nPrtr,IF_PHI,PHI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Bin_CheckBal(nbpop,dir,Iext,K,file,n,g,IF_Nk,IF_PE,nPrtr,Iprtr,IF_RING,Cff,Crec,IF_IEXT)
% Compares simulation to analytics of the balance state
% if : - file=IdvInputs, plots input distributions
%      - file=IdvRates, plots rate distributions
% Utils : - Import_Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(isempty(Iext))
        Iext = ExternalInput('Binary',nbpop,dir) ;
    end
    
    Iprtr=Iext(nPrtr);
    IF_PE=false;

    cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
    if strfind(dir,'FF')
        cl = {[0 0 0] [1 0 0] [0 0 1] [0.7 0.7 0.7]} ;
    end

    theta = 1 ;
    THRESHOLD = 0 ;
    try
        [u a b] = Bin_InputDist(nbpop,dir,Iext.*C0,K,1,false) ;
    catch
        u = 0 ;
        a = 0 ;
        b = 0 ;
    end

    nbN = nbNeuron(nbpop,n,IF_Nk) ;
    Cpt = CptNeuron(nbpop,nbN) ;
    nb=5 ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input Dist %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(file,'IdvInputs')
        data = Import_Data('Binary',dir,nbpop,n,K,g,file,IF_PE,nPrtr,Iprtr,IF_RING,Crec,Cff,IF_IEXT,IF_PHI,PHI) 
        try
            data = Import_Data('Binary',dir,nbpop,n,K,g,file,IF_PE,nPrtr,Iprtr,IF_RING,Crec,Cff,IF_IEXT,IF_PHI,PHI) ;
            for i=2:length(data(1,:))
                IdvInputs(i) = mean(data(:,i)) ;
            end
            tps = data(:,1) ;
            for i=1:nbpop
                for j=1:length(data(:,1))
                    PopInput(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1))) ;
                end
            end
        catch
            fprintf('FILE NOT FOUND \n')
            for i=1:n*10000
                IdvInputs(i) = nan ;
            end
            tps = [] ;
            PopInput = [] ;
        end
                
        fprintf('theory : ')
        fprintf('u ')
        fprintf('%.3f | ', u)
        fprintf('b ')
        fprintf('%.3f | ', b)
        fprintf('\n')
        
        figname=sprintf('InputDist') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

        for i=1:nbpop
            MeanInput(i) = mean( IdvInputs( Cpt(i)+1:Cpt(i+1) ) ) ;
            VarInput(i) = var( IdvInputs( Cpt(i)+1:Cpt(i+1) ) ) ;
            
            h = histogram(IdvInputs( Cpt(i)+1:Cpt(i+1) ),100,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i}) ;
            
            x = min(IdvInputs( Cpt(i)+1:Cpt(i+1) ) ):.1:max( IdvInputs( Cpt(i)+1:Cpt(i+1) ) ) ;
            patchline(x,normpdf(x,u(i),b(i)),'linestyle','-','edgecolor',cl{i},'edgealpha',1) 
        end

        xlabel('I_{syn}')
        ylabel('\rho(m)')

        hold off ;

        fprintf('Simuls : ')
        fprintf('u ')
        fprintf('%.3f | ', MeanInput)
        fprintf('b ')
        fprintf('%.3f | ', VarInput)
        fprintf('\n')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inset PopAvgInput evolution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % size(data)
        % length(tps)
        % size(PopInput)
        
        axes('Position',[.75 .7 .2 .2]) ; hold on ;
        box off
        for i=1:nbpop
            plot(tps,PopInput(i,:),'color',cl{i})
        end    
        xlabel('t (ms)')
        ylabel('<I_{syn}>_i')
        set(gca,'FontSize',8)
        hold off ;

        hold off ;

        axes('Position',[.75 .35 .2 .2]) ; hold on ;
        box off
        for i=1:nbpop
            for k=1:nb
                nId = randi([Cpt(i)+1 Cpt(i+1)]) ;
                patchline(tps,data(:,nId),'linestyle','-','edgecolor',cl{i},'edgealpha',.25) 
            end
        end
        xlabel('t (ms)')
        ylabel('I_{syn_i}')
        set(gca,'FontSize',8)
        hold off ;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rate Dist %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(file,'IdvRates')
        try
            data = Import_Data('Binary',dir,nbpop,n,K,g,file,IF_PE,nPrtr,Iprtr,IF_RING,Crec,Cff,IF_IEXT,IF_PHI,PHI) ;
            for i=2:length(data(1,:))
                IdvRates(i) = mean(data(:,i)) ;
            end
            tps = data(:,1) ;
            for i=1:nbpop
                for j=1:length(data(:,1))
                    PopRate(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1)+1)) ;
                end
            end
        catch
            fprintf('FILE NOT FOUND \n')
            for i=1:n*10000
                IdvRates(i) = nan ;
            end
            tps = [] ;
            PopRate = [] ;
        end
                
        fprintf('theory : ')
        fprintf('%.3f | ', QchAvgTF(u,a) )
        fprintf('\n')
        
        figname=sprintf('RateDist') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        
        for i=1:nbpop
            MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
            VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
            
            m = IdvRates( Cpt(i)+1:Cpt(i+1) ) ;
            m = m(m>THRESHOLD) ;
            % if any(m>=1)
            %     fprintf('!!! m>1 !!! \n')
            % end
            h = histogram(m,100,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i}) ; 
            
            x = min( IdvRates( Cpt(i)+1:Cpt(i+1) ) ):.01:max( IdvRates( Cpt(i)+1:Cpt(i+1) ) );
            % patchline(x,\rho(m)(x,u(i),b(i)),'linestyle','-','edgecolor',cl{i},'edgealpha',1) 
            % sum(i) = trapz(x,rho(x,u(i),a(i),b(i))) ;
            patchline(x,rho(x,u(i),a(i),b(i)),'linestyle','-','edgecolor',cl{i},'edgealpha',1)             
        end
        
        xlabel('m')
        ylabel('\rho(m)')

        hold off ;

        fprintf('Simuls : ')
        fprintf('%.3f | ', MeanRate)
        fprintf('\n')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inset PopAvgInput evolution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % size(data)
        % length(tps)
        % size(PopInput)
        
        axes('Position',[.65 .7 .2 .2]) ; hold on ;
        box off
        for i=1:nbpop
            plot(tps,PopRate(i,:),'color',cl{i})
        end    
        xlabel('t (ms)')
        ylabel('<m>_i')
        set(gca,'FontSize',8)
        hold off ;

        axes('Position',[.65 .35 .2 .2]) ; hold on ;
        box off
        for i=1:nbpop
            for k=1:nb
                nId = randi([Cpt(i)+1 Cpt(i+1)]) ;
                patchline(tps,data(:,nId),'linestyle','-','edgecolor',cl{i},'edgealpha',.25) 
            end
        end
        xlabel('t (ms)')
        ylabel('m_i')
        set(gca,'FontSize',8)
        hold off ;

        axes('Position',[.3 .7 .2 .2]) ; hold on ; 
        box off
        nb = 100 ;
        for i=1:nbpop
        
            Y = [] ;
            k=1 ;
            while k <= nbN(i)-nb+1
                Y = [ Y mean( IdvRates( k+Cpt(i):k+nb+Cpt(i) ) ) ] ;
                k = k + nb ;
            end
            X = linspace(-pi./2,pi./2,nbN(i)./nb) ;
            L = length(X)./2 ;
            patchline(X(1:L),Y(L+1:end),'linestyle','-','edgecolor',cl{i},'edgealpha',.25)
            patchline(X(L+1:end),Y(1:L),'linestyle','-','edgecolor',cl{i},'edgealpha',.25)
            
            xAxis = linspace(0,pi, nbN(i)./nb) ;
            SmoothedRates = SmoothSignal(Y) ;
            [ elmax argmax ] = max(SmoothedRates) ;
            phaseShift = xAxis(argmax) ;        
            % plot(xAxis-pi,SmoothedRates)
            % plot(xAxis,Y,'color',cl{i})
            
            xAxis = linspace(-pi*.5,pi*.5, nbN(i)./nb) ;
            C0 = NthFourierMoment(IdvRates(Cpt(i)+1:Cpt(i+1)),0,phaseShift) ;
            m1 = NthFourierMoment(IdvRates(Cpt(i)+1:Cpt(i+1)),1,phaseShift) ;
            m1alt = M1Component(IdvRates(Cpt(i)+1:Cpt(i+1))) ;
            fprintf('C0 %.3f m1 %.3f m1alt %.3f phaseshift %.3f\n',C0,m1,m1alt,phaseShift)
            
            plot(xAxis, C0+m1.*cos(2.*(xAxis-phaseShift)),'color',cl{i})
            
        end
        
        xlabel('\phi')
        ylabel('<m>_{\phi}')
        set(gca,'FontSize',8)
        hold off ;

    end
    
    %%%%%%%%%
    % Utils %
    %%%%%%%%%

    function out = normpdf(m,u,b)
        out = 1./sqrt(2.*pi.*b).*exp(-(m-u).^2./(2.*b) ) ;
    end

    function out = phi(x)
        out = exp(-x.^2./2)./sqrt(2.*pi);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = Phi(x)
        out = .5.*(1+erf( x./sqrt(2) ) ) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    function out = rho(m,u,a,b) % (b f(ax))^-1 = 1/a f^-1(x/b)
        out = sqrt((a-b)./b).*exp( -(1-u-sqrt(a-b).*sqrt(2).*erfcinv(2.*m) ).^2./(2.*b) ...
                                   + erfcinv(2.*m).^2 ) ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % function out = SmoothSignal(x,y)
    %     out = lowess([x.' y.'],.08) ;
    % end
    
    function out = SmoothSignal(x)
        windowSize = 5 ; 
        Tw = (1./windowSize).*ones(1,windowSize) ; 
        % Tw = ones(1,length(x))/length(x) ;
        out = filter(Tw,.08,x) ;
    end

    function out = M1Component(x)
        dPhi = pi./length(x) ;
        y = 1:length(x) ;
        out = 2.*abs(x*exp(-2j.*y.*dPhi).' )./length(x) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = NthFourierMoment(x, n, phaseShift)
        out = [] ;
        nPoints = length(x) ;
        if(n==0)
            preFactor = 1 ;
        else
            preFactor = 2 * n ;
        end
        y = linspace(0, pi, nPoints) - phaseShift ;
        out = preFactor.*x*cos(2.*n.*y).'./nPoints ;
    end

    %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%
    % Spatial Profile %
    %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%

    if strcmp(file,'Space')
        try
            data = Import_Data('Binary',dir,nbpop,n,K,g,'IdvRates',IF_PE,nPrtr,Iprtr,IF_RING,Crec,Cff,IF_IEXT,IF_PHI,PHI) ;
            for i=1:length(data(1,:))-1
                IdvRates(i) = mean(data(:,i+1)) ;
            end
            tps = data(:,1)./1000 ;
            for i=1:nbpop
                for j=1:length(data(:,1))
                    PopRate(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1))) ;
                end
            end
        catch
            if(strfind(dir,'2pop'))
                for i=1:n*2/3*10000
                    IdvRates(i) = nan ;
                end
            else
                for i=1:n*10000
                    IdvRates(i) = nan ;
                end
            end
            tps = [] ;
            PopRate = [] ;
        end
        
        figname=sprintf('Spatial_Profile') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        % [M0 M1] = Bin_Ring(nbpop,dir,Iext*C0,K,[Crec Crec],[0 Cff],1,0) ;
        M0 = 0 ;        
        RatesK = QchAvgTF(u,a) ;

        for i=1:nbpop
            MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
            VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ; 
            m = IdvRates(1,Cpt(i)+1:Cpt(i+1)) ;

            mshuffle = m ;
            mshuffle(1:length(m)/2) = m(length(m)/2+1:end) ;
            mshuffle(length(m)/2+1:end) = m(1:length(m)/2) ;
            
            idx = find(mshuffle>THRESHOLD) ;
            mshuffle = mshuffle(idx) ;
            
            X = linspace(-pi/2,pi/2,nbN(i)) ;
            % fit = smooth(idx,m,.01,'sgolay') ;
            fit = smooth(idx,mshuffle,.001,'lowess') ;
            X = X(idx) ;
            patchline(X,fit,'linestyle','-','edgecolor',cl{i},'edgealpha',.1) 
            fit = smooth(idx,mshuffle,.05,'lowess') ;
            plot(X,fit,'linewidth',2,'color',cl{i})
            
        end
        
        xlabel('\phi')
        ylabel('[m_i]_t')        
        set(gca,'xtick',[-pi/2:pi/4:pi/2],'xticklabel',{'-\pi/2','-\pi/4', '0', '\pi/4', '\pi/2'})
        xlim([-pi/2,pi/2])
        hold off ;

        MF = Bal_RatesMF_Kinf('Binary',nbpop,dir,Iext*C0,[]) ;
        fprintf('MF : ')
        fprintf('%.3f | ', MF)

        fprintf('Finite K : ')
        fprintf('%.3f | ', QchAvgTF(u,a) )
        fprintf('\n')
        
        fprintf('Ring Finite K : ')
        fprintf('%.3f | ', M0 )

        fprintf('Simuls : ')
        fprintf('%.3f | ', MeanRate)
        fprintf('\n')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(file,'MeanRates')
        try
            data = Import_Data('Binary',dir,nbpop,n,K,g,'IdvRates',IF_PE,nPrtr,Iprtr,IF_RING,Crec,Cff,IF_IEXT,IF_PHI,PHI) ;
            for i=2:length(data(1,:))
                IdvRates(i) = mean(data(:,i)) ;
            end
            tps = data(:,1)./1000 ;
            for i=1:nbpop
                for j=1:length(data(:,1))
                    PopRate(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1)+1)) ;
                end
            end
        catch
            fprintf(' FILE NOT FOUND \n')
            if(strfind(dir,'2pop'))
                for i=1:n*2/3*10000
                    IdvRates(i) = nan ;
                end
            else
                for i=1:n*10000
                    IdvRates(i) = nan ;
                end
            end
            tps = [] ;
            PopRate = [] ;
        end
        
        figname=sprintf('Mean Rates') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

        for i=1:nbpop
            MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1)+1 ) ) ;
            VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1)+1 ) ) ;
            if(i==1 | i==4)
                plot(tps,PopRate(i,:),'--','color',cl{i})
            else
                plot(tps,PopRate(i,:),'color',cl{i})
            end
        end   
        xlabel('t (ms)')
        ylabel('<m>_i')

        drawnow ;
        hold off ;

        MF = Bal_RatesMF_Kinf('Binary',nbpop,dir,Iext*C0,[]) ;
        fprintf('MF : ')
        fprintf('%.3f | ', MF)

        fprintf('Finite K : ')
        fprintf('%.3f | ', QchAvgTF(u,a) )
        fprintf('\n')

        fprintf('Simuls : ')
        fprintf('%.3f | ', MeanRate)
        fprintf('\n')
    end

end