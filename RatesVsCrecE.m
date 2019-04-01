function RatesVsCrecE(model,nbpop,dir,Iext,K,dIlim,IF_DATA,n,g,IF_Nk,IF_RING,Crec,Cff,IF_SAVE,Xlim,Ylim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function InputDistdI(model,nbpop,dir,Iext,K,dIlim,IF_DATA,n,g,IF_Nk,IF_RING,Crec,Cff,IF_SAVE,Xlim,Ylim)
% Plots population time average rate vs the 
% perturbed input Iext(2)+dI.
% utils : Rate_InputDist(nbpop,dir,I,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF_OPS = 0 ;
    ndI = 2 ;
    
    if(isempty(Iext))
        Iext = ExternalInput(model,nbpop,dir) ;
    end

    warning off ;

    if nargin<7
        IF_DATA = false ;
    end 

    if nargin<10
        IF_Nk = false ;
    end

    if nargin<11
        IF_RING = false ;
    end

    if nargin<14
        IF_SAVE = false ;
    end

    if strcmpi(model,'Rates')
        THRESHOLD = 0 ;
    else
        THRESHOLD = .1 ;
    end

    cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
    % cl = {[1 0 0] [0 0 1] [1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
    
    if( strcmp(IF_RING,'Ring') | strcmp(IF_RING,'Gauss') | strcmp(IF_RING,'Exp') )
        Cth = .125 ;
    else
        Cth = 100 ;
    end

    if IF_DATA
        Iprtr = [dIlim(1) dIlim(2) dIlim(3)] ;
        nbN = nbNeuron(nbpop,n,IF_Nk,[]) ;
        Cpt = CptNeuron(nbpop,nbN) ;
        I = Iprtr(1):Iprtr(2):Iprtr(3) ;
        nb=10 ;

        for i=1:length(I)

            Crec(1) = I(i) ;
                        
            try
                BL = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,1,ndI,Iext(ndI)) ;
                data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,1,ndI,Iext(ndI)+.1);
                for j=1:length(data(1,:))-1
                    IdvRates(i,j) = mean(data(:,j+1)) ;
                    IdvRatesBL(i,j) = mean(BL(:,j+1)) ;
                end
            catch
                fprintf('FILE NOT FOUND \n')
                for j=1:n*10000
                    IdvRates(i,j) = nan ;
                    IdvRatesBL(i,j) = nan ;
                end
            end

            fprintf(' Rates ')

            for j=1:nbpop

                X = linspace(-pi,pi,length( IdvRates(i,Cpt(j)+1:Cpt(j+1) ) ) ) ;
                Idv = IdvRates(i,Cpt(j)+1:Cpt(j+1) );
                IdvBL = IdvRatesBL(i,Cpt(j)+1:Cpt(j+1) );

                nbROI = length(find(abs(X)<=Cth)) ;
                idx = find(Idv>=THRESHOLD) ;

                if( length(idx)==0 )
                    m(j,i) = nan ;
                    mBL(j,i) = nan ;
                else

                    X = X(idx) ;
                    Idv = Idv(idx) ;
                    IdvBL = IdvBL(idx) ;

                    idx = find(abs(X)<=Cth) ; % Cff

                    if( length(idx)==0 )
                        m(j,i) = nan ;
                        mBL(j,i) = nan ;
                    else
                        
                        X = X(idx) ;
                        Idv = Idv(idx) ;
                        IdvBL = IdvBL(idx) ;
                        fprintf(' %d %d ',idx(1),idx(length(idx)))
                        
                        m(j,i) = sum( Idv )./nbROI ;
                        mBL(j,i) = sum( IdvBL )./nbROI ;

                        m(j,i) = m(j,i)./mBL(j,i) ;
                    end
                end
                fprintf('%.3f ', m(j,i) )
                fprintf('%.3f ', mBL(j,i) )
            end
            
            fprintf('\n')
        end

        for i=1:nbpop
            if(i==1 | i==2)
                figname=sprintf('%s_RatesVsCrecE_Cth%.3f',dir,2*Cth) ;
            else
                figname=sprintf('%s_RatesvsdI_SV_Cth%.3f',dir,2*Cth) ;
            end

            if( ishandle( findobj('type','figure','name',figname) ) )
                fig = findobj('type','figure','name',figname) ; 
                figure(fig); hold on ; 
            else
                fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
                xlabel('\sigma_E')
                ylabel('Norm. Rates')
            end
            
            r = m(i,:) ;
            Y = r(~isnan(r)) ;
            X = I(find(~isnan(m(i,:)) )) ;
            plot(X,Y,'o','MarkerEdgeColor',cl{i},'MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1)
            plot(X,Y,'-','Color',cl{i})

                            
            X = linspace(-pi,pi,length( IdvRates(1,Cpt(i)+1:Cpt(i+1) ) ) ) ;
            idx = find(abs(X)<=Cth) ; % Cff

            % for k=1:25
            %     nId = randi([idx(1)+Cpt(i) idx(end)+Cpt(i)]) ;
            %     if IdvRates(1,nId)>THRESHOLD
            %         patchline(I-Iext(2),IdvRates(:,nId)./IdvRates(1,nId),'linestyle','-','edgecolor',cl{i},'edgealpha',.1) 
            %     end
            % end
             
            % xlim([.01 100])
            % set(gca,'Xscale', 'log')

            plot([Cff,Cff],[0 2],'k--', 'linewidth', 1)
            plot([4*Crec(2) 4*Crec(2)],[0 2],'--b', 'linewidth', 1)
            ylim([0 1.5])

            if(IF_SAVE & (i==2 | i==4))
                
                figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff) ;
                fprintf('Writing %s \n',figdir)
                
                try
                    mkdir(figdir) ;
                end
                
                ProcessFigure(fig, fullfile(figdir,figname)) ;
            end
            
        end

        %ylim([0 2])
        
        hold off ;

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

end