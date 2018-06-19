function RatesVsCff(model,nbpop,dir,Iext,K,dIlim,IF_DATA,n,g,IF_Nk,IF_RING,Crec,Cff,IF_SAVE,Xlim,Ylim)
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
        Cth = .0625 ;
    else
        Cth = 100 ;
    end

    if strcmpi(model,'Rates') | IF_DATA == 0
        figname=sprintf('RatesVsCff') ;
        % figname=sprintf('RatesvsdI%.2f',Cff) ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

        X = [] ;
        Y = [] ;
        
        I = Iext ;
        for dI=dIlim(1):dIlim(2):dIlim(3)
            I(2) = Iext(2) + dI ;
            [u b] = RateInputDist(model,nbpop,dir,I,K,g,[],false) ;
            
            X = [X dI] ;
            Y = [Y QchAvgTF(u,b).'] ;

            fprintf('dI %.3f ',dI)
            fprintf('Rates ')
            fprintf('%.3f ',QchAvgTF(u,b))
            if(dI==dIlim(1))
                fprintf('\n')
            end
            fprintf('\r')
        end
        fprintf('\n')
        
        for i=1:nbpop
            plot(X,Y(i,:)./Y(i,1),'-','color',cl{i})
        end

        
        xlabel('\delta_{opto}')
        ylabel('Norm. Rates')
        % xlim([.01 100])
        % set(gca,'Xscale', 'log')
        % ylim([.01 10])
        % set(gca,'Yscale', 'log')
        
    end

    if(IF_SAVE & ~IF_DATA)

        figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff) ;
        fprintf('Writing %s \n',figdir)
        
        try
            mkdir(figdir) ;
        end

        ProcessFigure(fig, fullfile(figdir,figname)) ;
        hold off ;
    end

    if IF_DATA
        Iprtr = [dIlim(1) dIlim(2) dIlim(3)] ;
        nbN = nbNeuron(nbpop,n,IF_Nk,[]) ;
        Cpt = CptNeuron(nbpop,nbN) ;
        I = Iprtr(1):Iprtr(2):Iprtr(3) ;
        nb=10 ;

        if IF_OPS            
            OpsIdx = ImportData(model,nbpop,dir,'IdxFile',n,K,g,IF_RING,Crec,Cff,1,ndI,I(1)) ;
            OpsIdx = OpsIdx+1 ;
        end

        for i=1:length(I)+1
                        
            try
                if(i==1)
                    data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,2,1,ndI,Iext(ndI)) ;
                else
                    data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,I(i-1),1,ndI,Cff+Iext(ndI)) ;
                end
                for j=1:length(data(1,:))-1
                    IdvRates(i,j) = mean(data(:,j+1)) ;
                end
            catch
                fprintf('FILE NOT FOUND \n')
                for j=1:n*10000
                    IdvRates(i,j) = nan ;
                end
            end

            fprintf(' Rates ')

            for j=1:nbpop

                X = linspace(-pi,pi,length( IdvRates(i,Cpt(j)+1:Cpt(j+1) ) ) ) ;
                Idv = IdvRates(i,Cpt(j)+1:Cpt(j+1) );

                if(IF_OPS & j==2)

                    Nlength = Cpt(j)+1:Cpt(j+1) ;

                    IdvUp = Idv(OpsIdx) ;
                    IdvDown = Idv(setdiff(Nlength-Cpt(2),OpsIdx) );
                    
                    IdxUp = find(IdvUp>=THRESHOLD) ;
                    IdxDown = find(IdvDown>=THRESHOLD) ;
                    
                    if(i==1)
                        UpSize = length(IdxUp) ;
                        DownSize = length(IdxDown) ;
                    end

                    if( length(IdxUp)==0 )
                        mUp(i) = 0 ;
                    else
                        IdvUp = IdvUp(IdxUp) ;
                        mUp(i) = sum( IdvUp )./UpSize ;
                    end                   

                    fprintf('Up %.3f ', mUp(i) )
                    
                    if( length(IdxDown)==0 )
                        mDown(i) = 0 ;
                    else
                        IdvDown = IdvDown(IdxDown) ;
                        mDown(i) = sum( IdvDown )./DownSize ;
                    end                                                                              
                    fprintf('Down %.3f ', mDown(i) )
                end

                nbROI = length(find(abs(X)<=Cth)) ;
                idx = find(Idv>=THRESHOLD) ;

                if( length(idx)==0 )
                    m(j,i) = 0 ;
                else

                    X = X(idx) ;
                    Idv = Idv(idx) ;
                    idx = find(abs(X)<=Cth) ; % Cff

                    if( length(idx)==0 )
                        m(j,i) = 0 ;
                    else
                        
                        X = X(idx) ;
                        Idv = Idv(idx) ;
                        fprintf(' %d %d ',idx(1),idx(length(idx)))
                        
                        m(j,i) = sum( Idv )./nbROI ;
                    end
                end
                fprintf('%.3f ', m(j,i) )
            end
            
            fprintf('\n')
        end

        for i=1:nbpop
            if(i==1 | i==2)
                figname=sprintf('%s_RatesVsCff_Cth%.3f',dir,2*Cth) ;
                if(IF_OPS)
                    figname=sprintf('%s_RatesvsdI_EI_OPS',dir) ;
                end
            else
                figname=sprintf('%s_RatesvsdI_SV_Cth%.3f',dir,2*Cth) ;
                if(IF_OPS)
                    figname=sprintf('%s_RatesvsdI_SV_OPS',dir) ;
                end
            end
            if( ishandle( findobj('type','figure','name',figname) ) )
                fig = findobj('type','figure','name',figname) ; 
                figure(fig); hold on ; 
            else
                fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
                xlabel('\delta_{opto}')
                ylabel('Norm. Rates')
            end

            if(IF_OPS==0 | i==1)
                plot(I,m(i,2:end)./m(i,1),'o','MarkerEdgeColor',cl{i},'MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1)
                plot(I,m(i,2:end)./m(i,1),'-','Color',cl{i})
            end

            if(IF_OPS & i==2)
                plot(I-Iext(2),mUp./mUp(1),'^','MarkerEdgeColor','m','MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1) 
                plot(I-Iext(2),mUp./mUp(1),'-','Color','m')

                % plot(I-Iext(2),mDown./mDown(1),'v','MarkerEdgeColor','g','MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1)
                % plot(I-Iext(2),mDown./mDown(1),'-','Color','g')
                 
                plot(I-Iext(2),m(i,:)./m(i,1),'o','MarkerEdgeColor','b','MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1) 
                plot(I-Iext(2),m(i,:)./m(i,1),'-','Color','b')
            end

                            
            plot([4*Crec(1) 4*Crec(1)],[0 2],'--r', 'linewidth', 1)
            plot([4*Crec(2) 4*Crec(2)],[0 2],'--b', 'linewidth', 1)
            if(i==2)
                ylim([0 2])
            else
                ylim([0 1.5])
            end

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
            % ylim([.01 10])
            % set(gca,'Yscale', 'log')

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