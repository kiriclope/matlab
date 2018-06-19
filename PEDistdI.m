function PEDistdI(model,nbpop,dir,Iext,K,dIlim,IF_DATA,n,g,IF_Nk,IF_RING,Crec,Cff,IF_SAVE,Xlim,Ylim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function InputDistdI(model,nbpop,dir,Iext,K,dIlim,IF_DATA,n,g,IF_Nk,IF_RING,Crec,Cff,IF_SAVE,Xlim,Ylim)
% Plots population time average rate vs the 
% perturbed input Iext(2)+dI.
% utils : Rate_InputDist(nbpop,dir,I,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ndI = 5 ;
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
    
    if( strcmp(IF_RING,'Ring') | strcmp(IF_RING,'Gauss') | strcmp(IF_RING,'Exp') )
        Cth = 125 ;
    else
        Cth = 100 ;
    end

    if IF_DATA
        Iprtr = [Iext(2) dIlim(2) dIlim(3)] ;
        nbN = nbNeuron(nbpop,n,IF_Nk,[]) ;
        Cpt = CptNeuron(nbpop,nbN) ;
        I = Iprtr(1):Iprtr(2):Iprtr(3) ;
        I = [I(1) I(ndI)] ;
                
        %%%%%%%%%%%%%%%%%%%%
        % Individual Rates %
        %%%%%%%%%%%%%%%%%%%%
        
        for i=1:2
        
            try
                data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,1,2,I(i)) ;
                for j=1:length(data(1,:))-1
                    IdvRates(i,j) = mean(data(:,j+1)) ;
                end
            catch
                fprintf('FILE NOT FOUND \n')
                for j=1:n*10000
                    IdvRates(i,j) = nan ;
                end
            end
            
            for j=1:n*10000
                if IdvRates(1,j)>=THRESHOLD
                    IdxPE(i,j) = (IdvRates(1,j)-IdvRates(i,j))./(IdvRates(i,j)+IdvRates(1,j)) ;
                else
                    IdxPE(i,j) = nan ;
                end
            end
            
        end
        
        List = ['E','I','S','V'] ;
        ndI = 2 ;

        for i=1:nbpop
            figname=sprintf('PEIdxvsdI_%s',List(i)) ;
            fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

            meanIdx = [] ;
            Idx = IdxPE(ndI,Cpt(i)+1:Cpt(i+1)) ;
            idx = find(~isnan(Idx)) ;
            Idx = Idx(idx) ;
            
            X = linspace(-pi,pi,length( IdvRates(ndI,Cpt(i)+1:Cpt(i+1) ) ) ) ;
            X = X(idx) ;

            idx = find(abs(X)<=Cth) ; % Cff

            if( length(idx)==0 )
                Idx = nan ;
            else                
                X = X(idx) ;
                Idx = Idx(idx) ;
                
                meanIdx = [meanIdx mean(Idx) ] ;
                if(i==1)
                    histogram(Idx,100,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i}) ;
                    ylim([0 4])
                else
                    histogram(Idx,100,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i}) ;
                    ylim([0 6])
                end
                fprintf('mean Idx ')
                fprintf('%.3f ',meanIdx)
                fprintf('\n')
                
                xlabel('PE index')
                ylabel('pdf')
                xlim([-1 1])
                
                if(IF_SAVE)
                
                    figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff) ;
                    fprintf('Writing %s \n',figdir)
                    
                    try
                        mkdir(figdir) ;
                    end
                    
                    ProcessFigure(fig, fullfile(figdir,figname),1.25) ;
                end
                
                hold off ;        
            end
        end
        
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