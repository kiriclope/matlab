function PEDistdIOLD(model,nbpop,dir,Iext,K,dIlim,IF_DATA,n,g,IF_Nk,IF_RING,Crec,Cff,IF_SAVE,Xlim,Ylim)
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
    
    Cth = .125 ;

    if IF_DATA
        Iprtr = [Iext(2) dIlim(2) dIlim(3)] ;
        nbN = nbNeuron(nbpop,n,IF_Nk,[]) ;
        Cpt = CptNeuron(nbpop,nbN) ;
        I = Iprtr(1):Iprtr(2):Iprtr(3) ;
        nb=10 ;
        
        for i=1
            
            %%%%%%%%%%%%%%%%%%%%
            % Individual Rates %
            %%%%%%%%%%%%%%%%%%%%

            % data = Import_Data('Rate',dir,nbpop,n,K,g,'IdvRates',false,2,I(i),false,0,0,true) ;

            % for j=1:length(data(1,:))
            %     IdvRates(i,j) = mean(data(:,j)) ;
            % end
            
            try
                data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,1,2,I(ndI)) ;
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
            
            % fprintf(' Rates ')
            
            % for j=1:nbpop

            %     X = linspace(-pi,pi,length( IdvRates(i,Cpt(j)+1:Cpt(j+1) ) ) ) ;
            %     Idv = IdvRates(i,Cpt(j)+1:Cpt(j+1) );
            %     nbROI = length(find(abs(X)<=Cth)) ; % Cff

            %     idx = find(Idv>=THRESHOLD) ;

            %     if( length(idx)==0 )
            %         m(j,i) = 0 ;
            %         AvgIdxPE(j,i) = 1 ;
            %         VarIdxPE(j,i) = 0 ;
            %     else

            %         X = X(idx) ;
            %         Idv = Idv(idx) ;
            %         idx = find(abs(X)<=Cth) ; % Cff

            %         if( length(idx)==0 )
            %             m(j,i) = 0 ;
            %             AvgIdxPE(j,i) = 1 ;
            %             VarIdxPE(j,i) = 0 ;
            %         else
                        
            %             X = X(idx) ;
            %             Idv = Idv(idx) ;
            %             fprintf(' %d %d ',idx(1),idx(length(idx)))
                        
            %             m(j,i) = mean( Idv ) ;

            %             dum = IdxPE(i, idx(1)+Cpt(j):idx(end)+Cpt(j) ) ;
            %             dum = dum(~isnan(dum)) ;

            %             AvgIdxPE(j,i) = sum( dum )./nbROI ;
            %             VarIdxPE(j,i) = var( dum ) ;
            %         end
            %         fprintf('%.3f ', m(j,i) )
            %     end
            % end
                         
            % fprintf('\n')
        end

        % for i=1:nbpop
        %     if(i==1 | i==2)
        %         figname=sprintf('PEvsdI_EI') ;
        %     else
        %         figname=sprintf('PEvsdI_SV') ;
        %     end
        %     if( ishandle( findobj('type','figure','name',figname) ) )
        %         fig = findobj('type','figure','name',figname) ; 
        %         figure(fig); hold on ; 
        %     else
        %         fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        %         xlabel('I_{opto}')
        %         ylabel('< PE Idx>')
        %     end

        %     plot(I-Iext(2),AvgIdxPE(i,:),'-o','color',cl{i},'MarkerSize',1,'MarkerFaceColor','none','LineWidth', 1)
            
        %     if(IF_SAVE & (i==2 | i==4))
                
        %         figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff) ;
        %         fprintf('Writing %s \n',figdir)
                
        %         try
        %             mkdir(figdir) ;
        %         end
                
        %         ProcessFigure(fig, fullfile(figdir,figname),1.25) ;
        %     end
        % end

        % for i=1:nbpop
        %     if(i==1 | i==2)
        %         figname=sprintf('varPEvsdI_EI') ;
        %     else
        %         figname=sprintf('varPEvsdI_SV') ;
        %     end
        %     if( ishandle( findobj('type','figure','name',figname) ) )
        %         fig = findobj('type','figure','name',figname) ; 
        %         figure(fig); hold on ; 
        %     else
        %         fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        %         xlabel('I_{opto}')
        %         ylabel('<\delta PE Idx>')
        %     end
            
        %     plot(I-Iext(2),VarIdxPE(i,:),'-o','color',cl{i},'MarkerSize',1,'MarkerFaceColor','none','LineWidth', 1)
                        
        %     if(IF_SAVE & (i==2 | i==4))
                
        %         figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff) ;
        %         fprintf('Writing %s \n',figdir)
                
        %         try
        %             mkdir(figdir) ;
        %         end
                
        %         ProcessFigure(fig, fullfile(figdir,figname),1.25) ;
        %     end

        % end

        List = ['E','I','S','V'] ;
        ndI = 1 ;

        for i=1:nbpop
            figname=sprintf('PEIdxvsdI_%s',List(i)) ;
            fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

            meanIdx = [] ;
            Idx = IdxPE(ndI,Cpt(i)+1:Cpt(i+1)) ;
            idx = find(~isnan(Idx)) ;
            Idx = Idx(idx) ;
            
            X = linspace(-pi,pi,length( IdvRates(ndI,Cpt(j)+1:Cpt(j+1) ) ) ) ;
            X = X(idx) ;

            idx = find(abs(X)<=Cth) ; % Cff

            if( length(idx)==0 )
                Idx = nan ;
            else                
                X = X(idx) ;
                Idx = Idx(idx) ;
                
                meanIdx = [meanIdx mean(Idx) ] ;
                histogram(Idx,50,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i}) ;
                
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