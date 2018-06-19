function CVDistdI(model,nbpop,dir,Iext,K,dIlim,IF_DATA,n,g,IF_Nk,IF_RING,Crec,Cff,IF_SAVE,Xlim,Ylim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function InputDistdI(model,nbpop,dir,Iext,K,dIlim,IF_DATA,n,g,IF_Nk,IF_RING,Crec,Cff,IF_SAVE,Xlim,Ylim)
% Plots population time average rate vs the 
% perturbed input Iext(2)+dI.
% utils : Rate_InputDist(nbpop,dir,I,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
         
    Iprtr = [Iext(2) dIlim(2) dIlim(3)] ;
    nbN = nbNeuron(nbpop,n,IF_Nk,[]) ;
    Cpt = CptNeuron(nbpop,nbN) ;
    I = Iprtr(1):Iprtr(2):Iprtr(3) ;    
    meanCV = [] ;

    for i=1:length(I)
        CumSumSpk=[] ;
        try
            data = ImportData(model,nbpop,dir,'Raster',n,K,g,IF_RING,Crec,Cff,1,2,I(i)) ; 
            Spikes = sortrows(data) ; 
            Spikes(:,2) = Spikes(:,2)./1000. ;
        catch
            fprintf('FILE NOT FOUND\n')
        end

        maxIdx(1) = 1 ;
        
        % length( unique( Spikes(:,1) ) )
        
        for i=1:nbpop
            SpkPop = Spikes( find( Spikes(:,1)>=Cpt(i) & Spikes(:,1)<Cpt(i+1) ) , : ) ;
            excludeSpikes = SpkPop(:,1)<Cpt(i+1) & SpkPop(:,1)>Cpt(i+1)-0.9*nbN(i) ;
            SpkPop(excludeSpikes,:) = [] ; 
            maxIdx(i+1) =  maxIdx(i) + length( unique( SpkPop(:,1) ) ) - 1 ;
        end
        
        for i=1:nbpop
            excludeSpikes = Spikes(:,1)<Cpt(i+1) & Spikes(:,1)>Cpt(i+1)-0.9*nbN(i) ;
            Spikes( excludeSpikes, : ) = [] ;
        end

        % length( unique( Spikes(:,1) ) )
        reSizeIdx = 1:length( unique( Spikes(:,1) ) ) ;
        
        [nbSpk nidx]= hist(Spikes(:,1),unique( Spikes(:,1) ) ) ;
        CumSumSpk = [0 cumsum(nbSpk)] ;
        
        SpkTimes = {} ;
        for j=1:length(CumSumSpk)-1
            fprintf('# %d nbSpk %.3f ',nidx(j),nbSpk(j)) 

            SpkTimes = [ SpkTimes;{Spikes(CumSumSpk(j)+1:CumSumSpk(j)+nbSpk(j)-1,2).'} ] ; 
            ISI = SpkTimes{j}(2:nbSpk(j)-1)-SpkTimes{j}(1:nbSpk(j)-2) ; 
            CV(j) = CVfunc(ISI) ;
             
            fprintf('CV %.3f ',CV(j))
            fprintf('\r')
        end
        
        fprintf('\n')

        fprintf('[CV] ')

        CVI = [] ;
        for k=1:nbpop
            try            
                fprintf('%.3f ', mean( ~isnan( CV( maxIdx(k):maxIdx(k+1) ) ) ) )
                CVI = [CVI ; mean( ~isnan( CV( maxIdx(k):maxIdx(k+1) ) ) ) ] ;
            catch
                CVI = [CVI ; 0 ] ;
            end
        end

        meanCV = [meanCV  CVI] ;
        fprintf('\n')
        
    end

    cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;

    figname=sprintf('CVvsdI') ;
    fig = figure('Name',figname,'NumberTitle','off') ; hold on ;  
    xlabel('I_{opto}')
    ylabel('CV')
    
    for i=1:nbpop
        plot(I(1:end)-Iext(2), meanCV(i,:),'o','MarkerEdgeColor',cl{i},'MarkerSize',5,'MarkerFaceColor','none','LineWidth', 1) 
        plot(I-Iext(2), meanCV(i,:),'-','Color',cl{i}) 
    end
        
    % xlim([.01 100])
    % set(gca,'Xscale', 'log')
    % ylim([.01 10])
    % set(gca,'Yscale', 'log')
    %ylim([0 2])
    
    if(IF_SAVE) 
        figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff) ; 
        fprintf('Writing %s \n',figdir) 
        
        try
            mkdir(figdir) ;
        end
        
        ProcessFigure(fig, fullfile(figdir,figname)) ; 
    end
    hold off ;        

    function out = CVfunc(ISI)
        out = sqrt(var(ISI))./mean(ISI) ;
    end

    function out = CV2func(ISI)
        dISI = abs(ISI(1:end-1)-ISI(2:end)) ;
        sumISI = ISI(1:end-1)+ISI(2:end) ;
        out = 2.* mean( dISI )./mean( sumISI ) ;
    end

end