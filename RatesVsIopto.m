clear all ; 
GlobalVars

Iext = ExternalInput(model,nbpop,dir) ;    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

v_Iprtr = v_Iprtr(1):v_Iprtr(2):v_Iprtr(3) ;

J = ImportJab(model,nbpop,dir) ; 
Rates = linsolve(J,-Iext.') ; 

for i=1:length(v_Iprtr)   
    data = ImportData(model, nbpop, dir, 'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, Iext(prtrPop) + v_Iprtr(i) ) ;
    try
        for j=1:length(data(1,:))-1
            IdvRates(i,j) = mean(data(:,j+1)) ;
        end

        fprintf(' Rates ') 
        for j=1:nbpop 
            m(j,i) = 0 ;
            Rates = IdvRates(i, Cpt(j)+1:Cpt(j+1)) ; 
            if(i==1)
                BL = IdvRates(i, Cpt(j)+1:Cpt(j+1)) ; 
            end
            m(j,i) = ratesCutOff(Rates, BL, THRESHOLD, Cth, DIM, L) ; 
            fprintf('%.3f ', m(j,i) ) 
        end 
        fprintf('\n')
        
    catch
        for j=1:nbN(nbpop)
            IdvRates(i,j) = nan ;
        end

        fprintf(' Rates ') 
        for j=1:nbpop 
            m(j,i) = nan ;
            Rates = IdvRates(i, Cpt(j)+1:Cpt(j+1)) ;
            fprintf('%.3f ', m(j,i) )
        end 
        fprintf('\n')
        
    end

end

if(~FIGPERPOP)
    figtitle = sprintf('%s_RatesVsIopto',dir) ;
    if(IF_LOGSCALE)
        figtitle = sprintf('%s_RatesVsIopto_Log',dir) ;
    else
        figtitle = sprintf('%s_RatesVsIopto_LogX',dir) ;
    end
    if(IF_POWER~=0) 
        figtitle = sprintf('%s_%d',figtitle,IF_POWER) ; 
    end
    if( ishandle( findobj('type','figure','name',figtitle) ) )
        fig = findobj('type','figure','name',figtitle) ; 
        fig = figure(fig); hold on ; 
    else
        fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
        xlabel('I_{opto}') 
        if(IF_NORM)
            ylabel('Norm. Rates') 
        else
            ylabel('Rates') 
        end
        if(IF_POWER~=0) 
            xlabel('\Gamma_{opto} (mW/mm^2)') 
        end
    end
end

if(IF_POWER~=0) 
    if IF_POWER==1
        %v_Iprtr = Pinf ./ ( 1 + exp(-( v_Iprtr-I0 )/Iinf ) ) ;
        v_Iprtr = Pinf ./ ( 1 + exp(-( v_Iprtr*2.5 - I0 )/Iinf ) ) ;
        %v_Iprtr = exp( ( sqrt(K)*v_Iprtr*m0 -.4 ) / 0.2253 ) ;
    else
        v_Iprtr = P0 .* ( exp( v_Iprtr .* 20 ./ I0 ) - 1 ) ; 
    end
end

for i=1:nbpop 
    
    if(FIGPERPOP)
        if(i==1 || i==2) 
            figtitle = sprintf('%s_RatesVsIopto%s_EI',dir,popList(prtrPop)) ; 
            if(IF_LOGSCALE)
                figtitle = sprintf('%s_RatesVsIopto%sLogEI',dir,popList(prtrPop)) ;
            end
        else 
            figtitle = sprintf('%s_RatesVsIopto%s_SV',dir,popList(prtrPop)) ; 
            if(IF_LOGSCALE)
                figtitle = sprintf('%s_RatesVsIopto%sLogSV',dir,popList(prtrPop)) ;
            end
        end
        if(IF_POWER~=0) 
            figtitle = sprintf('%s_%d',figtitle,IF_POWER) ; 
            xlabel('\Gamma_{opto} (mW/mm^2)') 
        end
        fig = popPerFig(i,dir,figtitle) ;
        xlabel('I_{opto}') 
        if(IF_NORM)
            ylabel('Norm. Rates') 
        else 
            ylabel('Rates') 
        end 
        if(IF_POWER~=0) 
            xlabel('\Gamma_{opto} (mW/mm^2)') 
        end
    end
    
    if(IF_NORM)
        NormRates = m(i,:)./m(i,1) ;
    else
        NormRates = m(i,:) ;
    end

    id = find( ~isnan(NormRates) ) ; 
    id = IDX:length(v_Iprtr) ; 

    if(~IF_COUNTPIF && ~IF_ROBUST)
        plot(v_Iprtr(id), NormRates(id), '-','Color',cl{i}) 
        plot(v_Iprtr(id), NormRates(id), 'o','MarkerEdgeColor',cl{i},'MarkerSize',1,'MarkerFaceColor','none','LineWidth', 1) 
    elseif(IF_ROBUST)
        patchline(v_Iprtr(id), NormRates(id), ...
                  'linestyle','-','edgecolor',cl{i},'edgealpha',alp,'linewidth',1.5) 
        plot(v_Iprtr(id), NormRates(id), mk,'MarkerEdgeColor',cl{i},'MarkerSize',1,'MarkerFaceColor','none','LineWidth', 1)
    end
    
    if( (i==2 || i==4 ) && IF_NORM)
        plot(v_Iprtr(IDX:end),ones(1, length(v_Iprtr(IDX:end)) ),'--','Color','k') 

    %     X=[ 0.0955    0.1592    0.3183    0.4775    0.6366    1.0504    1.5915    2.5465    4.7746] ;

    %     if( strfind(dir,'L23') ) 
    %         YE  = [0.7606    0.6740    0.4651    0.3539    0.2710    0.1742    0.1076    0.0743    0.0615] ;
    %         YI = [1.2566    1.2614    1.4110    1.6964    2.0685    2.8008    3.1163    4.2371    4.9344] ;
    %     else
    %         YE=[  0.8850    0.8073    0.5802    0.4504    0.3465    0.1590    0.0874    0.0512    0.0306] ;
    %         YI=[  0.7772    0.7402    0.6295    0.6209    0.6190    0.6844    0.9313    1.4011    1.7021] ;
    %     end
    %     if(i==2)
    %         plot(X,YE,'r--')
    %         plot(X,YI,'b--')
    %     end
    end
    
    if(IF_IDVTraces) 
        countUp = 1 ;
        countDown = 1 ;
        countMax = 1 ;
        while countUp+countDown<nbIdv && countMax<100 
            nId = randi([Cpt(i)+1 Cpt(i+1)]) ; 
            countMax = countMax + 1 ;
            
            if IdvRates(1,nId)>THRESHOLD
                
                if IdvRates(4,nId)./IdvRates(1,nId) >= 1 && countUp<=nbIdv/2
                    countUp = countUp+1 ;
                    IdvNormRates = IdvRates(:,nId)./IdvRates(1,nId) ;
                    patchline(v_Iprtr(id), IdvNormRates(id), ...
                              'linestyle','-','edgecolor',cl{i},'edgealpha',.2,'linewidth',1.5) 
                    
                end                    
                if IdvRates(4,nId)./IdvRates(1,nId) <= 1 && countDown<=nbIdv/2
                    countDown = countDown+1 ;
                    IdvNormRates = IdvRates(:,nId)./IdvRates(1,nId) ;
                    patchline(v_Iprtr(id), IdvNormRates(id), 'linestyle','-','edgecolor',cl{i},'edgealpha',.2,'linewidth',1.5) 
                end
                
            end
        end
    end
    
    if(IF_COUNTPIF)
        countPif = 1 ; 
        nbPif = [30 10 10 10] ; 
        nbPif = [30 10 10 10]*2 
        normRatePif = 0 ; 
        while countPif<=nbPif(i) 
            nId = randi([Cpt(i)+1 Cpt(i+1)]) ; 
            if IdvRates(1,nId)>THRESHOLD 
                normRatePif = normRatePif + IdvRates(:,nId)./IdvRates(1,nId) ; 
                countPif = countPif + 1 ; 
            end 
        end
        
        normRatePif = normRatePif ./ nbPif(i) ;
        plot(v_Iprtr, normRatePif, '-o','Color',cl{i})
        
    end

    if(IF_LOGSCALE) 
        xlim([.01 100])
        set(gca,'Xscale', 'log')
        if(prtrPop==1)
            ylim([.1 10])
        else
            ylim([.01 10])
        end
        set(gca,'Yscale', 'log') 
    else
        xlim([.01 100])
        set(gca,'Xscale', 'log')
        if(i==2)
            ylim([0 3])
        elseif(i==4)
            ylim([0 2])
        end
    end

    if(IF_SAVE & (i==2 | i==4))
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
        fprintf('Writing %s \n',figdir)
        try
            mkdir(figdir) ;
        end
        ProcessFigure(fig, fullfile(figdir,figtitle)) ;
    end
    
end

hold off ;
