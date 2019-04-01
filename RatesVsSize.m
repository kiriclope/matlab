clear all ;
GlobalVars

Iext = ExternalInput(model,nbpop,dir) ;    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

v_Cff = v_Cff(1):v_Cff(2):v_Cff(3) ;

for i=1:length(v_Cff)   
    data = ImportData(model, nbpop, dir, 'IdvRates', N, K, g, IF_RING, Crec, v_Cff(i), IF_DATA, prtrPop, Iext(prtrPop) + Iprtr ) ;
    try
        for j=1:length(data(1,:))-1
            IdvRates(i,j) = mean(data(:,j+1)) ;
        end
    catch
        for j=1:nbN(nbpop)
            IdvRates(i,j) = nan ;
        end
    end

    fprintf(' Rates ') 
    for j=1:nbpop 
        m(j,i) = 0 ;
        Rates = IdvRates(i, Cpt(j)+1:Cpt(j+1)) ;
        [m(j,i) idx] = ratesCutOff(Rates, Rates, THRESHOLD, Cth, DIM, L) ;
        fprintf('%.3f ', m(j,i) )
    end    
    fprintf('\n')
end

for i=1:nbpop

    figtitle = sprintf('%s_RatesVsSize',dir) ;
    fig = popPerFig(i,dir,figtitle) ;
    xlabel('\delta_{opto}')
    ylabel('Norm. Rates')

    NormRates = m(i,:)./m(i,1) ;

    plot(v_Cff, NormRates, 'o','MarkerEdgeColor',cl{i},'MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1)
    plot(v_Cff, NormRates, '-','Color',cl{i})
        
    % for j=1:nbIdv
    %     nId = randi([idx(1)+Cpt(i) idx(end)+Cpt(i)]) ;
    %     if IdvRates(1,nId)>THRESHOLD
    %         IdvNormRates = IdvRates(:,nId)./IdvRates(1,nId) ;
    %         patchline(v_Cff, IdvNormRates, 'linestyle','-','edgecolor',cl{i},'edgealpha',.1) 
    %     end
    % end

    ylim([0 2])

    if( (i==2 || i==4 ) )
        plot(v_Cff,ones(1, length(v_Cff) ) ,'--','Color','k') 
    end

    if(IF_SAVE & (i==2 | i==4)) 
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
        fprintf('Writing %s \n',figdir)                
        try
            mkdir(figdir) ;
        end                
        ProcessFigure(fig, fullfile(figdir,figtitle)) ;
    end
    
    hold off ;
end

