clear all ;
GlobalVars

Iext = ExternalInput(model,nbpop,dir) ;    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

baseline = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, Iext(prtrPop) ) ;

try
    for i=1:length(baseline(1,:))-1
        IdvRates(i) = mean(baseline(:,i+1)) ;
    end
catch
    fprintf('Error Rates not found\n')
    return ;
end

prtr = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, Iext(prtrPop) + Iprtr) ; 

try
    for i=1:length(prtr(1,:))-1
        IdvRatesPrtr(i) = mean(prtr(:,i+1)) ;
    end
catch
    fprintf('Error prtr Rates not found\n')
    return ;
end

for i=1:nbpop
            
    figname=sprintf('SupIndex%s_Iprtr%.3f',popList(i),Iprtr) ; 
    fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
    
    MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ; 

    MeanRatePrtr(i) = mean( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRatePrtr(i) = var( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ; 
    
    Baseline = IdvRates( Cpt(i)+1:Cpt(i+1) ) ; 
    RatesPrtr = IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ; 

    [m idx a ROI1] = ratesCutOff(Baseline, Baseline, THRESHOLD, Cth, DIM, L) ; 
    [m idx b ROI2] = ratesCutOff(RatesPrtr, RatesPrtr, THRESHOLD, Cth, DIM, L) ; 
        
    ROI = intersect(ROI1,ROI2) ; 
    fprintf('length ROI %d ROI1 %d ROI2 %d\n',length(ROI),length(ROI1),length(ROI2)) ; 

    Baseline = Baseline(ROI1) ; 
    RatesPrtr = RatesPrtr(ROI2) ; 

    ROI1 = find(Baseline>THRESHOLD) ; 
    ROI2 = find(RatesPrtr>THRESHOLD) ; 
    ROI = intersect(ROI1,ROI2) ; 
    fprintf('length ROI %d ROI1 %d ROI2 %d\n',length(ROI),length(ROI1),length(ROI2)) ;

    Baseline = Baseline(ROI) ; 
    RatesPrtr = RatesPrtr(ROI) ; 
    
    SupIdx = ( RatesPrtr-Baseline ) ./ ( Baseline + RatesPrtr )  ;     
    histogram(SupIdx,30,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i},'LineWidth',1) ;
    %histogram(-log(SupIdx)/log(10),27,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i},'LineWidth',1) ;

    xlim([-1 1]) 
    ylim([0 5])
    xticks([0 1 2 3 4 5])
    xlabel('SI')
    ylabel('pdf')
    drawnow ;

    if(IF_SAVE)
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
        fprintf('Writing %s \n',figdir)
        try
            mkdir(figdir)
        end
        
        ProcessFigure(fig, fullfile(figdir,figname),1.5,[1.33*1.5,1.5]) ;
    end
    hold off ;
    pause(.2)

end


MF = BalRatesMF(model,nbpop,dir,Iext,[]) ;
fprintf('MF : ')
fprintf('%.3f | ', MF)
fprintf('\n')

fprintf('Simuls : ')
fprintf('%.3f | ', MeanRate)
fprintf('\n')

fprintf('Simuls : ')
fprintf('%.3f | ', MeanRatePrtr)
fprintf('\n')