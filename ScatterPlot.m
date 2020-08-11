clear all ;
GlobalVars

baseline = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_IEXT, prtrPop, Iext ) ;

try
    for i=1:length(baseline(1,:))-1
        IdvRates(i) = mean(baseline(:,i+1)) ;
    end
catch
    fprintf('Error Rates not found\n')
    return ;
end

prtr = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_IEXT, prtrPop, Iprtr) ; 

try
    for i=1:length(prtr(1,:))-1
        IdvRatesPrtr(i) = mean(prtr(:,i+1)) ;
    end
catch
    fprintf('Error prtr Rates not found\n')
    return ;
end

for i=1:2
            
    figname=sprintf('Scatter%s_Iprtr%.3f',popList(i),Iprtr) ;
    fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
    
    MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ; 

    MeanRatePrtr(i) = mean( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRatePrtr(i) = var( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ; 
    
    Baseline = IdvRates( Cpt(i)+1:Cpt(i+1) ) ;
    RatesPrtr = IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ; 

    % [m idx Rates ROI] = ratesCutOff(Baseline, Baseline, THRESHOLD, Cth, DIM, L) ;    
    % RatesPrtr = RatesPrtr(ROI) ;

    [m idx a ROI1] = ratesCutOff(Baseline, Baseline, THRESHOLD, Cth, DIM, L) ; 
    [m idx b ROI2] = ratesCutOff(RatesPrtr, RatesPrtr, THRESHOLD, Cth, DIM, L) ; 
        
    % ROI = intersect(ROI1,ROI2) ; 
    Baseline = Baseline(ROI1) ; 
    RatesPrtr = RatesPrtr(ROI2) ; 
    
    ROI1 = find(Baseline<.1) ; 
    Baseline(ROI1)=.1 ; 
    ROI2 = find(RatesPrtr<.1) ; 
    RatesPrtr(ROI2)=.1 ; 

    whos Baseline
    whos RatesPrtr
    
    plot([.1 100],[.1 100],'--k','Linewidth',.5) 

    Idx = randi([1,length(Baseline)],300,1) ;
    sc = scatter(Baseline(Idx),RatesPrtr(Idx),8,cl{i},'MarkerEdgeAlpha',.5, ...
                 'MarkerFaceColor', 'None','LineWidth',.5) ; 
    % sc = scatterhist(Baseline(Idx),RatesPrtr(Idx),'MarkerSize',1,'Color',cl{i}) ; 

     % [coeff,score,latent] = pca([Rates.',RatesPrtr.']) ;
     % sc = scatter(score(:,1),score(:,2),.1,cl{i},'MarkerEdgeAlpha',1) ;

     % set(sc,'SizeData',2) ;
    
    % xLimits = get(gca,'XLim');
    % yLimits = get(gca,'YLim');
    plot([.1 100],[.1 100],'--k','Linewidth',.5) 
    xlim([.1 100])
    ylim([.1 100])
    set(gca,'yscale','log')
    set(gca,'xscale','log')

    % set(sc(2),'xdir','reverse')
    % set(sc(3),'xdir','reverse')

    % set(sc(2),'ydir','reverse')
    % set(sc(3),'ydir','normal')

    % set(sc(2),'xscale','log')
    % set(sc(3),'xscale','log')

         
    %Aylim(yLimits)
    
    xlabel('Baseline (Hz)')
    ylabel('Perturbartion (Hz)')
    drawnow ;

    if(IF_SAVE)
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
        fprintf('Writing %s \n',figdir)
        try
            mkdir(figdir)
        end
        
        %ProcessFigure(fig, fullfile(figdir,figname),1.25,[1.33*1.25,1.25]) ;
        ProcessFigure(fig, fullfile(figdir,figname)) ;
    end
    pause(.2) ;
    hold off ;
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