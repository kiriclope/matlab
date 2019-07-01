clear all ; 
GlobalVars 

Iext = ExternalInput(model,nbpop,dir) ;
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

baseline = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_IEXT, prtrPop, Iext(prtrPop) ) ;

try
    for i=1:length(baseline(1,:))-1
        IdvRates(i) = mean(baseline(:,i+1)) ;
    end
catch
    fprintf('Error Rates not found\n')
    return ;
end

prtr = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_IEXT, prtrPop, Iext(prtrPop) + Iprtr);

try
    for i=1:length(prtr(1,:))-1
        IdvRatesPrtr(i) = mean(prtr(:,i+1)) ;
    end
catch
    fprintf('Error prtr Rates not found\n')
    return ;
end

for i=1:nbpop
            

    MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ; 

    MeanRatePrtr(i) = mean( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRatePrtr(i) = var( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ; 
    
    Baseline = IdvRates( Cpt(i)+1:Cpt(i+1) ) ; 
    RatesPrtr = IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ; 

    [m idx a ROI1] = ratesCutOff(Baseline, Baseline, THRESHOLD, Cth, DIM, L) ; 
    [m idx b ROI2] = ratesCutOff(RatesPrtr, RatesPrtr, THRESHOLD, Cth, DIM, L) ; 
    
    ROI = intersect(ROI1,ROI2) ; 
    %fprintf('length ROI %d ROI1 %d ROI2 %d\n',length(ROI),length(ROI1),length(ROI2)) ;

    Baseline = Baseline(ROI1) ; 
    RatesPrtr = RatesPrtr(ROI2) ; 
    
    ROI1 = find(Baseline<.1) ; 
    Baseline(ROI1)=.1 ; 
    ROI2 = find(RatesPrtr<.1) ; 
    RatesPrtr(ROI2)=.1 ; 

    sc = scatter(Baseline,RatesPrtr,8,cl{i},'MarkerEdgeAlpha',.5, ...
                 'MarkerFaceColor', 'None','LineWidth',.5) ; 

    % ROI1 = find(Baseline>=THRESHOLD) ; 
    % ROI2 = find(RatesPrtr>=THRESHOLD) ; 
    % ROI = intersect(ROI1,ROI2) ; 

    % Baseline = Baseline(ROI) ; 
    % RatesPrtr = RatesPrtr(ROI) ; 

    ROI1 = find(Baseline<.1) ; 
    Baseline(ROI1) = .1 ; 
    ROI2 = find(RatesPrtr<.1) ; 
    RatesPrtr(ROI2) = .1 ; 
    
    SupIdx = ( RatesPrtr-Baseline ) ./ ( Baseline + RatesPrtr )  ; 
    
    SupIdx( ROI1 ) = [] ; 

    propUp = round( length( find(SupIdx>0.05) ) / length(SupIdx) * 100, 0) ;
    propDn = round( length( find(SupIdx<-0.05) ) / length(SupIdx) * 100, 0 ) ; 
    propZo = round( length( find(SupIdx<=-0.95) ) / length(SupIdx) * 100, 0 ) ; 
    propNc = round( ( length(SupIdx) - length( find(SupIdx>=.05) ) - length( find(SupIdx<=-.05) ) ) / length(SupIdx) * 100, 0 ) ; 

    propSum = propUp + propDn + propNc ;
    fprintf('%s propUp %.0f propDn %.0f propZo %.0f propNc %.0f sum %.0f\n', popList(i) , propUp, propDn, propZo, propNc, propSum ) 
 
    % figtitle=sprintf('SupIndex%s',popList(i)) ; 
    % %figtitle=sprintf('SupIndex%s_Iprtr%.3f',popList(i),Iprtr) ; 

    % if( ishandle( findobj('type','figure','name',figtitle) ) )
    %     fig = findobj('type','figure','name',figtitle) ; 
    %     fig = figure(fig); hold on ; 
    % else
    %     fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ;
    %     xlabel('SI')
    %     ylabel('Probability')
    % end
    % histogram(SupIdx,30,'Normalization', 'probability' ,'DisplayStyle','stairs','EdgeColor',cl{i},'LineWidth',1) ;
    % %histogram(-log(SupIdx)/log(10),27,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i},'LineWidth',1) ;
   
    % xlim([-1 1])
    % if(i==1 || i==2) 
    %     ylim([0 .25]) 
    % else 
    %     ylim([0 .4]) 
    % end   
    % drawnow ;

    % if(IF_SAVE)
    %     figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
    %     fprintf('Writing %s \n',figdir)
    %     try
    %         mkdir(figdir)
    %     end
        
    %     ProcessFigure(fig, fullfile(figdir,figtitle)) ;
    %     %ProcessFigure(fig, fullfile(figdir,figtitle),1.5,[1.33*1.5,1.5]) ;
    % end
    % hold off ; 

    figtitle=sprintf('Pie_SupIndex%s_Iprtr%.3f',popList(i),Iprtr) ; 
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ;
    X = [propUp, propDn-propZo, propZo, propNc] ;
    labels = {'Up','Down','Zero','NC'} ;
    labels = {'','','',''} ; 
    %pie(X,labels) 
    p = pie(X) ;
    colormap([1 0 0;     %// red 
              0 0 1;      %// blue
              0 1 0; %// green
             ])

    for j=1:length(p) 
        p(j).EdgeColor = 'None' ;
    end
    axis off ;
    drawnow ;
   
    if(IF_SAVE)
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
        fprintf('Writing %s \n',figdir)
        try
            mkdir(figdir)
        end
        
        ProcessFigure(fig, fullfile(figdir,figtitle), 3, [3 3]) ;
        %ProcessFigure(fig, fullfile(figdir,figtitle),1.5,[1.33*1.5,1.5]) ;
    end
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