% clear all ; 
GlobalVars 

Iext = ExternalInput(model,nbpop,dir) ;
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

baseline = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_IEXT, prtrPop, Iext ) ;

try
    for i=1:length(baseline(1,:))-1
        IdvRates(i) = mean(baseline(:,i+1)) ;
    end
catch
    fprintf('Error Rates not found\n')
    return ;
end

prtr = ImportData(model, nbpop, dir,'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_IEXT, prtrPop, Iprtr);

try
    for i=1:length(prtr(1,:))-1
        IdvRatesPrtr(i) = mean(prtr(:,i+1)) ;
    end
catch
    fprintf('Error prtr Rates not found\n')
    return ;
end

for i=1:2
                
    MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ; 

    MeanRatePrtr(i) = mean( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRatePrtr(i) = var( IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ) ; 
    
    Baseline = IdvRates( Cpt(i)+1: Cpt(i+1) ); 
    RatesPrtr = IdvRatesPrtr( Cpt(i)+1: Cpt(i+1) ) ; 

    % [m idx a ROI1] = ratesCutOff(Baseline, Baseline, THRESHOLD, Cth, DIM, L) ; 
    % [m idx b ROI2] = ratesCutOff(RatesPrtr, RatesPrtr, THRESHOLD, Cth, DIM, L) ; 
    
    % ROI = intersect(ROI1,ROI2) ; 
    % %fprintf('length ROI %d ROI1 %d ROI2 %d\n',length(ROI),length(ROI1),length(ROI2)) ;

    % Baseline = Baseline(ROI1) ; 
    % RatesPrtr = RatesPrtr(ROI2) ;  

    ROI1 = find(Baseline<=.01) ; 
    Baseline(ROI1)=.01 ; 
    ROI2 = find(RatesPrtr<=.01) ; 
    RatesPrtr(ROI2)=.01 ; 

    ROI = intersect(ROI1,ROI2) ; 
    ZeroIdx = union(ROI,ROI2) ; 

    % SupIdx( ROI1 ) = [] ; 
    SupIdx = ( RatesPrtr-Baseline ) ; %./ ( Baseline + RatesPrtr ) ; 
    Ratio = SupIdx ./ Baseline ;

    % if( any( union(find(SupIdx<-0.95),ZeroIdx) ) ) 
    %     IdxZo = union(find(SupIdx<-0.95),ZeroIdx) ; 
    % else
    %     IdxZo = [] ;
    % end
    
    IdxZo = [] ;

    if(~isempty(ZeroIdx))
        IdxZo = ZeroIdx ;
    else
        IdxZo = [] ;
    end

    if(any(abs(Ratio)<=0.01))
        IdxNc = setdiff(find(abs(Ratio)<=0.01),IdxZo) ; 
    else
        IdxNc = [] ; 
    end
    
    if(any(SupIdx>0)) 
        IdxUp = setdiff(setdiff(find(SupIdx>0),IdxZo),IdxNc) ; 
    else
        IdxUp = [] ; 
    end
    if(any(SupIdx<0))
        IdxDn = setdiff(setdiff(find(SupIdx<0),IdxZo),IdxNc) ; 
    else
        IdxDn = [] ;
    end

    propSum=0 ;
    if(~isempty(SupIdx)) 
        propUp = round( length(IdxUp) / length(SupIdx) * 100 ) ;
        propDn = round( length(IdxDn) / length(SupIdx) * 100 ) ; 
        propZo = round( length(IdxZo) / length(SupIdx) * 100 ) ;
        propNc = round( length(IdxNc) / length(SupIdx) * 100 ) ; 
        propSum = propUp + propDn + propNc + propZo ;        
    else
        propUp = 0 ;
        propZo = 100 ;
        propDn = 0 ;
        propNc = 0 ;
        propSum = propUp + propDn + propNc + propZo ;
    end

    fprintf('%s propUp %.0f propDn %.0f propZo %.0f propNc %.0f sum %.0f\n', popList(i) , propUp, propDn, propZo, propNc, propSum ) 

    while(propSum~=100)
        if(propSum>100)
            if(propZo>=1)
                propZo = propZo - 1 ;
            elseif(propZo<1)
                if(propNc>=1)
                    propNc = propNc -1 ;
                end
            end
        else 
            if(propZo>=1)
                propZo = propZo + 1 ;
            elseif(propZo<1)
                if(propNc>=1)
                    propNc = propNc + 1 ;
                end
            end
        end

        propSum = propUp + propDn + propNc + propZo ;
    end

    fprintf('%s propUp %.0f propDn %.0f propZo %.0f propNc %.0f sum %.0f\n', popList(i) , propUp, propDn, propZo, propNc, propSum ) 


    Idx = randi([Cpt(i)+1 Cpt(i+1)],300,1) ; 
    Baseline = IdvRates(Idx) ; 
    RatesPrtr = IdvRatesPrtr(Idx) ; 

    % [m idx a ROI1] = ratesCutOff(Baseline, Baseline, THRESHOLD, Cth, DIM, L) ; 
    % [m idx b ROI2] = ratesCutOff(RatesPrtr, RatesPrtr, THRESHOLD, Cth, DIM, L) ; 
    
    % ROI = intersect(ROI1,ROI2) ; 
    % %fprintf('length ROI %d ROI1 %d ROI2 %d\n',length(ROI),length(ROI1),length(ROI2)) ;

    % Baseline = Baseline(ROI1) ; 
    % RatesPrtr = RatesPrtr(ROI2) ;  

    ROI1 = find(Baseline<=.01) ; 
    Baseline(ROI1)=.01 ; 
    ROI2 = find(RatesPrtr<=.01) ; 
    RatesPrtr(ROI2)=.01 ; 

    ROI = intersect(ROI1,ROI2) ; 
    ZeroIdx = union(ROI,ROI2) ; 

    % SupIdx( ROI1 ) = [] ; 
    SupIdx = ( RatesPrtr-Baseline ) ; %./ ( Baseline + RatesPrtr ) ; 
    Ratio = SupIdx ./ Baseline ;

    % if( any( union(find(SupIdx<-0.95),ZeroIdx) ) ) 
    %     IdxZo = union(find(SupIdx<-0.95),ZeroIdx) ; 
    % else
    %     IdxZo = [] ;
    % end
    
    IdxZo = [] ;

    if(~isempty(ZeroIdx))
        IdxZo = ZeroIdx ;
    else
        IdxZo = [] ;
    end

    if(any(abs(Ratio)<=0.01))
        IdxNc = setdiff(find(abs(Ratio)<=0.01),IdxZo) ; 
    else
        IdxNc = [] ; 
    end
    
    if(any(SupIdx>0)) 
        IdxUp = setdiff(setdiff(find(SupIdx>0),IdxZo),IdxNc) ; 
    else
        IdxUp = [] ; 
    end
    if(any(SupIdx<0))
        IdxDn = setdiff(setdiff(find(SupIdx<0),IdxZo),IdxNc) ; 
    else
        IdxDn = [] ;
    end
 
    figtitle=sprintf('Scatter%s_Iprtr%.2f',popList(i),Iprtr(prtrPop)) ; 

    % if( ishandle( findobj('type','figure','name',figtitle) ) )
    %     fig = findobj('type','figure','name',figtitle) ; 
    %     fig = figure(fig); hold on ; 
    % else
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ;
    xlabel('Baseline (Hz)')
    ylabel('Light On (Hz)')
    % end

    scatter(Baseline,RatesPrtr,20,cl{i},'filled','MarkerFaceAlpha',.5,'LineWidth',.5) ; 
    % scatter(Baseline(IdxUp),RatesPrtr(IdxUp),10,'r','filled','MarkerFaceAlpha',.5,'LineWidth',.5) ; 
    % scatter(Baseline(IdxDn),RatesPrtr(IdxDn),10,'b','filled','MarkerFaceAlpha',.5,'LineWidth',.5) ; 
    % scatter(Baseline(IdxNc),RatesPrtr(IdxNc),10,'g','filled','MarkerFaceAlpha',.5,'LineWidth',.5) ;
    % scatter(Baseline(IdxZo),RatesPrtr(IdxZo),10,'w','filled','MarkerFaceAlpha',.5,'MarkerEdgeColor','b','LineWidth',.1) ; 
    
    plot([.01 100],[.01 100],'--k','Linewidth',.5) 
    xlim([.01 100])
    ylim([.01 100])
    set(gca,'yscale','log')
    set(gca,'xscale','log')

    drawnow ;

    if(IF_SAVE)
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
        figdir = sprintf('%s/PieChart',figdir) ;

        fprintf('Writing %s \n',figdir)
        try
            mkdir(figdir)
        end
        
        ProcessFigure(fig, fullfile(figdir,figtitle)) ;
        %ProcessFigure(fig, fullfile(figdir,figtitle),1.5,[1.33*1.5,1.5]) ;
    end
    hold off ; 

    figtitle=sprintf('PieProp%s_Iprtr%.3f',popList(i),Iprtr(prtrPop)) ; 
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ;
    X = [propUp, propNc, propDn, propZo] ;
    labels = {'Up','NC','Down','Zero'} ;
    labels = {'','','',''} ; 
    %pie(X,labels) 
    p = pie(X) ; 
    pieColor = gray(4) ; 
    isProp = logical([propUp propNc propDn propZo]) ;
    for l=4:-1:1
        if(~isProp(l))
            pieColor(l,:) = [] ;
        end
    end
    colormap(pieColor) 

    % colormap([0 1 1; %// red 
    %           1 0 1; %// green
    %           0 1 0; %// blue
    %           1 1 1
    %          ]) 

    % if(propUp==0 && propNc~=0 && propDn~=0 && propZo~=0) 
    %     colormap([1 0 1; % // green
    %               0 1 0; % // blue
    %               1 1 1 
    %              ])
    % elseif(propUp==0 && propNc==0 && propDn~=0 && propZo~=0)
    %     colormap([1 0 1; % // green
    %               1 1 1 
    %              ])
    % elseif(propUp==0 && propNc==0 && propDn==0 && propZo~=0)
    %     colormap([1 1 1])
    % elseif(propUp~=0 && propNc==0 && propDn~=0 && propZo~=0)
    %     colormap([0 1 1; %// red 
    %               0 1 0; %// blue
    %               1 1 1
    %              ])    
    % elseif(propUp~=0 && propNc~=0 && propDn~=0 && propZo==0)
    %     colormap([0 1 1; %// red 
    %               1 0 1; %// green
    %               0 1 0; %// blue
    %              ])    
    % end
    
    % for j=1:length(p) 
    %     p(j).EdgeColor = 'None' ;
    % end
    axis off ; 
    drawnow ; 
   
    if(IF_SAVE)
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;

        figpath = sprintf('%s/PieChart',figdir) ; 
        fprintf('Writing %s \n',figpath)

        try
            mkdir(figdir)
        end
        
        ProcessFigure(fig, fullfile(figpath,figtitle), 3, [3 3]) ;
        %ProcessFigure(fig, fullfile(figdir,figtitle),1.5,[1.33*1.5,1.5]) ;
    end
    hold off ;
end


MF = BalRatesMF(model,nbpop,dir,Iext,[]) ;
fprintf('MF : ')
fprintf('%.3f | ', MF)
fprintf('\n')

fprintf('Norm. Rates: ')
fprintf('%.3f | ', MeanRatePrtr./MeanRate)
fprintf('\n')
