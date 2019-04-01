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
    catch
        for j=1:nbN(nbpop)
            IdvRates(i,j) = nan ;
        end
    end
end

if(IF_POWER~=0) 
    if IF_POWER==1
        v_Iprtr = P0 .* ( exp( v_Iprtr ./ I0 ) - 1 )  ; % profile Moi
    else
        v_Iprtr = P0 .* ( exp( v_Iprtr *20 ./ I0 ) - 1 ) ; 
    end
end

if(~FIGPERPOP)
    figtitle = sprintf('%s_SuppIdxVsIopto',dir) ; 
    if( ishandle( findobj('type','figure','name',figtitle) ) )
        fig = findobj('type','figure','name',figtitle) ; 
        fig = figure(fig); hold on ; 
    else
        fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
        xlabel('I_{opto}') 
        ylabel('SI') 
    end
end

for i=1:nbpop

    if(FIGPERPOP)
        if(i==1 || i==2) 
            figtitle = sprintf('%s_SuppIdxVsIopto_EI',dir) ; 
        else 
            figtitle = sprintf('%s_SuppIdxVsIopto_SV',dir) ; 
        end
        fig = popPerFig(i,dir,figtitle) ;
        xlabel('I_{opto}') 
        ylabel('SI') 
    end

    if(IF_POWER~=0) 
        figtitle = sprintf('%s_%d',figtitle,IF_POWER) ; 
        xlabel('\Gamma_{opto} (mW/mm^2)') 
    end
    
    countUp = 1 ;
    countDown = 1 ;
    while countUp+countDown<nbIdv 
        nId = randi([Cpt(i) Cpt(i+1)]) ; 
        if IdvRates(1,nId)>THRESHOLD 
            if IdvRates(4,nId)./IdvRates(1,nId) > 1 && countUp<=nbIdv/2
                countUp = countUp+1 ;
                SuppIdx = - ( ( IdvRates(1,nId) - IdvRates(:,nId) ) ./ ( IdvRates(1,nId) + IdvRates(:,nId) ) ) ; 
                patchline(v_Iprtr, SuppIdx, 'linestyle','-','edgecolor',cl{i},'edgealpha',.25,'linewidth',1.5) 
            end  
            if IdvRates(4,nId)./IdvRates(1,nId) < 1 && countDown<=nbIdv/2
                countDown = countDown+1 ; 
                %SuppIdx = ( IdvRates(1,nId) - IdvRates(:,nId) ) ./ ( IdvRates(1,nId) + IdvRates(:,nId) ) ;
                SuppIdx = - ( ( IdvRates(1,nId) - IdvRates(:,nId) ) ./ ( IdvRates(1,nId) + IdvRates(:,nId) ) ) ; 
                patchline(v_Iprtr, SuppIdx, 'linestyle','-','edgecolor',cl{i},'edgealpha',.25,'linewidth',1.5) 
            end            
        end
    end

    ylim([-1 1]) 
    xlim([0 4]) 

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
