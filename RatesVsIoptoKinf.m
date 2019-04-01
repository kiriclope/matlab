clear all ; 
GlobalVars 

IextBL = ExternalInput(model,nbpop,dir) ; 
v_Iprtr = v_Iprtr(1):v_Iprtr(2):v_Iprtr(3) ; 
FIGPERPOP = 0 ; 
IF_LOGSCALE = 0 ; 
IF_POWER = 0 ; 
J = ImportJab(model,nbpop,dir) ; 
Iext = IextBL ; 

for i=1:length(v_Iprtr) 
    Iext(prtrPop) = IextBL(prtrPop) + v_Iprtr(i) ; 
    [RatesMF DetJ]= BalRatesMF(model,nbpop,dir,Iext,J,0) ; 

    fprintf('I_opto ') 
    fprintf('%.3f ', v_Iprtr(i)) 

    fprintf('Det ') 
    fprintf('%.3f ', DetJ) 

    fprintf(' Rates ') 
    for j=1:nbpop 
        m(j,i) = RatesMF(j) ; 
        fprintf('%.3f ', m(j,i) ) 
    end 
    fprintf('\n') 
end 

if(~FIGPERPOP) 
    figtitle = sprintf('%s_RatesVsIopto%s_MF',dir,popList(prtrPop)) ; 

    if( ishandle( findobj('type','figure','name',figtitle) ) ) 
        fig = findobj('type','figure','name',figtitle) ; 
        fig = figure(fig); hold on ; 
    else
        if(IF_NORM) 
            figtitle = sprintf('%s_Norm',figtitle) ; 
        end
        fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
        xlabel('I_{opto} (\mu A . cm^{-2})') 
        if(IF_NORM) 
            ylabel('Norm. Activity') 
        else
            ylabel('Activity (Hz)') 
        end
    end
end

for i=nbpop:-1:1    

    if(FIGPERPOP)
        if(i==1 || i==2) 
            figtitle = sprintf('%s_RatesVsIopto%s_MF_EI',dir,popList(prtrPop)) ; 
        else 
            figtitle = sprintf('%s_RatesVsIopto%s_MF_SV',dir,popList(prtrPop)) ; 
        end
        fig = popPerFig(i,dir,figtitle) ;
        xlabel('I_{opto}') 
        ylabel('Norm. Activity') 
    end
    
    if(IF_NORM)
        NormRates = m(i,:) ./ m(i,1) ; 
    else 
        NormRates = m(i,:) ; 
        ylabel('Activity (Hz)') 
    end

    if( (i==2 || i==4 ) && IF_NORM) 
        plot(v_Iprtr(IDX:end),ones(1, length(v_Iprtr(IDX:end)) ),'--','Color','k')
    end

    if(i==1 && IF_NORM)
        if(FIGPERPOP || nbpop==2) 
            plot(v_Iprtr(IDX:end)  , NormRates(IDX:end), '-','Color',cl{i})
        else
            plot(v_Iprtr(IDX:40:end), NormRates(IDX:40:end), '+', ...
                 'MarkerEdgeColor',cl{i},'MarkerSize',6,'MarkerFaceColor', ...
                 cl{i},'LineWidth', 1) 
        end
    else
        plot(v_Iprtr(IDX:end)  , NormRates(IDX:end), '-','Color',cl{i}) 
    end
    
    drawnow ; 
    
    if(IF_LOGSCALE) 
        xlim([.01 100]) 
        set(gca,'Xscale', 'log') 
        ylim([.01 10]) 
        set(gca,'Yscale', 'log') 
    else
        if(IF_LOGSCALEX) 
            xlim([.01 100]) 
            set(gca,'Xscale', 'log')
            ylim([0 2])
        end
    end

    if(IF_SAVE & i==1)
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ; 
        fprintf('Writing %s \n',figdir) 
        try
            mkdir(figdir) ; 
        end
        ProcessFigure(fig, fullfile(figdir,figtitle)) ; 
    end

end

