% clear all ;
GlobalVars                              

if(any(IF_NORM*Iprtr)~=0) 
    Prtr = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iprtr,IF_Dij,Dij) ; 
    
    if(IF_PROP)
        BL = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,0,IF_IEXT,prtrPop,Iprtr,IF_Dij,Dij) ; 
    else
        BL = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext,IF_Dij,Dij) ;
    end
else
    Prtr = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iprtr,IF_Dij,Dij) ; 
    BL = Prtr ; 
end

if(IF_PROPWEAK)
    Cff= Cff/sqrt(K) ; 
end

Prtr(:,1) = [] ;
BL(:,1) = [] ;

IdvBL = mean(BL) ; 
IdvPrtr = mean(Prtr) ;

% IdvBL(1:3)
% IdvPrtr(1:3)

% IdvPrtr(find(IdvPrtr<=.01)) = .01 ;
% IdvBL(find(IdvBL<=.01)) = .01 ;

IdvNorm = IdvPrtr./IdvBL ; 

if(length(IdvNorm)==0)
    return 
end

if(~FIGPERPOP)

    if(IF_Dij)
        figtitle = sprintf(['SpatialProfile_' ...
                            'CrecEE%.4fCrecEI%.4f' ...
                            'CrecIE%.4fCrecII%.4f' ...
                            'Cff%.4f_IextI%.4f'], ... 
                           Crec(1)*Dij(1), Crec(2)*Dij(2), ...
                           Crec(1)*Dij(3), Crec(2)*Dij(4), ...
                           Cff,Iprtr(prtrPop)) ; 
    else
        figtitle = sprintf('SpatialProfile_CrecE%.4fCrecI%.4fCff%.4f_IextI%.4f',Crec(1),Crec(2),Cff,Iprtr(prtrPop)) ; 
    end    
    if( ishandle( findobj('type','figure','name',figtitle) ) ) 
        fig = findobj('type','figure','name',figtitle) ; 
        fig = figure(fig); hold on ; 
    else
        fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ;
        xlabel('X (mm)') 
        ylabel('Norm. Rates')  
    end
end

for i=1:nbpop 

    if(FIGPERPOP)
        figtitle = sprintf('%s_SpatialProfile_%s',dir,popList(i)) ; 
        if( ishandle( findobj('type','figure','name',figtitle) ) ) 
            fig = findobj('type','figure','name',figtitle) ; 
            fig = figure(fig); hold on ; 
        else
            fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ;
            xlabel('X (mm)') 
            ylabel('Norm. Rates') 
        end
    end

    Idx = Cpt(i)+1:Cpt(i+1) ;
    Ybl = IdvBL(Idx) ;
    Yprtr = IdvPrtr(Idx) ;

    fprintf('Rates BL ')
    fprintf('%.3f | ', mean(Ybl))
    
    fprintf(' Prtr ')
    fprintf('%.3f | ', mean(Yprtr) )

    if(IF_NORM)
        fprintf(' Norm ') 
        fprintf('%.3f | ', mean(Yprtr)./mean(Ybl) ) 
    end

    if(IF_PROP)

        if(IF_PROPWEAK)
            [mUp mDn] = ratesProp(Yprtr, Ybl, THRESHOLD, Cff*L/2/sqrt(K), DIM, L) ; 
        else 
            [mUp mDn] = ratesProp(Yprtr, Ybl, THRESHOLD, Cff*L/2, DIM, L) ; 
        end
        
        fprintf(' Up ') 
        fprintf('%.3f | ', mUp./mean(Ybl)) 
        fprintf(' Dn ') 
        fprintf('%.3f | ', mDn./mean(Ybl)) 
    end

    fprintf('\n') ; 

    X =  linspace(-L/2,L/2,nbN(i)) ; 

    fitBL = smooth(X,Ybl,.1,'lowess') ; 
    fitPrtr = smooth(X,Yprtr,.1,'lowess') ; 
    
    if(IF_NORM)
        plt = patchline(X,fitPrtr./fitBL,'linestyle','-', ...
                        'edgecolor',cl{i},'edgealpha',alp,'linewidth',2) ; 
    else
        plt = patchline(X,fitPrtr,'linestyle','-', ...
                        'edgecolor',cl{i},'edgealpha',alp,'linewidth',2) ; 
    end

    plot(X,ones(1, length(X) ),'--','Color','k') 

    if(~IF_PROP)
        plot(4*Crec(i) *[1 1],[0 2], '--','Color',cl{i}) 
        xlim([0 1]) 
    end
end

plot(4 * Cff *[1 1],[0 2], '--','Color','k') 


if(IF_SAVE) 
    ylim([0 2])
    
    if(IF_Dij)
        figdir = sprintf( ['./Space/%dpop/%s/N%dK%dg%.2f/%s/%sPrtr/' ...
                           'CrecEE%.4fCrecEI%.4f' ...
                           'CrecIE%.4fCrecII%.4f'], ...
                          nbpop,dir,N,K,1,IF_RING, IF_IEXT, ...
                          Crec(1)*Dij(1), Crec(2)*Dij(2), ...
                          Crec(1)*Dij(3), Crec(2)*Dij(4) ) ; 
    else
        figdir = sprintf(['./Space/%dpop/%s/N%dK%dg%.2f/%s/%sPrtr'],nbpop,dir,N,K,1,IF_RING,IF_IEXT) ; 
    end
    
    fprintf('Writing %s \n',figdir) 
    try
        mkdir(figdir) ;
    end
    ProcessFigure(fig, fullfile(figdir,figtitle)) ;
end
