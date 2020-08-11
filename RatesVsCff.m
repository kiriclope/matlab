%clear all ; 
GlobalVars

for i=1:length(v_Cff)   
    data = ImportData(model, nbpop, dir, 'IdvRates', N, K, g, IF_RING, ...
                      Crec, v_Cff(i), IF_IEXT, prtrPop, Iext(prtrPop) + Iprtr) ;
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

            % m(j,i) = ratesCutOff(Rates, BL, THRESHOLD, Cth, DIM,
            % L) ; 
            
            m(j,i) = mean(Rates) ;
            fprintf('%.3f ', round(m(j,i),3) ) 

            if(j==prtrPop)
                if(IF_PROPWEAK)
                    [mUp(i) mDn(i)] = ratesProp(Rates, BL, THRESHOLD, v_Cff(i)*L/2/sqrt(K), DIM, L) ; 
                else  
                    [mUp(i) mDn(i)] = ratesProp(Rates, BL, THRESHOLD, v_Cff(i)*L/2, DIM, L) ; 
                end 
                fprintf('%.3f %.3f ', round(mUp(i),3), round(mDn(i),3)) 
            end
            
        end 
        fprintf('\n')

        Gext = Iext ;
        Gext(prtrPop) = Iext(prtrPop) + v_Cff(i) ;
        fprintf('MF Rates ') 
        RatesMF = BalRatesMF(model,nbpop,dir,Gext*.01,J,0) ; 
        for j=1:nbpop 
            MFrates(j,i) = RatesMF(j)*1000 ; 
            fprintf('%.2f ', round(MFrates(j,i),2) ) 
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
    figtitle = sprintf('%s_RatesVsCff',dir) ;

    if( ishandle( findobj('type','figure','name',figtitle) ) )
        fig = findobj('type','figure','name',figtitle) ; 
        fig = figure(fig); hold on ; 
    else
        fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
            if(IF_PROP)
                if(IF_PROPWEAK)
                    xlabel('p / $\sqrt{K}$') 
                else
                    xlabel('p') 
                end
            else
                xlabel('Cff') 
            end

        if(IF_NORM)
            ylabel('Norm. Rates') 
        else
            ylabel('Rates') 
        end
    end
end

for i=1:nbpop
    
    if(FIGPERPOP)
        if(i==1 || i==2) 
            figtitle = sprintf('%s_RatesVsCff%s_EI',dir,popList(prtrPop)) ; 
        else 
            figtitle = sprintf('%s_RatesVsCff%s_SV',dir,popList(prtrPop)) ; 
        end

        if( ishandle( findobj('type','figure','name',figtitle) ) )
            fig = findobj('type','figure','name',figtitle) ; 
            fig = figure(fig); hold on ; 
        else

            fig = popPerFig(i,dir,figtitle) ;
            if(IF_PROP)
                if(IF_PROPWEAK)
                    xlabel('p / $\sqrt{K}$') 
                else
                    xlabel('p') 
                end
            else
                xlabel('Cff') 
            end
            
            if(IF_NORM)
                ylabel('Norm. Rates') 
            else 
                ylabel('Rates') 
            end 
        end
    end

    if(IF_NORM)
        NormRates = m(i,:)./m(i,1) ;
    else
        NormRates = m(i,:) ;
    end

    id = find( ~isnan(NormRates) ) ;
    id = 1:length(v_Cff) ; 

    
    patchline(abs(v_Cff(id)), NormRates(id), 'linestyle','-','edgecolor',cl{i},'edgealpha',alp,'linewidth',1.5) 
    plot(abs(v_Cff(id)), NormRates(id), mk,'MarkerEdgeColor', ...
         cl{i},'MarkerSize',mkSize,'MarkerFaceColor','none','LineWidth', alp) 
    
    if(i==prtrPop)
        NormUp = mUp(:)./mDn(1) ; 
        NormDn = mDn(:)./mDn(1) ; 
        
        patchline(abs(v_Cff(id)), NormUp(id), 'linestyle','-','edgecolor','g','edgealpha',alp,'linewidth',1.5) 
        plot(abs(v_Cff(id)), NormUp(id), '^','MarkerEdgeColor', ...
             'g','MarkerSize',mkSize,'MarkerFaceColor','none','LineWidth', alp) 
        
        patchline(abs(v_Cff(id)), NormDn(id), 'linestyle','-','edgecolor','m','edgealpha',alp,'linewidth',1.5) 
        plot(abs(v_Cff(id)), NormDn(id), 'v','MarkerEdgeColor', ...
             'm','MarkerSize',mkSize,'MarkerFaceColor','none','LineWidth', alp) 
    end            
    
    if( (i==2 || i==4 ) && IF_NORM)
        if(IF_MF_RATES)
            if(i==2)
                for j=1:2
                    MFpop = MFrates(j,:)./MFrates(j,1) ;
                    plot(abs(v_Cff(id)), MFpop(id), '--','Color',cl{j}) 
                end
            elseif(i==4)
                for j=3:4
                    MFpop = MFrates(j,:)./MFrates(j,1) ;
                    plot(abs(v_Cff(id)), MFpop(id), '--','Color',cl{j}) 
                end
            end

        end
        
        plot(abs(v_Cff(1:end)),ones(1, length(v_Cff(1:end)) ),'--','Color','k') 
        
    end
    
    if(IF_IDVTraces) 
        countUp = 1 ;
        countDown = 1 ;
        countMax = 1 ;
        while countUp+countDown<nbIdv && countMax<100 
            nId = randi([Cpt(i)+1 Cpt(i+1)]) ; 
            countMax = countMax + 1 ;
            
            if IdvRates(1,nId)>=0
                
                if IdvRates(4,nId)./IdvRates(1,nId) >= 1 && countUp<=nbIdv/2
                    countUp = countUp+1 ;
                    IdvNormRates = IdvRates(:,nId)./IdvRates(1,nId) ;
                    patchline(v_Cff(id), IdvNormRates(id), 'linestyle','-','edgecolor',cl{i},'edgealpha',.2,'linewidth',1.5) 
                    
                end                    
                if IdvRates(4,nId)./IdvRates(1,nId) <= 1 && countDown<=nbIdv/2
                    countDown = countDown+1 ;
                    IdvNormRates = IdvRates(:,nId)./IdvRates(1,nId) ;
                    patchline(v_Cff(id), IdvNormRates(id), 'linestyle','-','edgecolor',cl{i},'edgealpha',.2,'linewidth',1.5) 
                end
                
            end
        end
    end
    
    if(IF_SAVE & (i==2 || i==4)) 
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
        fprintf('Writing %s \n',figdir)
        try
            mkdir(figdir) ;
        end
        ProcessFigure(fig, fullfile(figdir,figtitle)) ;
    end
    
end

hold off ;
