clear all ;
GlobalVars

IextBL = ExternalInput(model,nbpop,dir) ; 
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

Iprtr = .05 ;
Iext = IextBL ; 
Iext(prtrPop) = IextBL(prtrPop) + Iprtr ; 


for i=1:100
    dirIdx = sprintf('%s_RND_%d',dir,i) ;
    
    try 
        BL = ImportData(model, nbpop, dirIdx, 'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, IextBL(prtrPop) ) ;  
        for j=1:length(BL(1,:))-1
            IdvRatesBL(i,j) = mean( BL(:,j+1) ) ;
        end
        
        for j=1:nbpop 
            mBL(j,i) = nan ; 
            RatesBL = IdvRatesBL(i, Cpt(j)+1:Cpt(j+1)) ; 
            [mBL(j,i) idx] = ratesCutOff(RatesBL, RatesBL, THRESHOLD, Cth, DIM, L) ;
        end                                 
        
        % fprintf(' RatesBL ') 
        % for j=1:nbpop         
        %     fprintf('%.3f ', mBL(j,i) )
        % end    
        % fprintf('\n')
        
        Prtr = ImportData(model, nbpop, dirIdx, 'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, Iext(prtrPop) ) ; 
        for j=1:length(Prtr(1,:))-1
            IdvRatesPrtr(i,j) = mean(Prtr(:,j+1)) ;
        end
        
        for j=1:nbpop 
            mPrtr(j,i) = nan ;
            RatesPrtr = IdvRatesPrtr(i, Cpt(j)+1:Cpt(j+1)) ; 
            [mPrtr(j,i) idx] = ratesCutOff(RatesPrtr, RatesPrtr, THRESHOLD, Cth, DIM, L) ;
        end 
        
        % fprintf(' RatesPrtr ') 
        % for j=1:nbpop 
        %     fprintf('%.3f ', mPrtr(j,i) ) 
        % end 
        % fprintf('\n') 
        
        if( all(mBL(:,i)>.01) && all(mPrtr(:,i)>.01) ) 
            for j=1:nbpop 
                khi(j,i) = ( mPrtr(j,i) - mBL(j,i) ) ./ Iprtr ; 
            end 
            
            % if strfind(dir,'L5') 
            %     if(khi(2,i)>0) 
            %         for j=1:nbpop 
            %             mBL(j,i) = nan ;
            %             khi(j,i) = nan ;
            %         end
            %     end
            % end

            if strfind(dir,'L23') 
                if(i==11 || i==39 || i==48 || i==51 || i==61 || i==77 ...
                   || i==85 || i==88 || i==96 ) 
                    for j=1:nbpop 
                        mBL(j,i) = nan ;
                        khi(j,i) = nan ;
                    end
                end
            end
            
            fprintf(' Khi ') 
            fprintf('%.3f ', khi(:,i) ) 
            fprintf('\n') 
        else
            for j=1:nbpop 
                mBL(j,i) = nan ;
                khi(j,i) = nan ;
            end            
        end 
        
    catch         
        for j=1:nbpop 
            khi(j,i) = nan ; 
        end 
    end 
    
end 

for i=1:nbpop 
    
    switch i 
      case 1
        figname='RatesDistE' ;
      case 2
        figname='RatesDistI' ;
      case 3
        figname='RatesDistS' ;
      case 4
        figname='RatesDistV' ;
    end
    
    if( ishandle( findobj('type','figure','name',figname) ) ) 
        fig = findobj('type','figure','name',figname) ; 
        figure(fig); hold on ; 
    else 
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        xlabel('Activity (Hz)') 
        ylabel('pdf') 
    end

    if strfind(dir,'L5LAST') 
        xlim([0 30]) 
    else
        xlim([0 15]) 
    end
    
    idx = ~isnan( mBL(i,:) ) ; 
    X = mBL(i,idx) ; 
    X=X(X~=0) ; 
    
    h = histogram( X , 10, 'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i},'EdgeAlpha',alp,'Linewidth',2) ; 
    drawnow ; 
    
    if(IF_SAVE) 
        figdir = sprintf('./Figures/RatesDist/%s',dir) ; 
        fprintf('Writing %s \n',figdir) 
        try 
            mkdir(figdir) 
        end 
        
        ProcessFigure(fig, fullfile(figdir,figname), 2.2, [1.33*2.2, 2.2]) ; 
    end        

    hold off ;
end

for i=1:nbpop
    switch i 
      case 1
        figname='ChiDistE' ;
      case 2
        figname='ChiDistI' ;
      case 3
        figname='ChiDistS' ;
      case 4
        figname='ChiDistV' ;
    end
    
    if( ishandle( findobj('type','figure','name',figname) ) ) 
        fig = findobj('type','figure','name',figname) ; 
        figure(fig); hold on ; 
    else
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        xlabel('\chi_{II}') 
        ylabel('pdf') 
    end
    
    idx = ~isnan( khi(i,:) ) ; 
    X = khi(i,idx) ;
    X=X(X~=0) ;

    h = histogram( X , 10, 'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i},'EdgeAlpha',alp,'Linewidth',2) ; 
    drawnow ; 


    if strfind(dir,'L23') 
        if(i==2)
            xlim([0 20])
        else
            xlim([-30 0])
        end
    end

    if strfind(dir,'L5LAST') 
        if(i==3)
            xlim([0 20])
        else
            xlim([-50 0])
        end
    end

    if strfind(dir,'S1L5') 
        if(i==4)
            xlim([0 30])
        else
            xlim([-30 0])
        end
    end
    
    
    if(IF_SAVE) 
        figdir = sprintf('./Figures/ChiDist/%s_Iprtr%.3f', dir, Iprtr) ; 
        fprintf('Writing %s \n',figdir) 
        try 
            mkdir(figdir) 
        end 
        
        ProcessFigure(fig, fullfile(figdir,figname), 2.2, [1.33*2.2, 2.2]) ; 
    end

    hold off ; 
end 