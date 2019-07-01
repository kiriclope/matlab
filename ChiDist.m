clear all ;
GlobalVars

IextBL = ExternalInput(model,nbpop,dir) ; 
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

Iprtr = .05 ;
Iext = IextBL ; 
Iext(prtrPop) = IextBL(prtrPop) + Iprtr ; 

file = '../LIF/Code/listL23.txt' ; 
data = [] ;
% try
%     data = importdata(file,'\n') ;
% catch
%     fprintf('DATA NOT FOUND\n') ;
% end

for i=1:100
    dirIdx = sprintf('%s/RND/%d',dir,i) ;
    
    try 
        % BL = ImportData(model, nbpop, dirIdx, 'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_IEXT, prtrPop, IextBL(prtrPop) ) ;  
        % for j=1:length(BL(1,:))-1
        %     IdvRatesBL(i,j) = mean( BL(:,j+1) ) ;
        % end
        
        % for j=1:nbpop 
        %     mBL(j,i) = nan ; 
        %     RatesBL = IdvRatesBL(i, Cpt(j)+1:Cpt(j+1)) ; 
        %     [mBL(j,i) idx] = ratesCutOff(RatesBL, RatesBL, THRESHOLD, Cth, DIM, L) ;
        % end                                 
        
        % fprintf(' RatesBL ') 
        % for j=1:nbpop         
        %     fprintf('%.3f ', mBL(j,i) )
        % end    
        % fprintf('\n')
        
        % Prtr = ImportData(model, nbpop, dirIdx, 'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_IEXT, prtrPop, Iext(prtrPop) ) ; 
        % for j=1:length(Prtr(1,:))-1
        %     IdvRatesPrtr(i,j) = mean(Prtr(:,j+1)) ;
        % end
        
        % for j=1:nbpop 
        %     mPrtr(j,i) = nan ;
        %     RatesPrtr = IdvRatesPrtr(i, Cpt(j)+1:Cpt(j+1)) ; 
        %     [mPrtr(j,i) idx] = ratesCutOff(RatesPrtr, RatesPrtr, THRESHOLD, Cth, DIM, L) ;
        % end 
        
        % fprintf(' RatesPrtr ') 
        % for j=1:nbpop 
        %     fprintf('%.3f ', mPrtr(j,i) ) 
        % end 
        % fprintf('\n') 

        BL = ImportData(model,nbpop,dirIdx,'Mean',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,IextBL(prtrPop)) ;        
        Prtr = ImportData(model,nbpop,dirIdx,'Mean',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ;        

        for j=1:nbpop 
            mBL(j,i) = mean(BL(:,j+1)) ;
            mPrtr(j,i) = mean(Prtr(:,j+1)) ;
            fprintf('%.2f ',mBL(j,i))
        end         
        fprintf('\n')
        
        if( all(mBL(:,i)>.01) && all(mPrtr(:,i)>.01) ) 
            for j=1:nbpop 
                khi(j,i) = ( mPrtr(j,i) - mBL(j,i) ) ;
            end 
                        
            if(any(data==i)) 
                for j=1:nbpop 
                    mBL(j,i) = nan ;
                    khi(j,i) = nan ;
                end
            end
            
            fprintf(' Rates ') 
            fprintf('%.3f ', mBL(:,i) ) 
            fprintf('\n') 

            fprintf(' Khi ') 
            fprintf('%.3f ', khi(:,i) ) 
            fprintf('\n') 
        else
            fprintf(' zero\n ') 
            
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
        ylabel('%') 
    end


    if strfind(dir,'L23') 
        if(i==2)
            xlim([0 15])
        else
            xlim([0 15])
        end
    end

    if strfind(dir,'L5') 
        if(i==3)
            xlim([0 20])
        else
            xlim([0 20])
        end
    end

    if strfind(dir,'S1L5') 
        if(i==4)
            xlim([0 15])
        else
            xlim([0 15])
        end
    end
    
    idx = ~isnan( mBL(i,:) ) ; 
    X = mBL(i,idx) ; 
    X=X(X~=0) ; 
    
    h = histogram( X , 8, 'Normalization', 'count' ,'DisplayStyle','stairs','EdgeColor',cl{i},'EdgeAlpha',alp,'Linewidth',2) ; 
    drawnow ; 
    
    if(IF_SAVE) 
        figdir = sprintf('./Paper/Robust/%s',dir) ; 
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
        xlabel('Activity Change(Hz)') 
        ylabel('%') 
    end
    
    idx = ~isnan( khi(i,:) ) ; 
    X = khi(i,idx) ;
    X=X(X~=0) ;

    h = histogram( X , 8, 'Normalization', 'count' ,'DisplayStyle','stairs','EdgeColor',cl{i},'EdgeAlpha',alp,'Linewidth',2) ; 
    drawnow ; 


    if strfind(dir,'L23') 
        if(i==2)
            xlim([0 2])
        else
            xlim([-2 0])
        end
    end

    if strfind(dir,'L5') 
        if(i==3)
            xlim([0 4])
        else
            xlim([-4 0])
        end
    end

    if strfind(dir,'S1L5') 
        if(i==4)
            xlim([0 3])
        else
            xlim([-3 0])
        end
    end

    % if strfind(dir,'S1L5') 
    %     if(i==4)
    %         xlim([0 30])
    %     else
    %         xlim([-30 0])
    %     end
    % end
    
    
    if(IF_SAVE) 
        figdir = sprintf('./Paper/Robust/%s_Iprtr%.3f', dir, Iprtr) ; 
        fprintf('Writing %s \n',figdir) 
        try 
            mkdir(figdir) 
        end 
        
        ProcessFigure(fig, fullfile(figdir,figname), 2.2, [1.33*2.2, 2.2]) ; 
    end

    hold off ; 
end 