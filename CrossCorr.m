function CrossCorr(Nfft,model,dir,x,y,AUTA_p,alpha,IF_RING) 

    cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
    
    if(nargin<8)
        if(nargin<1)
            Nfft = 1024 ;
        end
        IF_RING = 0 ;
    end

    str_Matrix = sprintf('../Cuda/Corr/xyCorr%d%d.dat',x,y) ;
    % if(IF_RING)
    %     str_Matrix = sprintf('../Cuda/Corr/LIF_xyCorr%d%d.dat',x,y) ;
    % else
    %     str_Matrix = sprintf('../Cuda/Corr/%s_xyCorr%d%d_%s_p%.2f.dat',model,x,y,dir,AUTA_p) ;
    % end
     
    fprintf('Import %s\n',str_Matrix) ;

    fMatrix = fopen(str_Matrix,'rb') ;
    M = fread(fMatrix,'*float') ;
    fclose(fMatrix) ;

    Nfft = 2.^nextpow2( floor( length(M)/10000) ) ;
    T = length(M)/Nfft ;

    M = reshape(M,Nfft,T).' ;
    
    figname=sprintf('%s_%s_AvgCrossCorr',model,dir) ; 
    if( ishandle( findobj('type','figure','name',figname) ) ) 
        fig = findobj('type','figure','name',figname) ; 
        fig = figure(fig); hold on ; 
    else
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
    end

    whos M ;

    AvgCorr = mean(M) ; 
    
    plot(0:Nfft/2-1, AvgCorr(1:Nfft/2), 'color', cl{x+1}) ; 
    patchline(0:Nfft/2-1, AvgCorr(1:Nfft/2),'linestyle','-','edgecolor',cl{x+1},'edgealpha',alpha,'linewidth',1.5)  
    xlabel('lag') 
    ylabel('<Cij>_N') 

   
    % if(IF_RING)
    %     str_Matrix = sprintf('../Cuda/Corr/AvgCorrRing.dat') ;
    % else 
    %     str_Matrix = sprintf('../Cuda/Corr/AvgCorr%d%d.dat',x,y) ;
    % end

    % fMatrix = fopen(str_Matrix,'rb') ;  
    % AvgCorr = fread(fMatrix,'*float') ; 
    % fclose(fMatrix) ;
     
    % plot(0:Nfft/2-1, AvgCorr(1:Nfft/2),'x','color',cl{x+1}) ;
    
    if(IF_RING)
       
        str_Idx = sprintf('../Cuda/Corr/xyIdxRing%d%d.txt',x,y) ;
        idx = importdata(str_Idx,' ') ; 
        
        X = zeros(T/2,Nfft/2) ; 
        counter = zeros(T/2,1) ; 
        for i=1:T/2
            for j=1:Nfft/2
                X(idx(i)+1,j) = X(idx(i)+1,j) + M(i,j) ; 
            end
            counter(idx(i)+1) = counter(idx(i)+1) + 1 ;
        end

        for i=1:T/2
            for j=1:Nfft/2
                X(i,j) = X(i,j) ./ counter(i) ;
            end
        end

        figname=sprintf('CrossCorrX') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ;

        for i=1:4
            plot(0:Nfft/2-1, X(i,1:Nfft/2),'x','color',cl{x+1}) ;
        end

        xlabel('lag') 
        ylabel('<Cij>_X')
        
        figname=sprintf('heatmap Cij_X') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
        
        xx=0:Nfft/2-1 ;
        imagesc(xx,xx,X(1:Nfft/2,:)) ;
        xlabel('lag') 
        ylabel('X')
        h = colorbar ;
        xlim([0 Nfft/2-1])
        ylim([0 Nfft/2-1])
        
    end

end