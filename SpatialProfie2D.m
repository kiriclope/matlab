function SpatialProfie2D(model,nbpop,dir,Iext,K,n,g,IF_RING,Crec,Cff,IF_Nk,IF_SAVE,Xlim,Ylim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Rate_InputDist_dI(nbpop,dir,Iext,K,dIlim,IF_DATA,n,k,g,IF_RING,Crec,Cff)
% Plots population time average rate vs the 
% perturbed input Iext(2)+dI.
% utils : Rate_InputDist(nbpop,dir,I,K,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    L = 2 ;
    SCALING = 1 ;
    nPrtr = 2 ;

    if( isempty(Iext) )
        Iext = ExternalInput(model,nbpop,dir) ;
        Iprtr = Iext(nPrtr) ;
    else
        dI = Iext ;
        Iext = ExternalInput(model,nbpop,dir) ;
        Iprtr = Iext(nPrtr) + dI ;
    end
    
    if nargin<8
        IF_RING = false ;
    end

    if nargin<11
        IF_Nk = false ;
    end

    if nargin<12
        IF_SAVE = false ;
    end

    THRESHOLD = 0 ; 
    popList = ['E' 'I' 'S' 'V'] ;

    if(strfind(dir,'2pop'))
        cl = {[1 0 0] [0 0 1] [0 1 1] } ;
        p = sscanf(dir,'2pop_p%f',[2,inf]) ;
    else
        p = [] ;
        cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
    end

    nbN = nbNeuron(nbpop,n,IF_Nk,p) ;
    Cpt = CptNeuron(nbpop,nbN) ;

    data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,1,nPrtr,Iprtr) ;
    baseline = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,1,nPrtr,Iext(nPrtr)) ;
    length(data(:,1))
    length(baseline(:,1))
    binLength = min( length(data(:,1)), length(baseline(:,1)) ) ;
    neuronLength = length(data(1,:))-1 ;

    % figname = 'Bump2D' ;
    % fig = figure('Name',figname,'NumberTitle','off') ;
    % set(gca,'nextplot','replacechildren') ; 

    for i=1:nbpop

        MeanRates = zeros( sqrt( nbN(i) ), sqrt( nbN(i) ) ) ;
        MeanBL = zeros( sqrt( nbN(i) ), sqrt( nbN(i) ) ) ;

        for k=1:binLength

            Rates = data(k,Cpt(i)+2:Cpt(i+1)+1) ;
            Rates = reshape( Rates, sqrt( length( Rates ) ), sqrt( length( Rates ) ) ) ;

            BL = baseline(k,Cpt(i)+2:Cpt(i+1)+1) ;
            BL = reshape( BL, sqrt( length( BL ) ), sqrt( length( BL ) ) ) ;

            MeanRates = MeanRates + Rates ;
            MeanBL = MeanBL + BL ;

            % M = interp2(Rates,'cubic') ;                                    
            % imagesc(M) ;

            % h = colorbar ;
            % xlabel('X')
            % ylabel('Y')
            
            % xlim([0 length(M(:,1))])
            % ylim([0 length(M(:,1))])

            % drawnow ;
        end

        MeanRates = MeanRates ./ MeanBL ;

        figname = sprintf('SpatialProfie_%s', popList(i));
        fig = figure('Name',figname,'NumberTitle','off') ;
        
        M = interp2(MeanRates./binLength,'cubic') ;          
        imagesc(M) ;
        
        h = colorbar ;
        xlabel('X (mm)')
        ylabel('Y (mm)')
        
        xlim([0 length(M(:,1))])
        ylim([0 length(M(:,1))])
        caxis([0 2]) 

        set(gca,'xtick',[0 length(M(:,1))/4 length(M(:,1))/2 3*length(M(:,1))/4 length(M(:,1))],'xticklabel',{-L/2, -L/4, 0, L/4, L/2})
        set(gca,'ytick',[0 length(M(:,1))/4 length(M(:,1))/2 3*length(M(:,1))/4 length(M(:,1))],'yticklabel',{L/2, L/4, 0, -L/4, -L/2})

    end

end