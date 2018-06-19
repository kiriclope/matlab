function SpatialProfie2D(model,nbpop,dir,Iext,K,n,g,ndI,dI,IF_RING,Crec,Cff,IF_Nk,IF_SAVE,Xlim,Ylim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Rate_InputDist_dI(nbpop,dir,Iext,K,dIlim,IF_DATA,n,k,g,IF_RING,Crec,Cff)
% Plots population time average rate vs the 
% perturbed input Iext(2)+dI.
% utils : Rate_InputDist(nbpop,dir,I,K,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    L = pi ;
    SCALING = 1.5 ;

    if nargin<10
        IF_RING = false ;
    end

    if nargin<13
        IF_Nk = false ;
    end

    if nargin<14
        IF_SAVE = false ;
    end

    THRESHOLD = 0 ; 
    popList = ['E' 'I' 'S' 'X'] ;
    if(strfind(dir,'2pop'))
        cl = {[1 0 0] [0 0 1] [0 1 1] } ;
        p = sscanf(dir,'2pop_p%f',[2,inf]) ;
    else
        p = [] ;
        cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
    end

    nbN = nbNeuron(nbpop,n,IF_Nk,p) ;
    Cpt = CptNeuron(nbpop,nbN) ;

    binLength = length(data(:,1)) ;
    neuronLength = length(data(1,:))-1 ;

    try
        data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ;
    end


    figname = 'Bump2D' ;
    fig = figure('Name',figname,'NumberTitle','off') ;
    set(gca,'nextplot','replacechildren') ; 

    for k=1:binLength

        for i=1:1

            Rates = data(k,Cpt(i)+2:Cpt(i+1)+1) ;
            Rates = reshape( Rates, sqrt( length( Rates ) ), sqrt( length( Rates ) ) );

            M = interp2(Rates,'cubic') ;                                    
            imagesc(M) ;

            h = colorbar ;
            xlabel('\theta_x ')
            ylabel('\theta_y ')
            
            xlim([0 length(M(:,1))])
            ylim([0 length(M(:,1))])

            drawnow ;
        end
        
    end

end