GlobalVars

Iext = ExternalInput(model,nbpop,dir) ;
Iprtr = .5 ;
Iprtr = Iext(prtrPop) + Iprtr ;
    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

data = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,1,nPrtr,Iprtr) ;
baseline = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,1,nPrtr,Iext(prtrPop)) ;

binLength = min( length(data(:,1)), length(baseline(:,1)) ) ;
neuronLength = length(data(1,:))-1 ;

for i=1:nbpop

    MeanRates = zeros( sqrt( nbN(i) ), sqrt( nbN(i) ) ) ;
    MeanBL = zeros( sqrt( nbN(i) ), sqrt( nbN(i) ) ) ;

    for j=1:binLength
        Rates = data(j,Cpt(i)+2:Cpt(i+1)+1) ;
        Rates = reshape( Rates, sqrt( length( Rates ) ), sqrt( length( Rates ) ) ) ;
        
        BL = baseline(j,Cpt(i)+2:Cpt(i+1)+1) ;
        BL = reshape( BL, sqrt( length( BL ) ), sqrt( length( BL ) ) ) ;

        MeanRates = MeanRates + Rates ;
        MeanBL = MeanBL + BL ;
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