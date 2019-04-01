clear all ;
GlobalVars

Iext = ExternalInput(model,nbpop,dir) ;
IextPrtr = Iext(prtrPop) + Iprtr ;

nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

data = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,IextPrtr);
baseline = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ;

binLength = min( length(data(:,1)), length(baseline(:,1)) ) ;
neuronLength = length(data(1,:))-1 ;

fprintf('Norm Rates ')

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
    MeanRates(find(MeanBL==0)) = 0 ; 

    fprintf('%.3f | ', mean(MeanRates(:)) ) ;

    figname = sprintf('SpatialProfie2D_I%s%.3f', popList(i),Iprtr);
    fig = figure('Name',figname,'NumberTitle','off') ;
    
    M = interp2(MeanRates, 'cubic') ;
    imagesc(M) ;
    
    xlabel('X (mm)')
    ylabel('Y (mm)')
    
    xlim([0 length(M(:,1))])
    ylim([0 length(M(:,1))])
    caxis([0 2]) 

    set(gca,'xtick',[0 length(M(:,1))/4 length(M(:,1))/2 3*length(M(:,1))/4 length(M(:,1))],'xticklabel',{-L/2, -L/4, 0, L/4, L/2})
    set(gca,'ytick',[0 length(M(:,1))/4 length(M(:,1))/2 3*length(M(:,1))/4 length(M(:,1))],'yticklabel',{L/2, L/4, 0, -L/4, -L/2})

    if(IF_SAVE)
        if(i==4)
            %h = colorbar ;
        end
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
        fprintf('Writing %s \n',figdir)
        try 
            mkdir(figdir) 
        end 
        ProcessFigure(fig, fullfile(figdir,figname), 2.5, [1.33*2.5, 2.5]) ;
    else
        h = colorbar ;
    end        
    hold off ;
    
end

fprintf('\n') ;