clear all ;
GlobalVars

Iext = ExternalInput(model,nbpop,dir) ; 
Iprtr = Iext(prtrPop) + Iprtr ;

nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

data = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,1,prtrPop,Iprtr) ; 

if(Iprtr~=0) 
    baseline = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,1,prtrPop,Iext(prtrPop)) ; 
else
    baseline = data ;
end

binLength = min( length(data(:,1)), length(baseline(:,1)) ) ;
neuronLength = length(data(1,:))-1 ;

fprintf('Norm Rates ')

for i=1:nbpop

    MeanRates = zeros(1, nbN(i) ) ;
    MeanBL = zeros(1,  nbN(i) ) ;
    
    for j=1:binLength
        Rates = data(j,Cpt(i)+2:Cpt(i+1)+1) ; 
        BL = baseline(j,Cpt(i)+2:Cpt(i+1)+1) ;

        MeanRates = MeanRates + Rates ; 
        MeanBL = MeanBL + BL ;
    end

    MeanRates = MeanRates ./ MeanBL ;
    MeanRates(find(MeanRates==Inf)) = 0 ;
    MeanRates(find(MeanRates==NaN)) = 0 ;

    fprintf('%.3f | ', mean(MeanRates(:)) ) ; 

    figname = sprintf('SpatialProfie_I%s%.3f', popList(i),Iprtr(prtrPop));
    fig = figure('Name',figname,'NumberTitle','off') ;
    
    
    h = colorbar ;
    xlabel('X (mm)')
    ylabel('Y (mm)')
    
    xlim([0 length(M(:,1))])
    ylim([0 length(M(:,1))])
    caxis([0 2]) 

    set(gca,'xtick',[0 length(M(:,1))/4 length(M(:,1))/2 3*length(M(:,1))/4 length(M(:,1))],'xticklabel',{-L/2, -L/4, 0, L/4, L/2})
    set(gca,'ytick',[0 length(M(:,1))/4 length(M(:,1))/2 3*length(M(:,1))/4 length(M(:,1))],'yticklabel',{L/2, L/4, 0, -L/4, -L/2})
    
end

fprintf('\n') ;