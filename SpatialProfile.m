clear all ;
GlobalVars

Iext = ExternalInput(model,nbpop,dir) ; 
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ; 
Cpt = CptNeuron(nbpop,nbN) ; 

if(Iprtr~=0) 
    Iprtr = Iext(prtrPop) + Iprtr ; 
    data = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iprtr) ; 
    baseline = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ; 
else
    data = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ;
    baseline = data ;
end

binLength = min( length(data(:,1)), length(baseline(:,1)) ) ;
neuronLength = length(data(1,:))-1 ; 

fprintf('Norm Rates ')

for i=1:nbpop

    MeanRates = zeros(1, nbN(i) ) ;     
    for j=1:binLength
        Rates = data(j,Cpt(i)+2:Cpt(i+1)+1) ; 
        BL = baseline(j,Cpt(i)+2:Cpt(i+1)+1) ; 
        MeanRates = MeanRates + Rates ; 
    end

    if(Iprtr~=0)
        MeanBL = zeros(1,  nbN(i) ) ;     
        for j=1:binLength
            BL = baseline(j,Cpt(i)+2:Cpt(i+1)+1) ; 
            MeanBL = MeanBL + BL ; 
        end
        MeanRates = MeanRates ./ MeanBL ;
    else
        MeanRates = MeanRates ./ binLength ;
    end
    
    MeanRates(find(MeanRates==Inf)) = 0 ;
    MeanRates(find(MeanRates==NaN)) = 0 ;
    MeanRates(find(MeanRates==-Inf)) = 0 ;

    fprintf('%.3f | ', mean(MeanRates(:)) ) ; 

    figtitle = sprintf('SpatialProfie_I%s%.3f', popList(i),Iprtr);

    if( ishandle( findobj('type','figure','name',figtitle) ) )
        fig = findobj('type','figure','name',figtitle) ; 
        fig = figure(fig); hold on ; 
    else
        fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ;
        xlabel('X (mm)') 
        ylabel('Rates') 
    end
    X = linspace(-L/2,L/2,nbN(i)) ; 
    idx = 1:nbN(i) ;

    fit = smooth(idx,MeanRates,.05,'lowess') ;
    plt = patchline(X,fit,'linestyle','-','edgecolor',cl{i},'edgealpha',.25,'linewidth',2) ; 
end

fprintf('\n') ;

