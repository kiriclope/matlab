clear all ; 
GlobalVars 

Iext = ExternalInput(model,nbpop,dir) ; 
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ; 
Cpt = CptNeuron(nbpop,nbN) ; 

Iprtr = 0 ;
Iext(prtrPop) = Iext(prtrPop) + Iprtr ;

for Idx=1:100
    dirIdx = sprintf('%s_RND_%d',dir,Idx) ; 

    try
        data = ImportData(model,nbpop,dirIdx,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ; 
        for i=2:length(data(1,:)) 
            IdvRates(i) = mean(data(:,i)) ;
        end
        tps = data(:,1)/1000 + 2 ;
        for i=1:nbpop 
            for j=1:length(data(:,1)) 
                PopRate(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1)+1)) ; 
            end 
        end 

        figname=sprintf('MeanRates_%s',dirIdx) ; 
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        
        for i=1:nbpop
            plot(tps,PopRate(i,:),'color',cl{i}) 
        end 
        xlabel('t (s)') 
        ylabel('Activities (Hz)') 
        
        drawnow ; 
        hold off ; 
        
        if(IF_SAVE)
            figdir = sprintf('./Figures/MeanRates/%s_Iprtr%.3f', dir, Iprtr) ; 
            fprintf('Writing %s \n',figdir)
            try 
                mkdir(figdir) 
            end 
            ProcessFigure(fig, fullfile(figdir,figname), 2.2, [1.33*2.2, 2.2]) ;
        end        
        hold off ; 
    catch 
        fprintf('ERROR \n') ;
    end
end