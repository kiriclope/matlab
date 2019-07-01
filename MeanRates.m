clear all ; 
GlobalVars 

Iext = ExternalInput(model,nbpop,dir) ; 
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ; 
Cpt = CptNeuron(nbpop,nbN) ; 

Iprtr = .0 ;
Iext(prtrPop) = Iext(prtrPop) + Iprtr ;
tps = [] ;
for Idx=[24 42 67 65 75 76 77 83 92 93 ]
    dirIdx = sprintf('%s/RND/%d',dir,Idx) ; 
    data = [] ;
    % data = ImportData(model,nbpop,dirIdx,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ;  
    data = ImportData(model,nbpop,dirIdx,'Mean',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ;
        
    try
        % for i=2:length(data(1,:))       
        %     IdvRates(i) = mean(data(:,i)) ;
        % end
        tps =  data(:,1)./1000 ;
        
        % for i=1:nbpop 
        %     for j=1:length(data(:,1)) 
        %         PopRate(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1)+1)) ; 
        %     end 
        % end         
        
        for i=1:nbpop 
            for j=1:length(data(:,1))
                PopRate(i,j) = data(j,i+1) ;
            end         
        end
    
        if(~isempty(data))
            figname=sprintf('MeanRates_%s',dirIdx) ; 
            fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
            
            for i=1:nbpop
                fprintf('%.2f %.2f\n',length(tps),length(PopRate(i,:))) 
                plot(tps,PopRate(i,1:length(tps)),'color',cl{i}) 
            end 
            xlabel('t (s)') 
            ylabel('Activities (Hz)') 
            xlim([0 10])
            drawnow ; 
            hold off ; 
        end
        
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