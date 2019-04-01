GlobalVars

Iext = ExternalInput(model,nbpop,dir) ;    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

v_Iprtr = v_Iprtr(1):v_Iprtr(2):v_Iprtr(3) ;
nbIdv = 25 ;

for i=1:length(v_Iprtr)   
    data = ImportData(model, nbpop, dir, 'MeanInputs', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, Iext(prtrPop) + v_Iprtr(i) ) ;
    try
        for j=1:length(data(1,:))-1
            IdvInputs(j) = mean(data(:,j+1)) ;
        end
    catch
        for j=1:nbN(nbpop)
            IdvInputs(j) = nan ;
        end
    end
        
    fprintf('Synaptic Inputs \n')
    for j=1:nbpop 
        for k=1:nbpop 
            avgIsyn(j,k,i) = 0 ;
            Isyn = IdvInputs( nbpop*Cpt(j)+1+(k-1)*nbN(j):nbpop*Cpt(j)+k*nbN(j) ) ; 
            
            avgIsyn(j,k,i) = mean(Isyn) ;
            %[avgIsyn(k,j,i) idx] = ratesCutOff(Isyn, THRESHOLD, Cth, DIM, L) ;
                         
            fprintf('%.3f ', avgIsyn(j,k,i) )
        end
        fprintf('\n')
    end
    fprintf('\n')
end

for i=1:nbpop

    figtitle = sprintf('%s_IsynVsIopto_%s',dir,popList(i)) ;
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
    xlabel('I_{opto}')
    ylabel('Norm. Isyn.')

    IsynE = zeros(1,length(v_Iprtr)) ;
    IsynI = zeros(1,length(v_Iprtr)) ;
    IsynNet = zeros(1,length(v_Iprtr)) ;

    for j=1:nbpop        
        normIsyn = zeros(1,length(v_Iprtr)) ;
        for k=1:length(v_Iprtr)
            normIsyn(k) = avgIsyn(i,j,k)./avgIsyn(i,j,1) ;
            IsynNet(k) = IsynNet(k) + avgIsyn(i,j,k) ;
            
            if(j==1)
                IsynE(k) = avgIsyn(i,j,k) ;
            else
                IsynI(k) = IsynI(k) + avgIsyn(i,j,k) ;
            end
        end
        % plot(v_Iprtr, normIsyn, 'o','MarkerEdgeColor',cl{j},'MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1)
        % plot(v_Iprtr, normIsyn, '-','Color',cl{j})
    end

    plot(v_Iprtr, IsynE./IsynE(1), 'o','MarkerEdgeColor','r','MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1)
    plot(v_Iprtr, IsynE./IsynE(1), '--','Color','r')

    plot(v_Iprtr, IsynI./IsynI(1), 'o','MarkerEdgeColor','b','MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1)
    plot(v_Iprtr, IsynI./IsynI(1), '--','Color','b')

    plot(v_Iprtr, IsynNet./IsynNet(1), 'o','MarkerEdgeColor','k','MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1)
    plot(v_Iprtr, IsynNet./IsynNet(1), '--','Color','k')
    
    % xlim([.01 100])
    % set(gca,'Xscale', 'log')
    % ylim([.01 10])
    % set(gca,'Yscale', 'log') 
         
    if(IF_SAVE & (i==2 | i==4))                
        figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff,IF_DATA) ;
        fprintf('Writing %s \n',figdir)                
        try
            mkdir(figdir) ;
        end                
        ProcessFigure(fig, fullfile(figdir,figname)) ;
    end
    
    hold off ;
end

