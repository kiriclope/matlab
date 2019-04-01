clear all ;
GlobalVars

IextBL = ExternalInput(model,nbpop,dir) ; 
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

v_Iprtr = v_Iprtr(1):v_Iprtr(2):v_Iprtr(3) ;
nbIdv = 25 ;

J = ImportJab(model,nbpop,dir) ;
Iext = IextBL ;

for i=1:length(v_Iprtr)
    Iext(2) = IextBL(2) + v_Iprtr(i) ;
    Rates = linsolve(J,-Iext.') ; 

    % [u b] = RateInputDist(model,nbpop,dir,Iext,100000,1,[],false) ;
    % Rates = QchAvgTF(u,b) ;

    fprintf('I_opto ') 
    fprintf('%.3f ', v_Iprtr(i)) 

    fprintf(' Rates ') 
    for j=1:nbpop 
        m(j,i) = Rates(j) ;
        fprintf('%.3f ', m(j,i) )
    end 
    fprintf('\n')
end

figtitle = sprintf('%s_RatesVsIopto_MF',dir) ; 
fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
xlabel('I_{opto}')
ylabel('Norm. Rates')

    % xlim([.01 100])
    % set(gca,'Xscale', 'log')
    % ylim([.01 10])
    % set(gca,'Yscale', 'log') 

for i=1:nbpop
    
    NormRates = m(i,:) ./ m(i,1) ; 
    % plot(v_Iprtr, NormRates, 'o','MarkerEdgeColor',cl{i},'MarkerSize',2,'MarkerFaceColor','none','LineWidth', 1)
    plot(v_Iprtr, NormRates, '-','Color',cl{i})            
end

if(IF_SAVE)
    figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
    fprintf('Writing %s \n',figdir)
    try
        mkdir(figdir) ;
    end
    ProcessFigure(fig, fullfile(figdir,figtitle)) ;
end

