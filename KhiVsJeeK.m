clear all ;
GlobalVars

IextBL = ExternalInput(model,nbpop,dir) ; 
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

v_Jab = v_Jab(1):v_Jab(2):v_Jab(3) ; 
nbIdv = 25 ; 

J = ImportJab(model,nbpop,dir) ; 
Iext = IextBL ; 
Iext(prtrPop) = IextBL(prtrPop) + Iprtr ; 

u=[] ;
b=[] ;
uPrtr=[] ;
bPrtr=[] ;

for i=1:length(v_Jab) 
    J(1,1) = v_Jab(i) ; 

    MF_RatesBL = linsolve(g.*J,-g.*IextBL.') ; 
    MF_RatesPrtr = linsolve(g.*J,-g.*Iext.') ; 
    
    [u b] = RateInputDist(model,nbpop,dir,g.*IextBL,K,g,J,false,u,b) ; 
    RatesBL = QchAvgTF(u,b) ; 
    
    [uPrtr bPrtr] = RateInputDist(model,nbpop,dir,g.*Iext,K,g,J,false,uPrtr,bPrtr) ; 
    RatesPrtr = QchAvgTF(uPrtr,bPrtr) ; 

    DET(i) = det(J) ;
    G = J ./ 2 ;
    G(1,1) = J(1,1) ./ 10 ; 
    
    if(nbpop>3)
        G(4,2) = J(4,2) ./ 10 ;
        G(3,4) = J(3,4) ./ 10 ;
    end
    eigen = eig(G) ;

    fprintf('Jab ') 
    fprintf('%.3f ', v_Jab(i)) 

    fprintf('detJ ') 
    fprintf('%.3f ', DET(i)) 

    fprintf(' Rates ') 
    for j=1:nbpop 
        mBL(j,i) = RatesBL(j) ; 
        MF_m(j,i) = MF_RatesBL(j) ; 
        fprintf('%.3f ', mBL(j,i) ) 
    end 
    fprintf('\n') 

    fprintf(' Khi ') 
    for j=1:nbpop 
        lbd(j,i) = eigen(j) ;
        khi(j,i) = ( RatesPrtr(j) - RatesBL(j) ) ./ Iprtr ; 
        MF_khi(j,i) = ( MF_RatesPrtr(j) - MF_RatesBL(j) ) ./ Iprtr ;
        fprintf('%.3f ', khi(j,i) ) 
    end 
    fprintf('\n') 

end

% figtitle = sprintf('%s_DetJVsJab_MF',dir) ; 
% fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
% xlabel('J_{EE}')
% ylabel('Det(J) ')

% plot(v_Jab, DET, '-','Color','k') 

% figtitle = sprintf('%s_LbdVsJab_MF',dir) ; 
% fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
% xlabel('J_{EE}')
% ylabel('Re \lambda ')

% for i=1:nbpop 
%     plot(v_Jab, real( lbd(i,:) ), '-','Color',cl{i}) 
% end

if(~IF_SAVE)
    figtitle = sprintf('%s_RatesVsJab',dir) ; 
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
    xlabel('J_{EE}')
    ylabel('Rates ')
end

for i=1:nbpop

    if(IF_SAVE)
        if(i==1 || i==2) 
            figtitle = sprintf('%s_RatesVsJab_EI',dir) ; 
        else 
            figtitle = sprintf('%s_RatesVsJab_SV',dir) ; 
        end
        fig = popPerFig(i,dir,figtitle) ;
        xlabel('J_{EE}')
        ylabel('Rates ')
    end 

    NormRates = mBL(i,:) ;
    plot(v_Jab, NormRates, 'o','MarkerEdgeColor',cl{i}, 'markersize', 1, 'MarkerFaceColor', 'w') 
    plot(v_Jab, MF_m(i,:), '--','Color',cl{i}) 
    
    if(IF_SAVE) 
        if(i==2 || i==4) 
            figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ; 
            fprintf('Writing %s \n',figdir) 
            try
                mkdir(figdir) ;
            end
            ProcessFigure(fig, fullfile(figdir,figtitle)) ;
        end
    end

end

% if(IF_SAVE)
%     figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
%     fprintf('Writing %s \n',figdir)
%     try
%         mkdir(figdir) ;
%     end
%     ProcessFigure(fig, fullfile(figdir,figtitle)) ;
% end


if(~IF_SAVE)
    figtitle = sprintf('%s_KhiVsJab',dir) ; 
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
    xlabel('J_{EE}')
    ylabel('\chi ')
    plot(v_Jab, zeros(1,length(v_Jab)), '--','Color','k')
end

% xlim([.01 100])
% set(gca,'Xscale', 'log')
% ylim([.01 10])
% set(gca,'Yscale', 'log') 

for i=1:nbpop

    if(IF_SAVE)
        if(i==1 || i==2) 
            figtitle = sprintf('%s_KhiVsJab_EI',dir) ; 
        else 
            figtitle = sprintf('%s_KhiVsJab_SV',dir) ; 
        end
        fig = popPerFig(i,dir,figtitle) ; 
        xlabel('J_{EE}') 
        ylabel('\chi / m_{Baseline}') 
        plot(v_Jab, zeros(1,length(v_Jab)), '--','Color','k') 
    end

    NormKhi = khi(i,:) ./ mBL(i,:) ; 
    plot(v_Jab, NormKhi, 'o','MarkerEdgeColor',cl{i}, 'markersize', 1, 'MarkerFaceColor','w') 
    plot(v_Jab, MF_khi(i,:) ./ MF_m(i,:) , '--','Color',cl{i}) 
    %ylim([-.5 .5])
    if(IF_SAVE)
        if(i==2 || i==4)
            figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
            fprintf('Writing %s \n',figdir)
            try
                mkdir(figdir) ;
            end
            ProcessFigure(fig, fullfile(figdir,figtitle)) ;
        end
    end

end

% plot([0.05 0.05],[-2 2], '--','Color','k')
% plot([0.4 0.4],[-2 2], '--','Color','k')

% if(IF_SAVE)
%     figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
%     fprintf('Writing %s \n',figdir)
%     try
%         mkdir(figdir) ;
%     end
%     ProcessFigure(fig, fullfile(figdir,figtitle)) ;
% end

