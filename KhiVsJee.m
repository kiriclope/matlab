clear all ;
GlobalVars

IextBL = ExternalInput(model,nbpop,dir) ; 
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;


v_Jab = .1:.1:1 ; 
nbIdv = 10 ;

J = ImportJab(model,nbpop,dir) ; 
Iext = IextBL ; 
Iext(prtrPop) = IextBL(prtrPop) + Iprtr ; 
 
for i=1:length(v_Jab) 
    
    BL = ImportData(model, nbpop, dir, 'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, IextBL(prtrPop), 0, 0 , 1 , v_Jab(i) ) ;

    try
        for j=1:length(BL(1,:))-1
            IdvRatesBL(i,j) = mean( BL(:,j+1) ) ;
         end
     catch
         for j=1:nbN(nbpop)
             IdvRatesBL(i,j) = nan ; 
         end
     end

     for j=1:nbpop 
         mBL(j,i) = 0 ;
         RatesBL = IdvRatesBL(i, Cpt(j)+1:Cpt(j+1)) ; 
         % mBL(j,i) = mean(RatesBL) ; 
         [mBL(j,i) idx] = ratesCutOff(RatesBL, RatesBL, THRESHOLD, Cth, DIM, L) ;
     end                                 

     fprintf('Jab ') 
     fprintf('%.3f ', v_Jab(i)) 

     fprintf(' RatesBL ') 
     for j=1:nbpop         
         fprintf('%.3f ', mBL(j,i) )
     end    
     fprintf('\n')

     Prtr = ImportData(model, nbpop, dir, 'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, Iext(prtrPop), 0, 0 , 1 , v_Jab(i)) ;
     try
         for j=1:length(Prtr(1,:))-1
             IdvRatesPrtr(i,j) = mean(Prtr(:,j+1)) ;
         end
     catch
         for j=1:nbN(nbpop)
             IdvRatesPrtr(i,j) = nan ; 
         end
    end

    for j=1:nbpop 
        mPrtr(j,i) = 0 ;
        RatesPrtr = IdvRatesPrtr(i, Cpt(j)+1:Cpt(j+1)) ;
        % mPrtr(j,i) = mean(RatesPrtr) ;
        [mPrtr(j,i) idx] = ratesCutOff(RatesPrtr, RatesPrtr, THRESHOLD, Cth, DIM, L) ;
    end    

    fprintf(' RatesPrtr ') 
    for j=1:nbpop         
        fprintf('%.3f ', mPrtr(j,i) )
    end 
    fprintf('\n')
        

    fprintf(' Khi ') 
    for j=1:nbpop 
        khi(j,i) = ( mPrtr(j,i) - mBL(j,i) ) ./ Iprtr ; 
        fprintf('%.3f ', khi(j,i) ) 
    end 
    fprintf('\n') 

end

if(~IF_SAVE)
    figtitle = sprintf('%s_RatesVsJab',dir) ; 
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
    xlabel('J_{EE}')
    ylabel('Rates ')
end
    % xlim([.01 100])
    % set(gca,'Xscale', 'log')
    % ylim([.01 10])
    % set(gca,'Yscale', 'log') 

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
    plot(v_Jab, NormRates, 'd','Color',cl{i},'markersize',1) 

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


figtitle = sprintf('%s_KhiVsJab',dir) ; 
if(~IF_SAVE)

    if( ishandle( findobj('type','figure','name',figtitle) ) ) 
        fig = findobj('type','figure','name',figtitle) ; 
        fig = figure(fig); hold on ; 
    else
        xlabel('J_{EE}')
        ylabel('Norm. \chi ')
        plot(v_Jab, zeros(1,length(v_Jab)), '--','Color','k')
    end
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
    plot(v_Jab, NormKhi, '+','Color',cl{i},'markersize',1) 
    
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

