clear all ; 
GlobalVars 

Iext = ExternalInput(model,nbpop,dir) ; 

%function RatesVsN(model,nbpop,dir,Iext,K,file,n,g,IF_Nk,IF_RING,Crec,Cff,IF_IEXT,nPrtr,IF_SAVE)

IF_IEXT = '' ; 
l_N = [4] ; 
l_K = [500 1000 2000] ; 
IDX=5 ;
% CheckArgs(nargin) ;

cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
mk = {'x' 'o' '*' '.' 'x'} ;


IF_DATA = 1 ;
nbPoints = 0 ;

if(IF_DATA)
    for i=1:length(l_N)
        
        nbN = nbNeuron(nbpop,l_N(i),IF_Nk,[]) ;
        Cpt = CptNeuron(nbpop,nbN) ;

        for j=1:length(l_K)            

            data = ImportData(model,nbpop,dir,'IdvRates',l_N(i),l_K(j),g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iprtr) ;
            
            if(~isempty(data))
                nbPoints = nbPoints + 1 ;

                for k=1:nbpop
                    
                    for l=1:length(data(:,1))
                        PopAvgData(k,l) = mean(data(l, Cpt(k)+2:Cpt(k+1)+1 ) ) ;
                    end
                    
                    TpsAvgData(i,k,j) = mean( PopAvgData(k,IDX:end) ) ; 
                    
                end
            else
                
                for k=1:nbpop
                    TpsAvgData(i,k,j) = nan ;
                end
                
            end
        end
    end
else    
    nbPoints = 2 ;
end

if(nbPoints>1)

    fprintf('Finite K Rates \n')
    for i=1:length(l_K) 
        [u b] = RateInputDist(model,nbpop,dir,Iext*.01,l_K(i),1,[],0) ;
        for j=1:nbpop
            FiniteRates(j,i) = QchAvgTF(u(j),b(j))*1000 ;
            fprintf('%.3f ' , FiniteRates(j,i) )
        end
        fprintf('\n')
    end

    fprintf('MFRates \n')
    MFRates = BalRatesMF(model,nbpop,dir,Iext*.01,[])*1000 ; 
    fprintf('%.3f ' , MFRates )
    fprintf('\n')

    if(IF_DATA)
        
        for i=1:length(l_N)
            fprintf('Simulations \n')
            for j=1:length(l_K)
                for k=1:nbpop
                    fprintf('%.3f ' , TpsAvgData(i,k,j) ) 
                end
                fprintf('\n')
            end
        end
    end
    
    figname=sprintf('%sVsN',dir) ;
    fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

    for i=1:nbpop 
        
        NumRatesK = FiniteRates(i,:) ;
        % plot(l_K , ones(1,length(l_K)) .* NumRatesK ,'o' ,'color', cl{i}, 'markerSize', 5)
        % plot(1./sqrt(l_K) , ones(1,length(l_K)) .* NumRatesK - MFRates(i) ,'o' ,'color', cl{i}, 'markerSize', 5)
        % plot(l_K , ones(1,length(l_K)) .* MFRates(i) , 'color', cl{i}, 'markerSize', 5,'linestyle','--')
        % plot(1./sqrt(l_K) , ones(1,length(l_K)) .* MFRates(i) , 'color', cl{i}, 'markerSize', 5,'linestyle','--')

        % plot(sqrt(l_K) , ones(1,length(l_K)) .* NumRatesK ,'d' ,'color', cl{i}, 'markerSize', 5)
        % plot(sqrt(l_K) , ones(1,length(l_K)) .* MFRates(i) , 'color', cl{i}, 'markerSize', 5,'linestyle','--')
        
    end
    
    if(IF_DATA)
        for i=1:length(l_N)
            for j=1:nbpop 
                RatesK =  squeeze(TpsAvgData(i,j,:) ) ; 
                % plot(sqrt(l_K) ,  RatesK, mk{i}, 'color', cl{j}, 'markerSize', 5) 
                % plot(l_K ,  RatesK, mk{i}, 'color', cl{j}, 'markerSize', 5)
                plot(1./sqrt(l_K) ,  RatesK - MFRates(j), mk{i}, 'color', cl{j}, 'markerSize', 5)
            end
        end
    end
    
    xlabel('1/ $\sqrt{K}$','Interpreter','latex')
    % xlabel('$\sqrt{K}$','Interpreter','latex')
    ylabel(' Rates','Interpreter','latex')
    xlim([0 1./sqrt(l_K(1))])

    % figname=sprintf('%sVsN_2',file) ;
    % fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
    
    % for i=1:nbpop
    %     RatesK =  TpsAvgData(i,:) ; 
    %     plot(sqrt(l_K) ,  sqrt(l_K) .* RatesK, mk{1}, 'color', cl{i}, 'markerSize', 5)
    %     % plot(sqrt(l_K) ,  FiniteRates(i), '*-', 'color', cl{i} )
    %     plot(sqrt(l_K) , sqrt(l_K) .* MFRates(i), 'o', 'color', cl{i}, 'markerSize', 5)
    % end
    
    % xlabel('$\sqrt{K}$','Interpreter','latex')
    % ylabel('$\sqrt{K}$ Rates','Interpreter','latex')

    if(IF_SAVE)
        
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
        fprintf('Writing %s \n',figdir)
        
        try
            mkdir(figdir) ;
        end
        
        ProcessFigure(fig, fullfile(figdir,figname)) ;
    end

end

% function CheckArgs(nargin)

%     if(nargin<8)
%         IF_Nk = false ;
%     end

%     if(nargin<10)
%         IF_RING = false ;
%         Crec = 0 ;
%         Cff = 0 ;
%     end

%     if(nargin<13)
%         IF_IEXT=0 ;
%     end

%     switch IF_IEXT
%       case 0
%         Iext = ExternalInput(model,nbpop,dir) ; 
%         nPrtr = 2 ;
%         dI = Iext(nPrtr) ;
%         Iprtr = dI ;

%       otherwise 
%         if( isempty(Iext) )
%             Iext = ExternalInput(model,nbpop,dir) ;
%             dI = Iext(nPrtr) ;
%             Iprtr = dI ;
%         else
%             dI = Iext ;
%             Iext = ExternalInput(model,nbpop,dir) ;
%             Iprtr = Iext(nPrtr) + dI ;
%         end
%     end

%     if(nargin<15)
%         IF_SAVE=0 ;
%     end

% end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


