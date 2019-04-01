function RatesVsN(model,nbpop,dir,Iext,K,file,n,g,IF_Nk,IF_RING,Crec,Cff,IF_IEXT,nPrtr,IF_SAVE)

    l_N = [8 16] ;
    l_K = [250 500 1000 2000] ;

    CheckArgs(nargin) ;

    cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
    mk = {'+' 'o' '*' '.' 'x'} ;


    nbPoints = 0 ;

    for i=1:length(l_N)

        nbN = nbNeuron(nbpop,l_N(i),IF_Nk,[]) ;
        Cpt = CptNeuron(nbpop,nbN) ;

        for j=1:length(l_K)            

            try
                data = ImportData(model,nbpop,dir,file,l_N(i),l_K(j),g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ;
                ERROR = 0 ;
            catch
                fprintf('DATA NOT FOUND\n') ;
                ERROR = 1 ;
            end

            if(ERROR==0)
                
                nbPoints = nbPoints + 1 ;

                for k=1:nbpop
                    
                    % for l=2:length(data(1,:))
                    %     IdvData(k) = mean(data(:,k)) ;
                    % end
                   
                    for l=1:length(data(:,1))
                        PopAvgData(k,l) = mean(data(l, Cpt(k)+2:Cpt(k+1)+1 ) ) ;
                    end
                    
                    TpsAvgData(k,i,j) = mean( PopAvgData(k,:) ) ;

                end
            else

                for k=1:nbpop
                    TpsAvgData(k,i,j) = nan ;
                end
                
            end
        end
    end

    if(nbPoints>1)

        fprintf('MFRates \n')
        for i=1:length(l_K)            
            [u b] = RateInputDist(model,nbpop,dir,Iext,l_K(i),1,[],0) ;
            for j=1:nbpop
                MFRates(j,i) = QchAvgTF(u(j),b(j)) ;
                fprintf('%.3f ' , MFRates(j,i) )
            end
            fprintf('\n')
        end

        fprintf('Simulations \n')
        for j=1:length(l_K)
            for i=1:nbpop
                fprintf('%.3f ' , TpsAvgData(i,j) )
            end
            fprintf('\n')
        end


        figname=sprintf('%sVsN',file) ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        
        for i=1:nbpop
            for j=1:length(l_N) 
                plot(sqrt(l_K) , sqrt(l_K) * TpsAvgData(i,j), mk{j}, 'color', cl{i}, 'markerSize', 5)
            end
            plot(sqrt(l_K) , sqrt(l_K) * MFRates(i), '--', 'color', cl{i} )
        end
        
        xlabel('$\sqrt{K}$','Interpreter','latex')
        ylabel('$\sqrt{K}$ Rates','Interpreter','latex')

        if(IF_SAVE)
            
            figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
            fprintf('Writing %s \n',figdir)
            
            try
                mkdir(figdir) ;
            end
            
            ProcessFigure(fig, fullfile(figdir,figname)) ;
        end

    end
    
    function CheckArgs(nargin)
        
        if(nargin<8)
            IF_Nk = false ;
        end
        
        if(nargin<10)
            IF_RING = false ;
            Crec = 0 ;
            Cff = 0 ;
        end
        
        if(nargin<13)
            IF_IEXT=0 ;
        end
        
        switch IF_IEXT
          case 0
            Iext = ExternalInput(model,nbpop,dir) ; 
            nPrtr = 2 ;
            dI = Iext(nPrtr) ;
            Iprtr = dI ;
            
          otherwise 
            if( isempty(Iext) )
                Iext = ExternalInput(model,nbpop,dir) ;
                dI = Iext(nPrtr) ;
                Iprtr = dI ;
            else
                dI = Iext ;
                Iext = ExternalInput(model,nbpop,dir) ;
                Iprtr = Iext(nPrtr) + dI ;
            end
        end
        
        if(nargin<15)
            IF_SAVE=0 ;
        end
    
    end 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = phi(x)
        out = exp(-x.^2./2)./sqrt(2.*pi);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = Phi(x)
        out = .5.*(1+erf( x./sqrt(2) ) ) ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = QchAvgTF(u,b)
        if(b>0)
            out = u.*Phi( u./sqrt(b) ) + sqrt(b).*phi( u./sqrt(b) ) ;
        else
            out = u ;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
