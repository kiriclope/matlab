function data = ImportData(model,nbpop,dir,file,n,k,g,IF_RING,Crec,Cff,IF_Prtr,nPrtr,Xprtr,IF_Dij,Dij)

% IF_Dij = 1 ; 
% Dij = [1.0 1.5 0.5 1.0] ; 
% Dij = [1.0 1.0 0.6 1.0] ; 
% Dij = [0.0 1.0 1.0 0.0] ; 
% Dij = [1.0,1.0,0.7, 1.0,1.0,1.0, 1.0,1.0,1.0] ; 

    if(nargin<12) 
        nPrtr = 2 ;
        if(nargin<10)
            IF_Prtr = '' ;
            Xprtr = [] ;
            if(nargin<8)
                IF_RING = '' ;
                Crec = [] ;
                Cff = [] ;
            end
        end
    end

    popList = ['E' 'I' 'S' 'V'] ;
    if(nbpop==1)
        popList = 'I' ;
    end

    if strcmp(model,'Rate') 
        dir = sprintf('%s/Threshold_Linear',dir) ; 
        separator = ';' ; 
    else 
        separator = ' ' ; 
    end 
    
    if strcmp(model,'Connectivity') 
        path = sprintf(['../%s/%dpop/N%d/K%d'],model,nbpop,n,k) ;
    else
        path = sprintf(['../%s/Simulations/%dpop/%s/N%d/K%d/g%.2f'],model,nbpop,dir,n,k,g) ; 
    end

    switch IF_RING
        
      case 'SHARED'
        add_path = sprintf(['Shared%.2f/Cluster%.2f'], Crec(1), Crec(2)) ;
        path = sprintf(['%s/%s'],path,add_path) ;
        
      otherwise
        
        if(~isempty(IF_RING))
            if(nbpop==1)
                RING_path = sprintf(['%s/CrecI%.4f'], IF_RING, Crec(1)) ;
            elseif(nbpop==2)
                if(IF_Dij)
                    RING_path = sprintf( ['%s/' ...
                                        'CrecEE%.4fCrecEI%.4f' ...
                                        'CrecIE%.4fCrecII%.4f'], IF_RING, ...
                                        Crec(1)*Dij(1), Crec(2)*Dij(2), ...
                                        Crec(1)*Dij(3), Crec(2)*Dij(4) ) ; 
                else
                    RING_path = sprintf(['%s/CrecE%.4fCrecI%.4f'], IF_RING, Crec(1), Crec(2)) ;
                end
            elseif(nbpop==3)
                if(IF_Dij)
                    RING_path = sprintf( ['%s/' ...
                                        'CrecEE%.4fCrecEI%.4fCrecES%.4f' ...
                                        'CrecIE%.4fCrecII%.4fCrecIS%.4f' ...
                                        'CrecSE%.4fCrecSI%.4fCrecSS%.4f'], IF_RING, ...
                                         Crec(1)*Dij(1), Crec(2)*Dij(2), Crec(3)*Dij(3), ...
                                         Crec(1)*Dij(4), Crec(2)*Dij(5), Crec(3)*Dij(6), ...
                                         Crec(1)*Dij(7), Crec(2)*Dij(8), Crec(3)*Dij(9) ) ;
                    
                else
                    RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.4f'], IF_RING,Crec(1), Crec(2), Crec(3)) ;
                end
            elseif(nbpop==4)
                RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.4fCrecV%.4f'], IF_RING, Crec(1), Crec(2), Crec(3), Crec(4)) ;
            end
            
            path = sprintf(['%s/%s'],path,RING_path) ;

            if(~isempty(Cff))
                path = sprintf(['%sCff%.4f'],path,Cff) ;
            end

        end
    end

    switch IF_Prtr
        
      case {'Cos','Delta','Gauss','DeltaGauss','Exp'}
        if(nPrtr>0)
            add_path = sprintf(['%sPrtr/Iext_%s%.4f'], IF_Prtr, popList(nPrtr), Xprtr(nPrtr)) ;
        else
            add_path = sprintf(['%sPrtr/Iext_All%.4f'], IF_Prtr, Xprtr(1)) ;
        end
      case 'Prop'
        add_path = sprintf(['Prop%.2f/Iext_%s%.4f'], Cff, popList(nPrtr), Xprtr(nPrtr)) ;

      case 'PropWeak'
        add_path = sprintf(['PropWeak%.2f/Iext_%s%.4f'], Cff, popList(nPrtr), Xprtr(nPrtr)) ;
        
      case 'PHI'
        add_path = sprintf(['PHI_%.2f'], Xprtr) ; 

      case 'TIMECOURSE'
        Iext = ExternalInput(model,nbpop,dir) ; 
        Xprtr1 = Xprtr(nPrtr) ; 
        add_path = sprintf(['DeltaPrtr/Iext_%s%.4f/TIMECOURSE'], popList(nPrtr), Xprtr1) ; 

      case 'JabLoop' 
        add_path = sprintf(['DeltaPrtr/Iext_%s%.4f/Je0%.4fJee%.4f'], ...
                           popList(nPrtr), Xprtr(1), Xprtr(2), Xprtr(3)) ;

      case 'AUTA'
        if(nPrtr==0)
            if(nbpop==1)
                add_path = sprintf(['AutaI%.2f'], Xprtr) ; 
            else
                add_path = sprintf(['AutaE%.2f'], Xprtr) ; 
            end
        elseif(nPrtr==1)
            add_path = sprintf(['AutaI%.2f'], Xprtr) ; 
        elseif(nPrtr==2) 
            add_path = sprintf(['AutaE%.2fI%.2f'], Xprtr(1), Xprtr(2) ) ; 
        end
        
      case 'BENCH'
        add_path = sprintf(['BenchMark/dt%.4f'], Xprtr) ;

      case 'COND'
        add_path = sprintf(['rho%.2f'], Xprtr) ;
        
      otherwise
        add_path = '' ;
    end
    
    path = sprintf(['%s/%s'],path,add_path) ;     
    if strcmp(model,'Connectivity') 
        file = sprintf(['%s/%s.dat'],path,file) ;
        data = file ; 
    else
        file = sprintf(['%s/%s.txt'],path,file) ;

        fprintf('Reading %s \n',file) 
        
        data = [] ;
        try
            data = importdata(file,separator) ; 
        catch
            fprintf('DATA NOT FOUND\n') 
        end
    end
    
end