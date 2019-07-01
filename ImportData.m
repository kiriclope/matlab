function data = ImportData(model,nbpop,dir,file,n,k,g,IF_RING,Crec,Cff,IF_Prtr,nPrtr,Xprtr)
    
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
    
    path = sprintf(['../%s/Simulations/%dpop/%s/N%d/K%d/g%.2f'],model,nbpop,dir,n,k,g) ; 

    switch IF_RING
        
      case 'SHARED'
        add_path = sprintf(['Shared%.2f/Cluster%.2f'], Crec(1), Crec(2)) ;
        path = sprintf(['%s/%s'],path,add_path) ;
        
      otherwise
        
        if(~isempty(IF_RING))
            if(nbpop==1)
                RING_path = sprintf(['%s/CrecI%.4f'], IF_RING, Crec(1)) ;
            elseif(nbpop==2)
                if(length(Crec)==nbpop)
                    RING_path = sprintf(['%s/CrecE%.4fCrecI%.4f'], IF_RING, Crec(1), Crec(2)) ;
                else
                    RING_path = sprintf(['%s/CrecEE%.4fCrecEI%.4fCrecIE%.4fCrecII%.4f'], IF_RING, Crec(1), Crec(2), Crec(3), Crec(4)) ;
                end
            elseif(nbpop==3)
                RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.4f'], IF_RING,Crec(1), Crec(2), Crec(3)) ;
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
        
      case {'Delta','Gauss','DeltaGauss'}
        % Iext = ExternalInput(model,nbpop,dir) ;
        % Xprtr1 = Xprtr + Iext(nPrtr) ;
        add_path = sprintf(['%sPrtr/Iext_%s%.4f'], IF_Prtr, popList(nPrtr), Xprtr) ;

      case 'PHI'
        add_path = sprintf(['PHI_%.2f'], Xprtr) ; 

      case 'TIMECOURSE'
        Iext = ExternalInput(model,nbpop,dir) ; 
        Xprtr1 = Xprtr + Iext(nPrtr) ; 
        add_path = sprintf(['DeltaPrtr/Iext_%s%.4f/TIMECOURSE'], popList(nPrtr), Xprtr1) ; 

      case 'JabLoop' 
        add_path = sprintf(['DeltaPrtr/Iext_%s%.4f/Jee%.4f'], Xprtr(1), Xprtr(2)) ;

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
        add_path = sprintf(['BenchMark/dt%.3f'], Xprtr) ;
         
      otherwise
        add_path = '' ;
    end
    
    path = sprintf(['%s/%s'],path,add_path) ;     
    file = sprintf(['%s/%s.txt'],path,file) ;
    fprintf('Reading %s \n',file) 

    data = [] ;
    try
        data = importdata(file,separator) ; 
    catch
        fprintf('DATA NOT FOUND\n') 
    end
    
end