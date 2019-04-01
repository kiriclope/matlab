function data = ImportData(model,nbpop,dir,file,n,k,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr,IF_PHI,PHI,IF_JabLoop,Jab)

    IF_TIMECOURSE = 0 ;

    popList = ['E' 'I' 'S' 'V'] ;
    if nargin<14
        IF_PHI=false ;
    end

    if nargin<16
        IF_JabLoop=false ;
    end
    
    if strcmp(model,'Rate') 
        dir = sprintf('%s/Threshold_Linear',dir) ; 
        separator = ';' ; 
    else 
        separator = ' ' ; 
    end 
    
    if(length(k)==1) 
        path = sprintf(['../%s/Simulations/%dpop/%s/N%d/K%d/g%.2f'],model,nbpop,dir,n,k,g) ;
    else
        path = sprintf(['../%s/Simulations/%dpop/%s/N%d/KE%dKI%d/g%.2f'],model,nbpop,dir,n,k(1),k(2),g) ;
    end 
    
    if strcmpi(IF_RING,'Ring') | strcmpi(IF_RING,'Gauss')  |  strcmpi(IF_RING,'Gauss2D')  | ...
            strcmpi(IF_RING,'Exp') | strcmpi(IF_RING,'Ring/Spec') | strcmpi(IF_RING,'Gauss/Spec')
        
        if(nbpop==1)
            if(isempty(Cff))
                RING_path = sprintf(['%s/CrecI%.4f'], IF_RING, Crec(1)) ;
            else
                RING_path = sprintf(['%s/CrecI%.4fCff%.4f'], IF_RING, Crec(1), Cff) ;
            end
        elseif(nbpop==2)

            if(length(Crec)==nbpop)
                RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCff%.4f'], IF_RING, Crec(1), Crec(2), Cff) ;
            else
                RING_path = sprintf(['%s/CrecEE%.4fCrecEI%.4fCrecIE%.4fCrecII%.4fCff%.2f'], IF_RING, Crec(1), Crec(2), Crec(3), Crec(4), Cff) ;
            end

        elseif(nbpop==3)

            if(length(Crec)==nbpop)
                RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.4fCff%.2f'], IF_RING,Crec(1), Crec(2), Crec(3), Cff) ;
            else
                RING_path = sprintf(['%s/CrecEE%.4fCrecEI%.4fCrecES%' ...
                                    '.2fCrecIE%.2fCrecII%.2fCrecIS%.2f' ...
                                    'CrecSE%.2fCrecSI%.2fCrecSS%.2fCff%.2f'], IF_RING, Crec(1), Crec(2), Crec(3), ... 
                                    Crec(4), Crec(5), Crec(6), Crec(7), Crec(8), Crec(9), Cff) ;
            end

        elseif(nbpop==4)

            if(length(Crec)==nbpop)
                if(IF_IEXT)
                    RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.4fCrecV%.4fCff%.4f'], IF_RING, Crec(1), Crec(2), Crec(3), Crec(4),Cff) ;
                else
                    RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.4fCrecV%.4f'], IF_RING, Crec(1), Crec(2), Crec(3), Crec(4)) ;
                end
            else
                RING_path = sprintf(['%s/CrecEE%.2fCrecEI%.2fCrecIE%.2fCrecII%.2fCff%.2f'], IF_RING, Crec(1), Crec(2), Crec(3), Crec(4), Cff) ;
            end

        end

        path = sprintf(['%s/%s'],path,RING_path) ;

    end

    switch IF_IEXT

      case 1
        if(nbpop>1)
            IEXT_path = sprintf(['DeltaPrtr/Iext_%s%.4f'],popList(nPrtr),Iprtr) ;
            %IEXT_path = sprintf(['Prtr_%s/Iext_%s%.4f'],popList(nPrtr),popList(nPrtr),Iprtr) ;
        else
            IEXT_path = sprintf(['DeltaPrtr/Iext_I%.4f'],Iprtr) ;
            % IEXT_path = sprintf(['Prtr_I/Iext_I%.4f'],Iprtr) ; 
        end                             
        path = sprintf(['%s/%s'],path,IEXT_path) ;
        
      case 2
        if(nbpop>1)
            IEXT_path = sprintf(['GaussPrtr/Iext_%s%.4f'],popList(nPrtr),Iprtr) ;
        else
            IEXT_path = sprintf(['GaussPrtr/Iext_I%.4f'],Iprtr) ; 
        end        
        path = sprintf(['%s/%s'],path,IEXT_path) ; 
        
      case 3 
        if(nbpop>1)
            IEXT_path = sprintf(['DeltaGaussPrtr/Iext_%s%.4f'],popList(nPrtr),Iprtr) ;
        else
            IEXT_path = sprintf(['DeltaGaussPrtr/Iext_I%.4f'],Iprtr) ; 
        end 
        path = sprintf(['%s/%s'],path,IEXT_path) ;
    end
    
    if IF_PHI
        PHI_path = sprintf(['PHI_%.2f'], PHI) ;
        path = sprintf(['%s/%s'],path,PHI_path) ;
    end

    if IF_TIMECOURSE
        path = sprintf(['%s/%s'],path,'TIMECOURSE') ;
    end

    if IF_JabLoop
        Jab_path = sprintf(['Jee%.4f'], Jab) ;
        path = sprintf(['%s/%s'],path, Jab_path) ;
    end
    
    file = sprintf(['%s/%s.txt'],path,file) ;
    fprintf('Reading %s \n',file)
    data = [] ;
    fid = fopen(file,'r') ;

    try
        data = importdata(file,separator) ; 
        fclose(fid) ;
    catch
        fprintf('DATA NOT FOUND\n') 
    end
    
end