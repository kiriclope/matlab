function data = ImportData(model,nbpop,dir,file,n,k,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr,IF_PHI,PHI)

    popList = ['E' 'I' 'S' 'X'] ;
    if nargin<14
        IF_PHI=false ;
    end
    
    if strcmp(model,'Rate')
        dir = sprintf('%s/Threshold_Linear',dir) ;
        separator = ';' ;
    end
    if strcmp(model,'Binary')
        separator = ';' ;
    end
    if strcmp(model,'IF') | strcmp(model,'Bump') | strcmp(model,'LIF') | strcmp(model,'Rates')
        separator = ' ' ;
    end

    if(length(k)==1)
        path = sprintf(['../%s/Simulations/%dpop/%s/N%d/K%d/g%.2f'],model,nbpop,dir,n,k,g) ;
    else
        path = sprintf(['../%s/Simulations/%dpop/%s/N%d/KE%dKI%d/g%.2f'],model,nbpop,dir,n,k(1),k(2),g) ;
    end 

    if IF_RING & strcmp(model,'Binary')
        
        if(nbpop>1 & length(Crec)==1)
            RING_path = sprintf(['Ring/Crec%.2fCff%.2f'], Crec, Cff) ;
        elseif(nbpop==1)
            RING_path = sprintf(['Ring/CrecI%.2fCff%.2f'], Crec(1), Cff) ;
        elseif(nbpop==2)
            if(length(Crec)==nbpop)
                RING_path = sprintf(['Ring/CrecE%.2fCrecI%.2fCff%.2f'], Crec(1), Crec(2), Cff) ;
            else
                RING_path = sprintf(['CrecEE%.2fCrecEI%.2fCrecIE%.2fCrecII%.2fCff%.2f'], Crec(1), Crec(2), Crec(3), Crec(4), Cff) ;
            end
        end
        
        path = sprintf(['%s/%s'],path,RING_path) ;
    end
    
    if strcmpi(IF_RING,'Ring') | strcmpi(IF_RING,'Gauss') | ...
            strcmpi(IF_RING,'Ring2D') | strcmpi(IF_RING,'Exp') | strcmpi(IF_RING,'Specific')
        if(nbpop==1)
            RING_path = sprintf(['%s/CrecI%.2fCff%.2f'], IF_RING, Crec(1), Cff) ;
        elseif(nbpop==2)
            if(length(Crec)==nbpop)
                RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCff%.4f'], ...
                                    IF_RING, Crec(1), Crec(2), Cff) ;
            else
                RING_path = sprintf(['%s/CrecEE%.4fCrecEI%.4fCrecIE%.4fCrecII%.4fCff%.2f'], IF_RING,Crec(1), Crec(2), Crec(3), Crec(4), Cff) ;
            end
        elseif(nbpop==3)
            if(length(Crec)==nbpop)
                RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.4fCff%.2f'], IF_RING,Crec(1), Crec(2), Crec(3), Cff) ;
            else
                RING_path = sprintf(['%s/' ...
                                    'CrecEE%.4fCrecEI%.4fCrecES%' ...
                                    '.2fCrecIE%.2fCrecII%.2fCrecIS%.2fCrecSE%.2fCrecSI%.2fCrecSS%.2fCff%.2f'], ...
                                    IF_RING,Crec(1), Crec(2), Crec(3), ...
                                    Crec(4), Crec(5), Crec(6), ...
                                    Crec(7), Crec(8), Crec(9), Cff) ;
            end
        elseif(nbpop==4)
            if(length(Crec)==nbpop)
                RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.4fCrecV%.4fCff%.2f'], IF_RING,Crec(1), Crec(2), Crec(3), Crec(4),Cff) ;
            else
                RING_path = sprintf(['%s/CrecEE%.2fCrecEI%.2fCrecIE%.2fCrecII%.2fCff%.2f'], IF_RING,Crec(1), Crec(2), Crec(3), Crec(4), Cff) ;
            end
        end

        path = sprintf(['%s/%s'],path,RING_path) ;
    end

    if IF_RING==3 & strcmp(model,'IF')
        % RING_path = sprintf(['Space2D/Crec%.2fCff%.2f'], Crec, Cff) ;
        if(nbpop==1)
            RING_path = sprintf(['Space2D/CrecI%.2fCff%.2f'], Crec(1), Cff) ;
        elseif(nbpop==2)
            if(length(Crec)==nbpop)
                RING_path = sprintf(['Space2D/CrecE%.2fCrecI%.2fCff%.2f'], Crec(1), Crec(2), Cff) ;
            else
                RING_path = sprintf(['Space2D/CrecEE%.2fCrecEI%.2fCrecIE%.2fCrecII%.2fCff%.2f'], Crec(1), Crec(2), Crec(3), Crec(4), Cff) ;
            end
        end
        path = sprintf(['%s/%s'],path,RING_path) ;
    end

    if IF_IEXT
        if(nbpop>1)
            IEXT_path = sprintf(['Prtr_%s/Iext_%s%.4f'],popList(nPrtr),popList(nPrtr),Iprtr) ;
        else
            IEXT_path = sprintf(['Prtr_I/Iext_I%.4f'],Iprtr) ; 
        end
        
        path = sprintf(['%s/%s'],path,IEXT_path) ;
    end

    if IF_PHI
        PHI_path = sprintf(['PHI_%.2f'], PHI) ;
        path = sprintf(['%s/%s'],path,PHI_path) ;
    end
    
    file = sprintf(['%s/%s.txt'],path,file) ;
    fprintf('Reading %s \n',file)
    data = [] ;
    fid = fopen(file,'r') ;
    data = importdata(file,separator) ;
    fclose(fid) ;

end