function figdir = FigDir(model,nbpop,dir,n,k,g,IF_RING,Crec,Cff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function FigDir(model,nbpop,dir,file,n,k,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

popList = ['E' 'I' 'S' 'X'] ;

if(length(k)==1)
    path = sprintf(['../%s/Figures/%dpop/%s/N%d/K%d/g%.2f'],model,nbpop,dir,n,k,g) ;
else
    path = sprintf(['../%s/Figures/%dpop/%s/N%d/KE%dKI%d/g%.2f'],model,nbpop,dir,n,k(1),k(2),g) ;
end 

if strcmpi(IF_RING,'Ring') | strcmpi(IF_RING,'Gauss') | strcmpi(IF_RING,'Ring2D')
    if(nbpop==1)
        RING_path = sprintf(['%s/CrecI%.2fCff%.2f'], IF_RING, Crec(1), Cff) ;
    elseif(nbpop==2)
        if(length(Crec)==nbpop)
            RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCff%.2f'],IF_RING, Crec(1), Crec(2), Cff) ;
        else
            RING_path = sprintf(['%s/CrecEE%.4fCrecEI%.4fCrecIE%.4fCrecII%.4fCff%.2f'], IF_RING, Crec(1), Crec(2), Crec(3), Crec(4), Cff) ;
        end
    elseif(nbpop==3)
        if(length(Crec)==nbpop)
            RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.2fCff%.2f'], IF_RING, Crec(1), Crec(2), Crec(3), Cff) ;
        else
            RING_path = sprintf(['%s/CrecEE%.2fCrecEI%.2fCrecES%.2fCrecIE%.2fCrecII%.2fCrecIS%.2fCrecSE%.2fCrecSI%.2fCrecSS%.2fCff%.2f'], ...
                                IF_RING, Crec(1), Crec(2), Crec(3), Crec(4), Crec(5), Crec(6), Crec(7), Crec(8), Crec(9), Cff) ;
        end
    elseif(nbpop==4)
        if(length(Crec)==nbpop)
            RING_path = sprintf(['%s/CrecE%.4fCrecI%.4fCrecS%.4fCrecV%.4fCff%.2f'], IF_RING, Crec(1), Crec(2), Crec(3), Crec(4),Cff) ;
        else
            RING_path = sprintf(['%s/CrecEE%.2fCrecEI%.2fCrecIE%.2fCrecII%.2fCff%.2f'], IF_RING,Crec(1), Crec(2), Crec(3), Crec(4), Cff) ;
        end
    end
    
    figdir = sprintf(['%s/%s'],path,RING_path) ;
else
    figdir = path ;
end


end