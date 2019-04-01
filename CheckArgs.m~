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