function J = ImportJab(model,nbpop,dir,IF_SYMB)
    
    if(nargin<4)
        IF_SYMB = false ;
    end

    Jparam = sprintf(['../%s/Parameters/%dpop/%s/Jparam.txt'], model,nbpop,dir) ;
    Jdata = importdata(Jparam,' ') ;

    if(IF_SYMB)
        C = Create_Cab(Jdata) ;
        J = sym('J%d%d',nbpop).*C ;       
    else
        J = Jdata ;
    end
end