function Tsyn = ImportTsyn(model,nbpop,dir)
    
    try
        Tparam = sprintf(['../%s/%dpop/Parameters/%s/Tparam.txt'],model,nbpop,dir);
        Tsyn = importdata(Tparam,' ') ;
    catch
        Tsyn = ones(nbpop) ;
        Tsyn(1:nbpop,1) = 3 ;
        Tsyn(1:nbpop,2:nbpop) = 2 ;
    end
end