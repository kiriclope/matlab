function Tsyn = ImportTsyn(model,nbpop,dir)
    
    try
        Tparam = sprintf(['../%s/Parameters/%dpop/%s/Tparam.txt'],model,nbpop,dir);
        Tsyn = importdata(Tparam,' ') ;
    catch
        Tsyn = ones(nbpop) ;
        Tsyn(1:nbpop,1) = 2 ;
        Tsyn(1:nbpop,2:nbpop) = 2 ;
    end
end