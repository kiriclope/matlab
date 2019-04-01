function NullCline(model,nbpop,dir)
    
    
    J = ImportJab(model,nbpop,dir,0) ; 
    Eqx = @(x) Iext.' + J * x.' ;

    
end