clear all ;
GlobalVars

Iext = .3 * ExternalInput(model,nbpop,dir) ;    
J = 300 * ImportJab(model,nbpop,dir) ;

Rates = linsolve(J,-Iext.') 

Tm = [20 10 10 10] ; 
Tsyn = ImportTsyn(model,nbpop,dir) ;

for i=1:nbpop
    for j=1:nbpop
        PSP(i,j) = Tm(i) ./ ( Tm(i)-Tsyn(i,j) ) * J(i,j) ;
    end
end
PSP