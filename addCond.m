clear all ; 
GlobalVars

Iext = ExternalInput(model,nbpop,dir) ;
J = ImportJab(model,nbpop,dir) ;
Rates = linsolve(J,-Iext.') 

J2pop = J(1:2,1:2) ;
I2pop = Iext(1:2) ;
Rates2pop = linsolve(J2pop,-I2pop.') 

Rates1pop = -Iext(2) / J(2,2) 

fprintf('Cond 1 1 0 sqrtK: ')
C1 = Iext(4) + J(4,1) * Rates2pop(1) - J(4,2) * Rates2pop(2) <0 || any(Rates2pop<0) ;
fprintf('%d \n',C1) 

fprintf('Cond 0 1 0 sqrtK: ')
C2 = Iext(4) - J(4,2) * Rates2pop(2) <0  ||  Iext(1) - J(1,2) * Rates1pop >0  || Rates1pop<0 ;
fprintf('%d \n',C2) 