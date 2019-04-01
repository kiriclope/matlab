clear all ;
GlobalVars 

Iext = 20 * ExternalInput(model,nbpop,dir) *.01 ; 
J = 200 * ImportJab(model,nbpop,dir) ; 

fprintf('det J %.2f \n', det(J) /20 ) ; 
Tm = [20 10 20 20] ; 
Tsyn = ImportTsyn(model,nbpop,dir) ; 

% Iext = 30*.01 ; 
% J = -30 ; 
% Tm = 20 ; 
% Tsyn = 2 ; 
% nbpop = 1 ; 

fprintf('Iext ')
fprintf('%.2f ',round(Iext,2)) 
fprintf('\n')

fprintf('J : \n')
for i=1:nbpop
    for j=1:nbpop
        fprintf('%.2f ', round( J(i,j), 2) ) 
    end
    fprintf('\n')
end
fprintf('\n')

Rates = linsolve(J,-Iext.') * 1000. ; 

[Rates DetJ]= BalRatesMF(model,nbpop,dir,Iext,J,0) ; 
Rates = Rates * 1000 ; 

fprintf('MF Rates : ') 
fprintf('%.2f ', round(Rates,2) ) 
fprintf('\n') 

[u b] = RateInputDist(model,nbpop,dir,Iext./20./.01,K,1,J/200,0) ; 

fprintf('Finite K Rates : ') 
fprintf('%.3f ', QchAvgTF(u,b) ) 
fprintf('\n') 

% Slopes(1) = J(2,3) * J(4,2) - J(2,2) * J(4,3) ; % chiEE
% Slopes(2) = J(1,2) * J(4,3) - J(1,3) * J(4,2) ; % chiEI

% Slopes(3) = - J(2,3) * J(4,1) + J(2,1) * J(4,3) ; chiIE
% Slopes(4) = - J(1,1) * J(4,3) + J(1,3) * J(4,1) ; chiII

% Slopes(5) = - J(1,2) * J(4,1) + J(1,1) * J(4,2) ; chiSI

% fprintf('Slopes : ') 
% fprintf('%.3f ', Slopes ) 
% fprintf('\n') 

% fprintf('Norm Slopes : ') 
% fprintf('%.3f %.3f', Slopes(2)./Rates(1), Slopes(4)./Rates(2) ) 
% fprintf('\n') 

fprintf('PSPs : \n')

for i=1:nbpop
    for j=1:nbpop
        TMAX = Tm(i) * Tsyn(i,j) * log( Tm(i)/Tsyn(i,j) ) ./ ( Tm(i)-Tsyn(i,j) ) ;
        MAX = exp(-TMAX/Tm(i)) - exp(-TMAX/Tsyn(i,j) ) ;
        PSP(i,j) = Tm(i) ./ ( Tm(i)-Tsyn(i,j) ) * J(i,j) * MAX ./ sqrt(K) ; 
        fprintf('%.2f ', round(PSP(i,j),2 ))        
    end
    fprintf('\n')
end
