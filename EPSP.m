clear all ;
GlobalVars 

Iext = 20 * ExternalInput(model,nbpop,dir) *.01 ; 
J = g*20 * ImportJab(model,nbpop,dir) ; 

fprintf('det J %.2f \n', det(J) /20 ) ; 
Tm = [20 10 20 20] ; 
Tsyn = ImportTsyn(model,nbpop,dir) ; 

fprintf('Iext ')
fprintf('%.2f ',round(Iext,2)*100 ) 
fprintf('\n')

fprintf('J : \n')
for i=1:nbpop
    for j=1:nbpop
        fprintf('%.2f ', round( J(i,j), 2) ) 
    end
    fprintf('\n')
end
fprintf('\n')

[Rates DetJ]= BalRatesMF(model,nbpop,dir,Iext,J,0) ; 
Rates = Rates * 1000 ; 

fprintf('MF Rates : ') 
fprintf('%.2f ', round(Rates,2) ) 
fprintf('\n') 

[u b] = RateInputDist(model,nbpop,dir,Iext/20,K,1,J/20,0) ; 

fprintf('Finite K Rates : ') 
fprintf('%.2f ', round(QchAvgTF(u,b) *1000,2) )
fprintf('\n') 

if(nbpop==4)

    % Slopes(1) = ( J(1,2) * J(4,3) - J(1,3) * J(4,2) )  / det(J) ; ...
    % % chiEI 
    % Slopes(2) = - ( J(1,1) * J(4,3) - J(1,3) * J(4,1) ) / det(J)  ; % chiII 
    % Slopes(3) = -( J(1,2) * J(4,1) - J(1,1) * J(4,2) ) /det(J) ;
    % %chiSI 
    Slopes(1) = ( J(1,2) * J(4,3) - J(1,3) * J(4,2) ) / Rates(1) / det(J) ; ...
    % chiEI 
    Slopes(2) = - ( J(1,1) * J(4,3) - J(1,3) * J(4,1) ) / Rates(2) / det(J)  ; % chiII 
    Slopes(3) = -( J(1,2) * J(4,1) - J(1,1) * J(4,2) ) / Rates(3) / det(J) ; ...
    % chiSI 
    
    Slopes(4) = ( J(2,3) * J(4,2) - J(2,2) * J(4,3) ) / Rates(1) / det(J); % chiEE 
    Slopes(5) = - ( J(2,3) * J(4,1) - J(2,1) * J(4,3) ) / Rates(2) / det(J) ; ...
    % chiIE 
    Slopes(6) = - ( J(2,1) * J(4,2) - J(2,2) * J(4,1) ) / Rates(3) / det(J) ; ...
    % chiSE 
    
    fprintf('Norm. Slopes I: ') 
    fprintf('%.2f ', round(Slopes(1:3) * 1000,2) )
    fprintf('\n') 
    
    fprintf('Norm. Slopes E: ') 
    fprintf('%.2f ', round(Slopes(4:6) * 1000,2) )
    fprintf('\n') 

end
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
