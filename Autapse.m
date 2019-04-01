clear all ;
options = optimset('Display','off') ;

I = 1 ; 
J = 2 ;
K = 50 ;
Pb_Aut = .5 ;

m = [] ;
for Pb_Aut=0:.1:1 
    f = @(m) m-erfc( - sqrt(K) * (I/J-m) / sqrt(m) + ( (I/J-m) / sqrt(m) /2 + sqrt(m) ) * Pb_Aut ) ; 
    flag=0 ;
    while flag<=0
        x0 = rand() ;
        [x,fval,flag] = fsolve(f,x0,options) ;
    end
    m = [m x] ;
    fprintf('Pb_Aut %.3f m %.3f\n', Pb_Aut, m) ;
end

figure(1) ; hold on ;
plot([0:.1:1], m)
xlabel('\kappa')
ylabel('m')

function [m_E m_I] = RatesEqs(m) 
    
    mu_E = sqrt(K) *( I_E + J(1,1) * m(1) + J(1,2) * m(2) ) ;
    mu_I = sqrt(K) *( I_I + J(2,1) * m(1) + J(2,2) * m(2) ) + J(2,2) * Kappa * m(2) ; 

    var_E = J(1,1).^2 * m(1) + J(1,2).^2 * m(2) ;
    var_I = J(2,1).^2 * m(1) + J(2,2).^2 * (1 + Kappa / sqrt(K) ) * m(2) ; 
    
    m_E = erfc( -mu_E/sqrt(var_E) ) ;
    m_I = erfc( -mu_I/sqrt(var_I) ) ;
end

for Pb_Aut=0:.1:1 
    flag=0 ;
    while flag<=0
        x0 = rand(1,2) ;
        [x,fval,flag] = fsolve(@RatesEqs,x0,options) ;
    end
    m = [m x] ;
    fprintf('Kappa %.3f m %.3f \n', Kappa, m) ;
end
