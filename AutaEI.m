function [] = AutaEI()

    warning off ;
    options = optimset('Display','off') ;

    J = [1 -1.5 ; 1 -1] ; 
    I_E = 0.5 ; 
    I_I = 0.25 ; 

    K=100 ;

    m = [] ;    
    Kappa = 0:.05:1 ;
    for i=1:length(Kappa)
        flag=0 ;
        while flag<=0 
            x0 = [rand(1,2) Kappa(i)] ;
            [x,fval,flag] = fsolve(@RatesEqs,x0,options) ;
        end
        m = [m x.'] ;
        fprintf('flag %d Kappa %.3f m_E %.3f m_I %.3f \n', flag, Kappa(i), x(1), x(2) ) ; 
    end
    fprintf('\n')

    figure(1) ; hold on ;
    plot(Kappa,m(1,:),'r')
    plot(Kappa,m(2,:),'b')

    ylabel('m')
    xlabel('\kappa')
    
function Eq = RatesEqs(X)
    mu_E = sqrt(K) *( I_E + J(1,1) * X(1) + J(1,2) * X(2) ) ;
    mu_I = sqrt(K) *( I_I + J(2,1) * X(1) + J(2,2) * X(2) ) + J(2,2) * X(3) * X(2) ;

    var_E = J(1,1).^2 * X(1) + J(1,2).^2 * X(2) ;
    var_I = J(2,1).^2 * X(1) + J(2,2).^2 * (1 + X(3) / sqrt(K) ) * X(2) ; 
    
    if(var_E>0 && var_I>0) 
        Eq = [ X(1) - erfc( -mu_E/sqrt(var_E) ) ; X(2) - erfc( -mu_I/sqrt(var_I) ) ; 0 ] ;
    else 
        Eq = [1 ; 1 ; 1] ;
    end
end

end