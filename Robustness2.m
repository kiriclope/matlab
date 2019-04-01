clear all ; 
GlobalVars

RNDMAX = 100 ;

Xmin = 0.9 ; 
Xmax = 1.1 ; 
DX = Xmax - Xmin ;

cpt = 0 ;
meanJee = 0 ;

Tsyn = ImportTsyn(model,nbpop,dir) ; 

Iext0 = ExternalInput(model,nbpop,dir) * 20; 
J0 = ImportJab(model,nbpop,dir) * 20 ; 
Jstar = J0(4,1) * J0(1,3) / J0(4,3) ; 

for i=1:nbpop 
    % rnd = DX * rand()  + Xmin ; 
    % Iext(i) = Iext(i) * rnd ; 

    for j=1:nbpop 

        J = J0 ; 
        Iext = Iext0 ;   
        cpt = 0 ;
        meanJee = 0 ;

        rnd = DX * rand() + Xmin ; 
        J(i,j) = J(i,j) * rnd ; 
        
        for k=1:RNDMAX
            
            Rates = linsolve(J,-Iext.') ; 
            Jee = J(1,1) ;
            Jee_Star = J(4,1) * J(1,3) / J(4,3) ; 
            meanJee = meanJee + Jee_Star ;
    
    % if( Jee-Jee_Star>0 ) 
    %     cpt = cpt + 1 ;
    % end
    
    % fprintf('RND %d Jstar %.3f DetJ %.3f Rates ', k, Jee_Star, det(J)*10 ) 
    % fprintf('%.3f ', Rates )
    % fprintf('\n')

    if( all(Rates>0) ) 
        lbd = eig( J ./ Tsyn ) ; 
        Relbd = real(lbd) ; 
        if( all(Relbd<0) )
            cpt = cpt + 1 ; 
            % dirRND = sprintf('%s_RND_%d',dir,k) ; 
            % WriteParam(model,nbpop,dirRND,Iext,J,Rates) ; 
        end
    end
        end

        fprintf('%d %d Balanced Proportion %.3f\n',i,j, cpt/RNDMAX) 
    end
end
