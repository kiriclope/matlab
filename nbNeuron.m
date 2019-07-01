function nbN = nbNeuron(nbpop,n,IF_Nk,p)

    nbN = zeros(nbpop,1) ;
    N = n*10000 ;
    %N= n*2500 ;

    for i=1:nbpop 
        nbN(i)=N./nbpop ;
    end

    if IF_Nk
        N = 76800 ;
        if(nbpop==2)
            % nbN(1)= N*75/100 ;
            % nbN(2)= N*25/100 ;

            nbN(1)= 57600 ; 
            nbN(2)= 19200 ;             
        end
        if(nbpop==3)
            nbN(1)= N*70/100 ;
            nbN(2)= N*15/100 ;
            nbN(3)= N*15/100 ;

            if(~isempty(p))
                nbN(1)= N./3 ; 
                nbN(2)= N./3*p ; 
                nbN(3)= N./3*(1-p) ; 
            end
        end
        if(nbpop==4)
            nbN(1)= N*70/100 ;
            nbN(2)= N*10/100 ;
            nbN(3)= N*10/100 ;
            nbN(4)= N*10/100 ; 

            % nbN(1)= 32400 ;             
            % nbN(2)= 3600 ;
            % nbN(3)= 3600 ;
            % nbN(4)= 3600 ;

            nbN(1)= 57600 ; 
            nbN(2)= 6400 ; 
            nbN(3)= 6400 ; 
            nbN(4)= 6400 ; 
        end 
    end

    % fprintf('nbN ')
    % fprintf('%d | ', nbN)
    % fprintf('\n')

end                                     