function Cpt = CptNeuron(nbpop,nbN)
    Cpt = zeros(nbpop+1,1) ;
   
    for i=1:(nbpop+1)
        for j=1:i-1
            Cpt(i) =  Cpt(i) + nbN(j) ;
        end
    end

    % fprintf('Cpt ')
    % fprintf('%d | ',Cpt)
    % fprintf('\n')

end