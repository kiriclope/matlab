function [r J] = STP_Solver(nbpop,dir,IF_WRITE)

    if(nargin<3)
        IF_WRITE=0 ;
    end

    switch nbpop

      case 1
        Trec = .20E3 ;
        Tfac = 1.5E3 ;
        U = 0.3 ;

      case 2
        
        Iext = [ 1 .5 ] ;
        Trec = [.2 0; 0 0] ;
        Tfac = [.6 0; 0 0] ;
        U = [.05 1; 1 1] ; 

        % Trec = [ .1 .1 ;.8 .8 ] ;
        % Tfac = [ .5 .5 ;.003 .003 ] ;
        % U = [ .03 .03 ; .5 .5 ] ;
          
      case 3

        Iext = [ 1 .5 0 ] ;
        Trec = [ 0 0 0 ; 0 0 0 ; .2 .2 0 ]  ;
        Tfac = [ 0 0 0 ; 0 0 0 ; .6 .6 0] ;
           U = [ 1 1 1 ; 1 1 1 ; .2 .05 1 ] ;
                
      case 4

        Iext = [ 1 .5 0 1] ;
        Trec = [ 0 0 0 0 ; 0 0 0 0 ; .1 0 0 0 ; 0 0 0 0] ;
        Tfac = [ 0 0 0 0 ; 0 0 0 0 ; .5 0 0 0 ; 0 0 0 0] ;
        U = [ 1 1 1 1 ; 1 1 1 1 ; .03 1 1 1 ; 1 1 1 1] ; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    J = ImportJab('STP',nbpop,'Test',0) ; 
    C = CreateCab(J) ; 

    r0 = rand(1,nbpop) ; 
    r = -r0 ;
    options = optimset('Display','off') ; 
    while(any(r<0))
        Iext = .3 .*rand(1,nbpop) ;
        r0 = rand(1,nbpop)+.1 ;
        
        J = rand(nbpop,nbpop)+.1 ;
        J = double( (rand(nbpop,nbpop)+.1).*C ) ;
        if(nbpop>=2) 
            J(:,2:nbpop) = -J(:,2:nbpop) ; 
        end 

        if(nbpop>=3) 
            %J(nbpop,nbpop) = 0 ; 
            Iext(3) = 0 ; 
        end                             
        
        r = fsolve(@SelfCstEq,r0,options) ; 

        for i=1:nbpop
            fprintf('%.3f | ', r(i) ) ; 
        end            
        fprintf('\r') ; 

        if(nbpop>=2) 
            if( r(1)>r(2) || r(1)<1 ) 
                r(1) = -10 ; 
            end 
        end 
    end 

    fprintf('\n') ; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('Iext ') 
    for i=1:nbpop 
        fprintf('%.3f ',Iext(i)) 
    end 
    fprintf('\n') 

    fprintf('Synaptic strength J\n') 
    for i=1:nbpop 
        for j=1:nbpop 
            fprintf('%.3f ',J(i,j)) 
        end 
        fprintf('\n') 
    end 

    fprintf('Rates\n') 
    for i=1:nbpop 
        fprintf('%.3f ',r(i)) 
    end 
    fprintf('\n') 

    if(IF_WRITE) 

        path = sprintf(['../STP/Parameters/%dpop/%s'],nbpop,dir) ; 
        try 
            mkdir(path) ; 
        end 

        fprintf('Writing Parameters to : ') ; 
        filename = sprintf(['%s/Param.txt'], path) ; 
        Jparam = sprintf(['%s/Jparam.txt'], path) ; 
        disp(filename) 
        file = fopen(filename,'w') ; 
        Jfile = fopen(Jparam,'w') ; 
        
        fprintf(file, 'Iext ') ; 
        fprintf(file, '%.3f ', Iext) ; 
        fprintf(file, '\n') ; 
        
        fprintf(file, 'Connectivity Matrix \n') ; 
        for i=1:nbpop 
            for j=1:nbpop 
                fprintf(file, '%.3f ', J(i,j)) ; 
                fprintf(Jfile, '%.3f ', J(i,j)) ; 
            end 
            fprintf(file, '\n') ; 
            fprintf(Jfile, '\n') ; 
        end 
        fclose(Jfile) ; 
        
        fprintf(file,'MF Rates ') ; 
        fprintf(file,'%.3f | ', r) ; 
        fprintf(file,'\n') ; 
        
        fclose(file) ; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Eq = SelfCstEq(r)
        Eq = [] ;
        for i=1:nbpop
            Eq_i =  Iext(i) ;
            for j=1:nbpop                
                u_ij = U(i,j) * ( 1 + Tfac(i,j) * r(j) ) / ( 1 + U(i,j) * Tfac(i,j) * r(j) ) ; 
                x_ij = 1 / ( 1 + Trec(i,j) * u_ij * r(j) ) ; 
                
                Eq_i = Eq_i + J(i,j) * u_ij * x_ij * r(j) ;
            end
            Eq = [Eq ; Eq_i] ;
        end
    end

end