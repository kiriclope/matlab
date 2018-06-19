function [] = BalGenParam(model,nbpop,dir,IF_WRITE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = Bal_GenParam(model,nbpop,dir,IF_WRITE)
% Generates balance set of parameters [Iext J] and if IF_WRITE 
% writes J and steady-state infos to file Jparam.txt and Param.txt
% Utils : ImportJab(model,nbpop,dir)
%         CreateCab(J) 
%         BalRatesMF()
%         BalConditions()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    warning off ;
    if(nargin<4)
        IF_WRITE = false ;
    end

    rng('shuffle') ;
    
    J = ImportJab(model,nbpop,dir) ;
    C = CreateCab(J) ;
    
    ntrial = 0 ;
    lbdSign = 0 ;

    while( lbdSign==0 )

        ntrial = ntrial + 1 ;

        % Initialize variables
        Rates = zeros(1,nbpop) ;
        Det = 0 ;
        Cond = zeros(1,nbpop) ;
        Slopes = zeros(1,nbpop) ;
        Relbd = zeros(1,nbpop) ;
        Imlbd = zeros(1,nbpop) ;
        
        % Create from random array Iext and Jab
        rd_array = round(rand(nbpop+1,nbpop)+.1,3) ;
        Iext = rd_array(1,:) ;
        
        if(nbpop>=3)
            Iext(3) = 0 ;
        end
        
        rd_array(2:end,2:end) = -rd_array(2:end,2:end) ;
        J = rd_array(2:end,1:end).*C ;
        
        fprintf('Trial %d, Parameters\n',ntrial)
        
        % Compute Rates and Det
        fprintf('Computing MF Rates and Det ...\n')
        try
            [Rates Det] = BalRatesMF(model,nbpop,dir,Iext,J) ;
            fprintf('Rates ')
            fprintf('%.3f ',Rates)
            fprintf('\n')
            fprintf('Det ')
            fprintf('%.3f ',Det)
            fprintf('\n')
        catch
            fprintf('\nError Computing Rates and Det !\n')
        end
        
        % Check Rates and Det signs
        RateSign = all(Rates>0) ;
        DetSign = Det*(-1).^nbpop>0 ;
        
        if( RateSign & DetSign ) 

            fprintf('Iext ')
            fprintf('%.3f ',Iext)
            fprintf('\n')
            
            fprintf('Jab\n')
            for i=1:nbpop
                for j=1:nbpop
                    fprintf('%.3f ',J(i,j))
                end
                fprintf('\n')
            end

            % Compute MF susceptibility
            fprintf('Computing MF Slopes ...\n')
            try
                Slopes = BalPrtrSlopesMF(model,nbpop,dir,Iext,J,true) ;
                fprintf('Slopes ')
                fprintf('%.3f ',Slopes)
                fprintf('\n')
            catch
                fprintf('Error Computing Slopes\n')
                Slopes = NaN(1,nbpop) ;
            end
            
            if( ~any(isnan(Slopes)) )
                fprintf('Checking Balance Conditions ...')

                if(nbpop==2)
                    if(Rates(2)>Rates(1) & Iext(1)>Iext(2))
                        try
                            Cond = BalConditions(model,nbpop,dir,Iext,J,0,0) ;
                        catch
                            fprintf('\nError Checking Balance Conditions !\n')
                        end
                    else
                        Cond = 0 ;
                    end
                end

                if(nbpop==3)
                    if( Rates(2)>Rates(1) )
                        try 
                            Cond = BalConditions(model,nbpop,dir,Iext,J,0,0) ;
                        catch
                            fprintf('\nError Checking Balance Conditions !\n')
                        end
                    else
                        Cond = 0 ;
                    end
                end

                if(nbpop==4)
                    if( Rates(2)>Rates(1) & Slopes(1)<0 & Slopes(2)>=0 )
                        try
                            Cond = BalConditions(model,nbpop,dir,Iext,J,0,0) ;
                        catch
                            fprintf('\nError Checking Balance Conditions \n!')
                        end
                    else
                        Cond = 0 ;
                    end
                end
            end

            if( all(Cond==true) )
                
                fprintf(' true\n')
                fprintf('Checking Steady-State Stability ...\n')
                
                try 
                    [Relbd Imlbd]= BalStab(model,nbpop,dir,double(Iext),nbpop,Inf,1,double(J));
                    fprintf('Relbd ')
                    fprintf('%.3f ',Relbd)
                    fprintf(' Imlbd ')
                    fprintf('%.3f ',Imlbd)
                    fprintf('\n')
                catch
                    fprintf('Error Computing Eigenvalues !\n')
                end
                try
                    lbdSign = all(Relbd<0) ;
                catch
                    lbdSign = 0 ;
                end
                
            else
                fprintf(' false\n')
                lbdSign = 0 ;
            end
            
        else 
            fprintf('Rates or Det have the wrong sign !\n')
        end
    end

    if(IF_WRITE)
        path = sprintf(['../%s/Parameters/%dpop/%s'],model,nbpop,dir) ;
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
        
        fprintf(file,'Det %.3f\n', Det) ;
        fprintf(file,'MF Rates ') ;
        fprintf(file,'%.3f | ', Rates) ;
        fprintf(file,'\n') ;
        fprintf(file,'MF Slopes ') ;
        fprintf(file,'%.3f | ', Slopes) ;
        fprintf(file,'\n') ;
        fprintf(file,'Stability ') ;
        fprintf(file,'Relbd ') ;
        fprintf(file,'%.3f ',Relbd) ;
        fprintf(file,'Imlbd ') ;
        fprintf(file,'%.3f ',Imlbd) ;
        fprintf(file,'\n') ;
        
        fclose(file) ;
    end

end