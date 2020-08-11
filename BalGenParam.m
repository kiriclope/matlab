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

    m0 = .01 ;
    ntrial = 0 ;
    lbdSign = 0 ;
    
    IF_SLOPES = 0 ; 
    SLOPES_SIGNS = [-1 -1 1 -1 1] ; 

    IF_COND = 1 ; 
    IF_STAB = 1 ; 
    IF_FINITE_K = 0 ; 
    IF_INPUTSOM = 1 ; 

    IS_BALANCE = 0 ; 
    while( IS_BALANCE == 0 ) 

        ntrial = ntrial + 1 ;

        %% Initialize variables
        RatesMF = zeros(1,nbpop) ;
        RatesK = zeros(1,nbpop) ;
        Det = 0 ;

        IS_COND = 1 ;

        Slopes = zeros(1,nbpop) ; 
        IS_SLOPES = 1 ; 

        Relbd = zeros(1,nbpop) ;
        Imlbd = zeros(1,nbpop) ; 
        
        %% Create from random array Iext and Jab
        
        % rd_array = round(.25+1.75*rand(nbpop+1,nbpop),2) ; 
        rd_array = round(.25+1.5*rand(nbpop+1,nbpop),2) ; 
        Iext = rd_array(1,:) * m0 ; 


        if(~IF_INPUTSOM && nbpop>2)
            Iext(3) = 0 ; 
            %Iext(4) = 0 ; 
        end

        if(IF_INPUTSOM==-1 && nbpop>2)        
            Iext(3) = -Iext(3) ;        
        end

        rd_array(2:end,2:end) = -rd_array(2:end,2:end) ; % standard balance

        %% Bump EI    EI interacting 
        % rd_array(2:end,2) = -rd_array(2:end,2) ; % EI to EI
        % rd_array(2:end,4) = -rd_array(2:end,4) ; % EI to EI

        J = rd_array(2:end,1:end).*C ;
        
        fprintf('Trial %d, Parameters\n',ntrial)        
        %% Compute Rates and Det
        fprintf('Computing MF Rates and Det ...\n')
        try
            [RatesMF Det] = BalRatesMF(model,nbpop,dir,Iext,J) ; 

            fprintf('Rates MF ')
            fprintf('%.3f ',RatesMF*1000)
            fprintf('\n')

            fprintf('Det ')
            fprintf('%.3f ',Det)
            fprintf('\n')
        catch 
            fprintf('\nError Computing MF Rates and Det !\n') 
        end
        
        %% Check Rates and Det signs
        InhibRates = [ RatesMF(2:nbpop) ];
        
        if(nbpop>=3)
            notPVrates = RatesMF(3:nbpop) ;
        else
            notPVrates = [RatesMF(1)] ; 
        end

        % RateSign = ( all(RatesMF>0) && all( RatesMF(1) < InhibRates ) && all( RatesMF(2) > notPVrates ) ) ; 
        RateSign = ( all(RatesMF>0) && all(RatesMF(2)>notPVrates) && all( InhibRates > RatesMF(1) ) && RatesMF(1) < 5./1000 ) ; 
        % RateSign = ( all(RatesMF>0) && all( RatesMF(1) < RatesMF(2) ...
        %                                     ) && all(RatesMF<20/1000)) ; 
        DetSign = Det*(-1).^nbpop>0 ; 
        
        if( RateSign && DetSign ) 

            fprintf(' TRUE \n')
            IS_BALANCE = 1 ; 
            
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
        else
            IS_BALANCE = 0 ; 
            fprintf(' #### #### ####\n') 
            fprintf(' #### WRONG RATES OR DET SIGNS #### \n') 
            fprintf(' #### #### ####\n') 
        end
        
        if(IF_FINITE_K && IS_BALANCE) 
            try 
                [u b] = RateInputDist(model,nbpop,dir,double(Iext),3000,1,double(J),0) ; 
                RatesK = QchAvgTF(u,b) ;
                fprintf('Rates finite K ')
                fprintf('%.3f ',RatesK)
                fprintf('\n')
            catch
                fprintf('\nError Computing Finite K Rates !\n')
            end
            
            switch nbpop
              case 2
                if( RatesK(2)<RatesK(1) )
                    RateSign = 0 ;
                end
                
              case 3
                if( RatesK(2)<RatesK(1) || RatesK(3)<RatesK(1) )
                    RateSign = 0 ;
                end
                
              case 4
                if( RatesMF(2)<RatesMF(1) || RatesMF(3)<RatesMF(1) ...
                    || RatesMF(4)<RatesMF(1) || RatesMF(3)>RatesMF(2) ...
                    || RatesMF(4)>RatesMF(2) )
                    %if( RatesK(2)<RatesK(1) )
                    RateSign = 0 ;          
                end
            end 
        end
        
        if(IF_SLOPES && IS_BALANCE )
            %% Compute MF susceptibility 
            fprintf('Computing MF Slopes ...\n')
            IS_SLOPES = 0 ;
            % try
            %     Slopes = BalPrtrSlopesMF(model,nbpop,dir,Iext,J,true) ;
            %     fprintf('Slopes ')
            %     fprintf('%.3f ',Slopes)
            %     fprintf('\n')
            % catch
            %     fprintf('Error Computing Slopes\n')
            %     Slopes = NaN(1,nbpop) ;
            % end
            if(nbpop==4)                        

                if strfind(dir,'S1')
                    fprintf('%%%%%%%% S1L5 %%%%%%%%%%\n') 
                    Slopes(1) = J(1,3) * J(4,4) - J(1,4) * J(4,3) ; % chiEI = chiII 
                    Slopes(2) = J(2,4) * J(4,3) - J(2,3) * J(4,4) ; % chiEE = chiIE 
                    Slopes(3) = 1 ; 
                    Slopes(4) = 1 ; 
                    Slopes(5) = 1 ; 
                    
                    %% chiEI = chiII = Jes Jxx - Jex Jxs 
                    %% chiEE = chiIE = Jix Jxs - Jis Jxx 
                elseif strfind(dir,'L23') 
                    fprintf('%%%%%%%% ALML23 %%%%%%%%%%\n') 
                    Slopes(1) = J(1,2) * J(4,3) - J(1,3) * J(4,2) ; % chiEI
                    Slopes(2) = - ( J(1,1) * J(4,3) - J(1,3) * J(4,1) ) ; % chiII
                    Slopes(3) = J(2,3) * J(4,2) - J(2,2) * J(4,3) ; % chiEE
                    Slopes(4) = - ( J(2,3) * J(4,1) - J(2,1) * J(4,3) ) ; % chiIE
                    Slopes(5) = 1 ;
                    
                elseif strfind(dir,'L5')
                    fprintf('%%%%%%%% ALML5 %%%%%%%%%%\n')
                    Slopes(1) = J(1,2) * J(4,3) - J(1,3) * J(4,2) ; % chiEI 
                    Slopes(2) = - ( J(1,1) * J(4,3) - J(1,3) * J(4,1) ) ; % chiII 
                    Slopes(3) = J(2,3) * J(4,2) - J(2,2) * J(4,3) ; % chiEE 
                    Slopes(4) = - ( J(2,3) * J(4,1) - J(2,1) * J(4,3) ) ; % chiIE 
                    Slopes(5) = 1 % -( J(1,2) * J(4,1) - J(1,1) * J(4,2) ) ; % chiSI 
                    
                    %% chiEI = Jei Jvs - Jes Jvi 
                    %% chiII = Jee Jvs - Jes Jve
                    %% chiEE = Jis Jvi - Jii Jvs
                    %% chiIE = Jis Jve - Jie Jvs                    
                    %% chiSI = Jei Jve - Jee Jvi
                else
                    Slopes(1) = J(1,3) * J(3,2) - J(1,2) * J(3,3) ; 
                    Slopes(2) = -( J(1,3) * J(3,1) - J(1,1) * J(3,3) ) ; 
                    Slopes(3) =  -( J(1,1) * J(3,2) - J(1,2) * J(3,1) ) ;
                    Slopes(4) = 1 ;
                    Slopes(5) = 1 ;
                    %% chiEI = Jes Jsi - Jei Jss ;
                    %% chiII = Jes Jse - Jee Jss ;
                    %% chiSI = Jee Jsi - Jei Jse  ;
                end           
            else
                Slopes(1) = -1 ;
                Slopes(2) = -1 ;
                Slopes(3) = 1 ;
                Slopes(4) = 1 ;
                Slopes(5) = 1 ;
            end

            IS_SLOPES = ( ~any(isnan(Slopes)) && all(sign(Slopes)==SLOPES_SIGNS) ) ; 
            
            if(IS_SLOPES) 
                fprintf(' TRUE \n')
                IS_BALANCE = 1 ;
            else
                fprintf(' #### #### ####\n') 
                fprintf(' #### WRONG SLOPES SIGNS #### \n') 
                fprintf(' #### #### ####\n') 
                IS_BALANCE = 0 ;
            end 
        end 
        
        if(IF_COND && IS_BALANCE) 
            fprintf('Checking Balance Conditions ...') 
            IS_COND = 0 ;
            try 
                %Cond = BalConditions(model,nbpop,dir,Iext,J,0,0) ; 
                IS_COND = all( BalCond(model,nbpop,dir,Iext,J,0) == 1 ) ;
            catch
                fprintf('\nError Checking Balance Conditions !\n') 
                IS_COND = 0 ;
            end

            if(IS_COND)
                fprintf(' TRUE \n')            
                IS_BALANCE = 1 ; 
            else 
                fprintf(' #### #### ####\n') 
                fprintf(' #### WRONG BAL COND #### \n') 
                fprintf(' #### #### ####\n') 
                IS_BALANCE = 0 ;
            end
        end

        if( IF_STAB && IS_BALANCE)
            fprintf('Checking Steady-State Stability ...\n')                    
            lbdSign = 0 ;
            try 
                [Relbd Imlbd]= BalStab(model,nbpop,dir,double(Iext),nbpop,3000,1,double(J));
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

            if(lbdSign)
                fprintf(' TRUE \n')
                IS_BALANCE = 1 ;
            else
                fprintf(' #### #### ####\n') 
                fprintf(' #### UNSTABLE #### \n') 
                fprintf(' #### #### ####\n') 
                IS_BALANCE = 0 ;
            end
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
        fprintf(file, '%.3f ', Iext / m0) ; 
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
        fprintf(file,'%.3f | ', RatesMF*1000) ;
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = phi(x)
        out = exp(-x.^2./2)./sqrt(2.*pi);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = Phi(x)
        out = .5.*(1+erf( x./sqrt(2) ) ) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = QchAvgTF(u,b)
        if(b>0)
            out = u.*Phi( u./sqrt(b) ) + sqrt(b).*phi( u./sqrt(b) ) ;
        else
            out = u ;
        end
    end

end