function [Cond C_Matrix] = BalCond(model,nbpop,dir,Iext,J,IF_DISPLAY)
           
    warning off ;

    % only want 3 optional inputs at most 
    if(nargin<4 || isempty(Iext) ) 
        Jparam = sprintf(['../%s/Parameters/%dpop/%s/Jparam.txt'],model,nbpop,dir);
        Jdata = importdata(Jparam) ;
        Iext = ExternalInput(model,nbpop,dir) ;
        J = Jdata ;
    end

    if(nargin<6)
        IF_DISPLAY=1
    end

    IF_SYM = 1 ;
    if(IF_SYM)
        % if(nbpop>=3)
        %     Iext(3) = 0 ;
        % end
                 
        C = CreateCab(Jdata) ;
        J = sym('J%d%d',nbpop).*C ;
        J(:,1)= ones([1 nbpop]) ;
        
        assume(J(:,2:nbpop)<=0) ;
        
        Iext = sym('I',[1 nbpop]) ;
        assume(Iext>=0) ;
        
        Rates = sym('m%d',[1 nbpop]) ;
        assume(Rates>=0) ;
        
        Isym = sym('I%d',[1 nbpop]) ; 
        Jsym = sym('J%d%d',nbpop) ;
    
        assume(Isym>0) ;  
        assume(Jsym>0) ;
        % assume(Jsym(:,2)<0) ;
        % assume(Jsym(:,3)>0) ;
        % assume(Jsym(:,4)<0) ;
        
        J = Jsym.*C ;
        J(:,2:nbpop) = -J(:,2:nbpop) ;
        Iext = Isym ;
        % Iext(3) = 0 ;
    end

    Conditions = (de2bi(1:2^(nbpop)-2 )) ;
    Solutions = cell(1,length(Conditions) ) ;
        
    for i=1:length(Conditions)
        
        if(IF_DISPLAY)
            fprintf('(%d) Condition, ', i)
            fprintf('%d ', Conditions(i,:) ) 
            fprintf(':\n')
        end 

        Rates = sym('m%d',[1 nbpop]) ;
        assume(Rates>=0) ;
        
        Isym = sym('I%d',[1 nbpop]) ;
        
        Rates = (Rates.*Conditions(i,:)).' ;
        Inputs = Iext.' + J*Rates ;
        
        nonInputs = Inputs( find( Conditions(i,:).'==0 ) ) ;
        Inputs = Inputs( find( Conditions(i,:).' ) ) ;

        Rates = Rates(find(Conditions(i,:).')) ;
        Isym = Isym(find(Conditions(i,:).')) ;

        %% Input always negative
        if( any( isAlways( Inputs<=0 ) ) )
            
            if(IF_DISPLAY)
                fprintf('Absurd, ')
                fprintf('%s silent, ', Rates( find(Inputs<=0) ) )
                fprintf('since ')
                fprintf('%s is always <=0, ', Isym( find(Inputs<=0) ) )
                fprintf('\n')
            end
            Cond(i) = 1 ;
            
            %% Input always positif 
        elseif( any( isAlways( Inputs>=0 ) ) )
                        
            Solution =  Inputs( find( isAlways(Inputs>=0) ) );
            syms m1 m2 m3 m4 ;

            m1 = 0 ; 
            m2 = 0 ;
            m3 = 0 ;
            m4 = 0 ;
            
            %% Saturation due to Iext = O(sqrtK)
            if( subs(Solution)~=0 )   
                
                if(IF_DISPLAY)
                    fprintf('Absurd, ')
                    fprintf('%s saturates, ', Rates( find(Solution~=0) ) )
                    fprintf('since ')                
                    fprintf('%s is always >0, ', Isym( find(Solution~=0) ) )
                    fprintf('\n')
                end
                Cond(i) = 1 ;

                %% mk = O(1/sqrtK) = eps/sqrtK
            else 
               
                if(IF_DISPLAY)
                    fprintf('mk = O(1/sqrtK) = eps/sqrtK\n')
                end
                syms m1 m2 m3 m4 positive;
                m1 = 0 ; 
                
                Inputs = subs(Inputs(find(~isAlways(subs(Inputs)<=0))) ) ;
                
                m =  sym('m%d',[1 nbpop]) ;
                m1 = 0 ; 
                m = nonzeros( subs(m( find( Conditions(i,:) ) ) ) ) ;
                
                [A, b] = equationsToMatrix( Inputs==0, m) ;
                
                Rates = linsolve(A,b) ;
                
                %% eps>0
                if( ~any( isAlways(Rates<=0) ) ) 
                    if(IF_DISPLAY)
                        fprintf('eps>0\n')
                    end
                    %% Input always positif

                    % Cond(i) = 0 ; % TO DELETE
                    
                    if( any( isAlways( nonInputs>=0 ) ) )
                        Solution = nonInputs( find( isAlways( nonInputs>=0 ) ) ) ;

                        if(IF_DISPLAY)
                            fprintf('Solution is always true, ')
                            fprintf('%s is always >0, ', Isym(find( isAlways( nonInputs>=0 ) ) ) )
                            fprintf('\n')
                        end
                        Cond(i) = 1 ;

                    else                        
                        nonInputs = subs(nonInputs) ;
                        nonInputs = subs(nonInputs, m, Rates ) ; 
                        Solution = [-Rates ; nonInputs] ;
                        if(IF_DISPLAY)
                            fprintf('Balance conditions :\n')
                            fprintf('%s >0 \n ', Solution)
                        end
                        % Solutions{i} = [Solutions{i} ; -Rates ; nonInputs] ; 
                        Solutions{i} = [-Rates ; nonInputs] ; 
                        % Cond(i) = any( isAlways( Solutions{i}>0 ) ) ; 
                        Cond(i) = any( isAlways( Solutions{i}>0 ) ) ; 
                    end
                    
                    if(IF_DISPLAY)
                        fprintf('\n')
                    end
                    %% eps<=0
                else
                    if(IF_DISPLAY)                        
                        fprintf('Solution is always true, eps<0 \n')
                    end
                    Cond(i) = 1 ;
                end
                
            end
            
            %% General Case
        else
            
            syms m1 m2 m3 m4 positive;           
            m = sym('m%d',[1 nbpop]) ;

            m = m( find( Conditions(i,:) ) ) ;
            [A, b] = equationsToMatrix( Inputs == 0, m) ;

            Rates = linsolve(A,b).' ;
            
            %% Rate always negatif
            if( any( isAlways( nonzeros(Rates)<=0 ) ) )
                
                if(IF_DISPLAY)
                    fprintf('Solution is always true, ')
                    fprintf('%s is always <0, ', m( find(isAlways(nonzeros(Rates)<=0 )) ) )
                    fprintf('\n')
                end
                Cond(i) = 1 ;
            else
                    
                nonInputs = nonzeros( subs(nonInputs, m, Rates) ) ;
                
                %% Input always positif
                if( any( isAlways( nonInputs>=0 ) ) )                         
                    if(IF_DISPLAY)
                        fprintf('Solution is always true, %s', nonInputs(find(isAlways( nonInputs>=0 ))) )
                        fprintf('\n')             
                    end
                    Cond(i) = 1 ;
                else
                    R = -nonzeros(Rates) ;
                    Solution = [R ; nonInputs] ;

                    if(IF_DISPLAY)
                        fprintf('Balance conditions :\n')
                        fprintf('%s >0 \n', Solution)
                        fprintf('\n')
                    end

                    if( ~any(abs(Solution) == Inf ) )
                        % Solutions{i} = [Solutions{i} ; Solution] ;
                        Solutions{i} = [Solution] ;
                        Cond(i) = any( isAlways( Solutions{i}>0 ) ) ; 
                    else
                        Cond(i) = 1 ; 
                    end
                end                   
            end
        end

    end
    
    fprintf('\n')
    Solutions = Solutions(~cellfun('isempty',Solutions)) ;
    
    for i=1:1:length(Solutions)        
        [Num Den] = numden(Solutions{i}) ;
        Solutions{i} = simplify(Num.*Den) ;
    end

    Solutions = Solutions.' ;    
    
    C_Matrix = [Solutions{:}] ;
    C_Matrix(isnan(C_Matrix))=0 ;

    % if(nbpop>3)                         % 
    %     fprintf('add Cond \n')
    %     J2pop = J(1:2,1:2) ;
    %     I2pop = Iext(1:2) ;
    %     Rates2pop = linsolve(J2pop,-I2pop.') ;
    
    %     Rates1pop = -Iext(2) / J(2,2) ;

    %     fprintf('Cond 1 1 0 sqrtK: ')
    %     C1 = Iext(4) + J(4,1) * Rates2pop(1) + J(4,2) * Rates2pop(2) <0 || any(Rates2pop<0) ;
    %     fprintf('%d \n',C1) 
    
    %     fprintf('Cond 0 1 0 sqrtK: ')
    %     C2 = Iext(4) + J(4,2) * Rates2pop(2) <0  ||  Iext(1) + J(1,2) * Rates1pop >0  || Rates1pop<0 ;
    %     fprintf('%d \n',C2) 
    %     Cond = [ Cond, C1 , C2 ] ; 
    % end

end