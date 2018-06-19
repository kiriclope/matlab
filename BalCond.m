function [] = BalCond(model,nbpop,dir)
   
    warning off ;

    % only want 3 optional inputs at most    
    Jparam = sprintf(['../%s/Parameters/%dpop/%s/Jparam.txt'],model,nbpop,dir);
    Jdata = importdata(Jparam) ;
    
    C = CreateCab(Jdata) ;
    J = sym('J%d%d',nbpop).*C ;
    J(:,1)= ones([1 nbpop]) ;

    assume(J(:,2:nbpop)<=0) ;

    Iext = sym('I',[1 nbpop]) ;
    assume(Iext>=0) ;

    Rates = sym('m%d',[1 nbpop]) ;
    assume(Rates>=0) ;
    
    if(nbpop>=3)
        Iext(3) = 0 ;
    end

    Iext = ExternalInput(model,nbpop,dir) ;
    J = Jdata ;

    Conditions = (de2bi(1:2^(nbpop)-2 )) ;
    Solutions = cell(1,length(Conditions) ) ;
        
    for i=1:length(Conditions)

        pause ;
        fprintf('(%d) Condition, ', i)
        fprintf('%d ', Conditions(i,:) ) 
        fprintf(':\n')
                
        Rates = sym('m%d',[1 nbpop]) ;
        assume(Rates>=0) ;
        
        Isym = sym('I%d',[1 nbpop]) ;
        
        Rates = (Rates.*Conditions(i,:)).' ;
        Inputs = Iext.' + J*Rates ;
        
        nonInputs = Inputs( find( Conditions(i,:).'==0 ) ) 
        Inputs = Inputs( find( Conditions(i,:).' ) ) 

        Rates = Rates(find(Conditions(i,:).')) ;
        Isym = Isym(find(Conditions(i,:).')) ;

        % Input always negative
        if( any( isAlways( Inputs<=0 ) ) )
            
            fprintf('Absurd, ')
            fprintf('%s silent, ', Rates( find(Inputs<=0) ) )
            fprintf('since ')
            fprintf('%s is always <=0, ', Isym( find(Inputs<=0) ) )
            fprintf('\n')
            
        % Input always positif
        elseif( any( isAlways( Inputs>=0 ) ) )
                        
            Solution =  Inputs( find( isAlways(Inputs>=0) ) );
            syms m1 m2 m3 m4 ;

            m1 = 0 ; 
            m2 = 0 ;
            m3 = 0 ;
            m4 = 0 ;
            
            % Saturation due to Iext = O(sqrtK)
            if( subs(Solution)~=0 )   
                
                fprintf('Absurd, ')
                fprintf('%s saturates, ', Rates( find(Solution~=0) ) )
                fprintf('since ')                
                fprintf('%s is always >0, ', Isym( find(Solution~=0) ) )
                fprintf('\n')

            % mk = O(1/sqrtK) = eps/sqrtK
            else 
               
                fprintf('mk = O(1/sqrtK) = eps/sqrtK\n')
                syms m1 m2 m3 m4 positive;
                m1 = 0 ; 
                
                Inputs = subs(Inputs(find(~isAlways(subs(Inputs)<=0))) ) 
                
                m =  sym('m%d',[1 nbpop]) ;
                m1 = 0 ; 
                m = nonzeros( subs(m( find( Conditions(i,:) ) ) ) ) ;
                
                [A, b] = equationsToMatrix( Inputs==0, m) ;
                
                Rates = linsolve(A,b) 
                
                % eps>0
                if( ~any( isAlways(Rates<=0) ) )
                    fprintf('eps>0\n')
                    % Input always positif

                    if( any( isAlways( nonInputs>=0 ) ) )
                        Solution = nonInputs( find( isAlways( nonInputs>=0 ) ) ) ;

                        fprintf('Solution is always true, ')
                        fprintf('%s is always >0, ', Isym(find( isAlways( nonInputs>=0 ) ) ) )
                        fprintf('\n')
                        
                    else                        
                        nonInputs = subs(nonInputs) ;
                        nonInputs = subs(nonInputs, m, Rates ) ;
                        Solution = [-Rates ; nonInputs] ;
                        fprintf('Balance conditions :\n')
                        fprintf('%s >0 \n ', Solution)
                        Solutions{i} = [Solutions{i} ; -Rates ; nonInputs] ;                        
                    end

                    fprintf('\n')

                % eps<=0
                else
                    fprintf('Solution is always true, eps<0 \n')
                end
                
            end
            
        % General Case
        else
            
            syms m1 m2 m3 m4 positive;           
            m = sym('m%d',[1 nbpop]) ;

            m = m( find( Conditions(i,:) ) ) ;
            [A, b] = equationsToMatrix( Inputs == 0, m) ;

            Rates = linsolve(A,b).' ;
            
            % Rate always negatif
            if( any( isAlways( nonzeros(Rates)<=0 ) ) )
                
                fprintf('Solution is always true, ')
                fprintf('%s is always <0, ', m( find(isAlways(nonzeros(Rates)<=0 )) ) )
                fprintf('\n')
                
            else
                    
                nonInputs = nonzeros( subs(nonInputs, m, Rates) )
                
                % Input always positif
                if( any( isAlways( nonInputs>=0 ) ) )                         
                    fprintf('Solution is always true, %s', nonInputs(find(isAlways( nonInputs>=0 ))) )
                    fprintf('\n')             
                else
                    R = -nonzeros(Rates) ;
                    Solution = [R ; nonInputs] ;

                    fprintf('Balance conditions :\n')
                    fprintf('%s >0 \n', Solution)
                    fprintf('\n')
                    
                    Solutions{i} = [Solutions{i} ; Solution] ;
                        
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
                 
end