function C_Matrix = BalConditions(model,nbpop,dir,Iext,J,IF_SYM,IF_RETURN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function C_Matrix = Bal_Cond(nbpop,dir,Iext,J,IF_SYM,IF_RETURN)
% Check the balanced conditions for a given set of Jij
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    warning off ;

    if(nargin<4 | isempty(Iext))
        Iext = ExternalInput(model,nbpop,dir) ;
    end
    
    if(nargin<5 | isempty(J)) 
        J = ImportJab(model,nbpop,dir) ;
    end
    
    C = CreateCab(J) ;
    
    if(nargin<6 | IF_SYM)
        IF_SYM = 1 ;
        
        Isym = sym('I%d',[1 nbpop]) ; 
        Jsym = sym('J%d%d',nbpop) ;
        
        assume(Isym>0) ;
        assume(Jsym(:,1)>0) ;
        assume(Jsym(:,2:end)<0) ;

        J = Jsym.*C ;
        Iext = Isym ;
    end
    
    if(nargin<7 | IF_RETURN) 
        IF_RETURN = 1 ;
    end

    if(nbpop>=3) % No external input onto SOM
        Iext(3) = 0 ;
    end
    
    Conditions = sortrows(de2bi(1:2^(nbpop)-2 )) ; 
    Solutions = cell(1,length(Conditions) ) ; 
    Cond = zeros(nbpop,1) ; 
    
    for i=1:length(Conditions) 
        Rsym = sym('m%d',[1 nbpop]) ; 
        assume(Rsym>0) 

        syms m1 m2 m3 m4 ; 
        assume(m1>0) ; 
        assume(m2>0) ; 
        assume(m3>0) ; 
        assume(m4>0) ; 
        
        Cond = Conditions(i,:) ; 

        if(IF_RETURN) 
            fprintf('(%d) Condition : ', i) 
            fprintf('%d ', Cond ) 
            fprintf('| ')
        end

        RatesCond = (Rsym.*Cond).' ;
        InputCond = Iext.' + J*RatesCond ;
        
        IdxBal = find(Cond~=0) ;
        InputBal = InputCond(IdxBal) ;
        IdxNeg = find(Cond==0) ;
        InputNeg = InputCond(IdxNeg) ;

        IdxIextZ = find(isAlways(Iext==0)) ;
        
        if(nbpop>=3)
            % me=0 & mi=0 & ms = O(1) & Input SOM == 0 => me = O(sqrt(K))
            if( any(isAlways(InputBal>=0) ) & Cond(1)==0 & Cond(2)==0 & Cond(3)==1 )
                
                IdxHas = find(~has( InputCond(IdxIextZ), [m2 m3 m4])) ;

                if(Cond(IdxHas)>=0)
                    
                    if(IF_RETURN) 
                        fprintf('me is order 1/sqrt(K) \n')
                    end

                    InputBalTmp = subs(InputBal, m1, 0) ;
                    InputNegTmp = subs(InputNeg, m1, 0) ;

                    InputBalTmp(1) = InputNegTmp(1) ;
                    InputNegTmp =  nonzeros(InputNegTmp) ;

                    if( any( isAlways(InputNegTmp>=0) ) )
                        IdxTmp = find(isAlways(InputNegTmp>=0)) ;

                        if(IF_RETURN) 
                            fprintf('Always True : non-balanced input %s always positive \n' , InputNegTmp(IdxTmp) )
                        end

                    else
                        % fprintf('InputBal \n')
                        % fprintf('%s \n',InputBalTmp)
                        [A, b] = equationsToMatrix( InputBalTmp==0, Rsym) ;
                        Rsym = linsolve(A,b) ;
                        % fprintf('Rsym ')
                        % fprintf('%s ',Rsym)
                        % fprintf('\n')
                        m1 = Rsym(1) ;
                        m2 = Rsym(2) ;
                        if(nbpop>=3)
                            m3 = Rsym(3) ;
                            if(nbpop>=4)
                                m4 = Rsym(4) ;
                            end
                        end
                        
                        Rsym = nonzeros(Rsym) ;
                        
                        if( any( isAlways(Rsym<=0) ) )
                            IdxTmp = find( isAlways(Rsym<=0) ) ;

                            if(IF_RETURN)
                                fprintf('But always False : balanced Rates %s always negative \n', Rsym(IdxTmp) )
                            end

                        else
                            if( any( isAlways(Rsym>=0) ) )
                                IdxTmp = find( ~isAlways(InputNegTmp<=0) ) ;
                                % fprintf('InputNeg \n')
                                % fprintf('%s \n',InputNegTmp)
                                Solution = subs(InputNegTmp(IdxTmp)) ; 
                            else
                                Solution = [ -Rsym subs(InputNegTmp) ] ;
                            end
                            [m,n] = size(Solution) ;
                            Solution = nonzeros(reshape(Solution,m*n,1)) ;

                            if(IF_RETURN)
                                fprintf('%s >0 | ', Solution)
                                fprintf('\n')
                            end

                            if(IF_SYM)
                                Solutions{i} = [Solutions{i} ; Solution ] ; 
                            else
                                Solutions{i} = [Solutions{i} ; any(Solution>0) ] ; 
                            end

                        end
                    end

                else
                    IdxTmp = find(isAlways(InputNeg>=0)) ;

                    if(IF_RETURN)     
                        fprintf('Always True : non-balanced input %s always positive \n' , InputNeg(IdxTmp) )
                    end

                end

            elseif( any( isAlways(InputBal<=0) ) )
                IdxTmp = find(isAlways(InputBal<=0)) ; 
                
                if(IF_RETURN)
                    fprintf('Always True : balanced input %s always negative \n', InputBal(IdxTmp) )
                end

            else

                if(IF_RETURN)     
                    fprintf('Non balanced solution exists unless :\n')
                end
                
                % fprintf('InputBal \n')
                % fprintf('%s \n',InputBal)
                [A, b] = equationsToMatrix( InputBal==0, Rsym) ;

                Rsym = linsolve(A,b) ; 
                % fprintf('Rsym ')
                % fprintf('%.2f ',double(Rsym))
                % fprintf('\n')

                m1 = Rsym(1) ;
                m2 = Rsym(2) ;
                if(nbpop>=3)
                    m3 = Rsym(3) ;
                    if(nbpop>=4)
                        m4 = Rsym(4) ;
                    end
                end

                Rsym = nonzeros(Rsym) ;

                if( any( isAlways(Rsym<=0) ) )
                    IdxTmp = find( isAlways(Rsym<=0) ) ;
                    if(IF_RETURN)
                        fprintf('But always False : balanced Rates %s always negative \n', Rsym(IdxTmp) )
                    end

                elseif( any(Rsym==nan) | any(Rsym==Inf))  
                    if(IF_RETURN) 
                        fprintf('But always False : non-balanced state requires fine-tuning \n')
                    end
                else
                    if( any( isAlways(Rsym>=0) ) )
                        IdxTmp = find( ~isAlways(InputNeg<=0) ) ;
                        % fprintf('InputNeg \n')
                        % fprintf('%s \n',InputNeg)
                        
                        Solution = subs(InputNeg(IdxTmp)) ;
                    else
                        Solution = [ -Rsym subs(InputNeg) ] ;
                    end
                    [m,n] = size(Solution) ;
                    Solution = nonzeros(reshape(Solution,m*n,1)) ;
                    
                    if(IF_RETURN)     
                        fprintf('%s >0 | ', Solution)
                        fprintf('\n')
                    end

                    if(IF_SYM)
                        Solutions{i} = [Solutions{i} ; Solution ] ; 
                    else
                        Solutions{i} = [Solutions{i} ; any(Solution>0) ] ; 
                    end
                    
                end
            end
        end
        C_Matrix = [Solutions{:}].' ;
        C_Matrix(isnan(C_Matrix))=0 ;
        
    end