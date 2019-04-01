function [Rates Det] = BalRatesMF(model,nbpop,dir,Iext,J,DisplayOn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Rates Det] = BalRatesMF(model,nbpop,dir,Iext,J,DisplayOn)
% Computes Balance state rates in the large N, large K limit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if( nargin<6 )
        DisplayOn = 0 ; 
        if(nargin<5 || isempty(J))
            J = ImportJab(model,nbpop,dir) ; 
            if(nargin<4 || isempty(Iext))
                Iext = ExternalInput(model,nbpop,dir) ; 
            end
        end
    end

    Det = det(J) ; 
    Rates = RatesMF_Func(nbpop,Iext,J) ; 
    
    if(DisplayOn)
        fprintf('Rates ') 
        fprintf('%.3f ',Rates) 
        fprintf('\n')
        fprintf('Det ')
        fprintf('%.3f ',Det)
        fprintf('\n')
    end

    function [out] = RatesMF_Func(nbpop,Iext,J)
    
        Rates = cell(1,nbpop) ;
        
        for i=1:nbpop+1
            Rates{i} = zeros(1,nbpop) ;
        end

        if(nbpop>1)
            Rates{1}(2) = linsolve(J(2,2),-Iext(2).') ;
        end

        if(nbpop>2)
            if(nbpop==4)
                if strfind(dir,'S1L5') 
                    Rates{1}(4) = linsolve(J(4,4),-Iext(4).') ;
                    Rates{1}(2) = linsolve(J(2,2), - ( Iext(2) - J(2,4)./ J(4,4) * Iext(4) ).') ; 
                else
                    Rates{1}(4) = 0 ;
                    Rates{1}(2) = linsolve(J(2,2),-Iext(2).') ;
                end
            else
                Rates{1}(2) = linsolve(J(2,2),-Iext(2).') ;
            end

            G = J ;
            I = Iext ;
            
            if strfind(dir,'S1L5') 
                G(3,:) = [] ;
                G(:,3) = [] ;
                I(3) = [] ;
                R = linsolve(G,-I.') ;
                if(all(R>0)) 
                    % fprintf('YOUPI') ; 
                    % fprintf('%.3f ', R) ;
                    % fprintf('\n') ;
                end
                Rates{nbpop+1} = [ R(1) R(2) 0 R(3) ] ;
            else

                if(nbpop==4)
                    G(4,:) = [] ;
                    G(:,4) = [] ;
                    I(4) = [] ;
                end

                G(3,:) = [] ;
                G(:,1) = [] ;
                I(3) = [] ;

                R = linsolve(G,-I.') ; 
                if(all(R>0)) 
                    % fprintf('YOUPI') ; 
                    % fprintf('%.3f ', R) ;
                    % fprintf('\n') ;
                end
                Rates{nbpop+1} = [0 R(1) R(2) 0] ; 
            end
        end

        if(nbpop>1)
            for i=2:nbpop
                Rates{i} = linsolve(J(1:i,1:i),-Iext(1:i).') ; 
                Rates{i}(i+1:nbpop) = 0 ; 
            end
        end

        out = Rates{nbpop} ; 
        nb = nbpop ; 
        while(any(out<=0) & nb>1)
            nb = nb-1 ;
            out = Rates{nb} ;

            if strfind(dir,'S1L5') 
                if(nb~=nbpop)
                    alt = Rates{nbpop+1} ;
                    alt(3) = [] ; 
                    if(all(alt>0)) 
                        % fprintf('YOUPI 2 %d ', all(alt>0) ) ;
                        % fprintf('%.3f ',alt) ; 
                        % fprintf('\n') ;
                        out = Rates{nbpop+1} ;
                    else
                        % fprintf('YOUPI 3' ) ;
                        % fprintf('\n') ;
                        out = Rates{1} ; 
                    end
                end
            else
                if(nb~=nbpop)
                    alt = Rates{nbpop+1} ;
                    if(nbpop==4)
                        alt(4) = [] ; 
                        D = J(1,3) * J(2,2) - J(1,2) * J(2,3) ;
                    end
                    alt(1) = [] ; 
                    
                    if( all(alt>0) && D>0 ) 
                        fprintf('YOUPI 2 %d ', all(alt>0) ) ;
                        fprintf('%.3f ',alt) ; 
                        fprintf('\n') ; 
                        out = Rates{nbpop+1} ; 
                    else
                        % fprintf('YOUPI 3' ) ;
                        % fprintf('\n') ;
                        out = Rates{1} ;
                    end
                end
            end
        end        
    end
end