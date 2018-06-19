function [Rates Det] = BalRatesMF(model,nbpop,dir,Iext,J,DisplayOn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Rates Det] = BalRatesMF(model,nbpop,dir,Iext,J,DisplayOn)
% Computes Balance state rates in the large N, large K limit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(nargin<4 | isempty(Iext))
        Iext = ExternalInput(model,nbpop,dir) ;
    end
    if(nargin<5 | isempty(J))       
        J = ImportJab(model,nbpop,dir) ;
    end
    if(nargin<6)
        DisplayOn = 0 ;
    end

    Det = det(J) ;
    Rates = linsolve(J,-Iext.') ;
    
    if(DisplayOn)
        fprintf('Rates ')
        fprintf('%.3f ',Rates)
        fprintf('\n')
        fprintf('Det ')
        fprintf('%.3f ',Det)
        fprintf('\n')
    end
end