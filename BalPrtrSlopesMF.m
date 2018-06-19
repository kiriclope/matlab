function dmdI = BalPrtrSlopesMF(model,nbpop,dir,Iext,J,IF_NUMERIC) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function dmdI = BalPrtrSlopesMF(model,nbpop,dir,Iext,J,IF_NUMERIC)
% Compute susceptibility to a perturbation of the 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    if(nargin<5 | isempty(J) )
        J = ImportJab(model,nbpop,dir) ;
        if(~IF_NUMERIC)
            J = sym('J%d%d',nbpop) ;
            C = CreateCab(J) ;
            J = J.*C ;
        end
    end

    Isym = sym('I',[1 nbpop]) ;
    
    if nbpop>2 
        Isym(3) = 0 ;
    end

    dm = sym('dm',nbpop) ;
    m = sym('m',nbpop) ;

    m = linsolve(J,-Isym.') ;

    syms I1 I2 I3 I4
    for i=1:nbpop
        Isym = sym('I',[1 nbpop]) ;
        Isym(Isym~=Isym(i)) = 0 ; 
        Isym(i) = 1 ;
        for j=1:nbpop
            if(nbpop==2)
                dm(i,j) = subs(m(j),[I1 I2], Isym) ;
            elseif(nbpop==3)
                dm(i,j) = subs(m(j),[I1 I2 I3], Isym) ;
            elseif(nbpop==4)
                dm(i,j) = subs(m(j),[I1 I2 I3 I4], Isym) ;
            end
        end
    end
    if(IF_NUMERIC)
        if(nbpop==2)
            m = subs(m,[I1 I2],Iext) ;
            dmdI = subs(dm(2,:),[I1 I2],Iext) ;
        elseif(nbpop==3)
            m = subs(m,[I1 I2 I3],Iext) ;
            dmdI = subs(dm(2,:),[I1 I2 I3],Iext) ;
        elseif(nbpop==4)
            m = subs(m,[I1 I2 I3 I4],Iext) ;
            dmdI = double(subs(dm(2,:),[I1 I2 I3 I4],Iext)) ;
        end
        for i=1:nbpop
            if(m(i)>0)
                dmdI(i) = dmdI(i)./m(i) ;
            else
                dmdI(i) = nan ;
            end
        end
    end
end
