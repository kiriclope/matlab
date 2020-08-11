clear all ; 
GlobalVars 

Iext = ExternalInput(model,nbpop,dir) ;
J = ImportJab(model,nbpop,dir) ;

N=128;
K=1; 

IF_IEXT = '' ; 
l_DT = [.1 .075 .05 .025 .01 .0075 .005 .0025 .001] ;
l_dir = {'EULER','EULERI','RK2','RK2I'} ;
IDX=1 ;

cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
mk = {'x' '+' 's' 'd'} ; 

nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ; 

MFrates = linsolve(J,-Iext.') ; 
nbPoints = 0 ;

for k=1:length(l_dir) 
    nbPoints = 0 ;
    for i=1:length(l_DT) 
        
        dir = l_dir{k} ; 

        file = sprintf('/BenchMark/dt%.4f/Mean',l_DT(i)) ; 
        data = ImportData(model,nbpop,dir,file,N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iprtr) ;
 
        if(~isempty(data)) 
            nbPoints = nbPoints + 1 ; 
            popRate(k,i) = mean(data(IDX:end,2) ) ; 
            % for j=1:nbpop 
            %     popRate(k,j,i) = mean(data(IDX:end,j+1) ) ; 
            % end
        else
            popRate(k,i) = nan ; 
            % for j=1:nbpop 
            %     popRate(k,j,i) = nan ; 
            % end 
        end
                
    end
    
end

if(nbpop==2)
    fprintf('MF %2.f %2.f\n', MFrates(1), MFrates(2)) 
end

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');

fun = @(Coeff,xdata) Coeff(1) ./ (1 + Coeff(2) .* xdata) ; 
fun2 = @(Coeff,xdata) Coeff(1) ./ (1 + Coeff(2) .* xdata.^2) ; 

lb = []; 
ub = []; 

% fun = @(Coeff,xdata) abs( (Coeff(1)-xdata) ./ ( Coeff(2) .* xdata) ) ; 
% fun2 = @(Coeff,xdata) sqrt( abs( (Coeff(1)-xdata) ./ ( Coeff(2) .* xdata) ) ) ; 

if(nbPoints>1)
    figname=sprintf('BenchMark') ; 
    fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
    
    Eqn = 'a/(1+b*x)' ;
    % Eqn = 'abs((a-x)/(b*x))' ;
    for k=1:length(l_dir)
        if(k==4) 
            Eqn = 'a/(1+b*x.^2)' ; 
            % Eqn = 'sqrt(abs((a-x)/(b*x)))' ; 
        end
        fprintf('%s \n', l_dir{k}) 
        for j=1:nbpop 
            Rates = [] ;
            Rates = popRate(k,:) ; 
            % for i=1:length(popRate(k,j,:))
            %     Rates(i) = popRate(k,j,i) ;
            % end
            Rates
            % FIT = fit(l_DT',Rates',Eqn,options) 
            % FIT = fit(Rates',l_DT',Eqn) 
            % Coeff = coeffvalues(FIT) ; 
            if(k==4)
                Coeff = lsqcurvefit(fun2,[50 .1],l_DT,Rates,lb,ub,options) ;
                % Coeff = lsqcurvefit(fun2,[50 .1],Rates,l_DT,lb,ub,options) ; 
            else
                Coeff = lsqcurvefit(fun,[50 .1],l_DT,Rates,lb,ub,options) ; 
                % Coeff = lsqcurvefit(fun,[50 .1],Rates,l_DT,lb,ub,options) ; 
            end
            % Coeff = lsqcurvefit(fun,[50 .1],l_DT,Rates,lb,ub,options) ; 
            
            fprintf('a %.5f b %.5f \n',Coeff(1),Coeff(2)) 
            plot(l_DT(1:end-1), abs(Coeff(1)-Rates(1:end-1))./ ...
                 Rates(1:end-1),'color',cl{j},'marker',mk{k},'markerSize',5) 
            % plot(l_DT(1:end-1), abs(Rates(end)-Rates(1:end-1))./ ...
            %      Rates(1:end-1),'color',cl{j},'marker',mk{k},'markerSize',5) 
            % plot(l_DT(1:end-1), abs(1./Rates(1:end-1)-1./Rates(end)),'color',cl{j},'marker',mk{k},'markerSize',5) 
        end
    end
    xlabel('dt (ms)')
    ylabel('\epsilon_f  (Hz)')
    set(gca,'xscale','log')
    set(gca,'yscale','log')

end
