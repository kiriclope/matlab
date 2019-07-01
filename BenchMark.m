clear all ; 
GlobalVars 

Iext = ExternalInput(model,nbpop,dir) ;
J = ImportJab(model,nbpop,dir) ;

IF_IEXT = '' ; 
l_DT = [.1 .05 .01 .005 .001] ; 
l_dir = {'EULER', 'RK2','EULERI','RK2I'} ;
IDX=5 ;

cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
mk = {'x' 'o' 's' 'd'} ; 

nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ; 

MFrates = linsolve(J,-Iext.') ; 
nbPoints = 0 ;

for k=1:length(l_dir) 
    nbPoints = 0 ;
    for i=1:length(l_DT) 
        
        dir = l_dir{k} ; 

        file = sprintf('/BenchMark/dt%.3f/Mean',l_DT(i)) ;
        data = ImportData(model,nbpop,dir,file,N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iprtr) ; 
 
        if(~isempty(data)) 
            nbPoints = nbPoints + 1 ;                
            for j=1:nbpop 
                popRate(k,j,i) = mean(data(IDX:end,j+1) ) ; 
            end
        else
            for j=1:nbpop 
                popRate(k,j,i) = nan ; 
            end        
        end
                
    end
    
end

fprintf('MF %2.f %2.f\n', MFrates(1), MFrates(2)) 

if(nbPoints>1)
    figname=sprintf('BenchMark') ; 
    fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
    
    Eqn = 'a/(1+b*x)' ;
    for k=1:length(l_dir) 
        if(k==4) 
            Eqn = 'a/(1+b*x.^2)' ; 
        end
        fprintf('%s \n', l_dir{k}) 
        for j=1:nbpop 
            f = coeffvalues(fit(l_DT',squeeze(popRate(k,j,:)),Eqn)) ; 
            fprintf('a %.2f b %.2f \n',f(1),f(2)) 
            plot(l_DT, squeeze(abs(f(1)-popRate(k,j,:))./popRate(k,j,:)),'-','color',cl{j},'marker',mk{k},'markerSize',5) 
            %plot(l_DT, squeeze(abs(MFrates(j)-popRate(k,j,:))./ ...
            %popRate(k,j,:)),'-','color',cl{j},'marker',mk{k},'markerSize',5) 
            % plot(l_DT, squeeze(popRate(k,j,:)),'-','color',cl{j},'marker',mk{k},'markerSize',5) 
        end
    end
    xlabel('dt (ms)')
    ylabel('\epsilon_R  (Hz)')

    set(gca,'xscale','log')
    set(gca,'yscale','log')
    % figname=sprintf('BenchMark2') ; 
    % fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
    
    % for j=1:nbpop 
    %     plot(l_DT, abs(popRate(j,:)-popRate2(j,:)),'-o','color',cl{j},'markerSize',5)  
    % end

    % xlabel('dt (ms)')
    % ylabel('\Delta Rates (Hz)')

end
