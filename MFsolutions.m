clear all ;
GlobalVars ;

IextBL = ExternalInput(model,nbpop,dir) ; 
Iext=IextBL;

v_Iprtr = v_Iprtr(1):v_Iprtr(2):v_Iprtr(3) ;

J = ImportJab(model,nbpop,dir) ; 
Conditions = (de2bi(1:2^(nbpop)-2 )) ; 

Baseline = linsolve(J,-IextBL.') ; 

figtitle = sprintf('%s_RatesVsIopto%s_MF',dir,popList(prtrPop)) ; 

if( ishandle( findobj('type','figure','name',figtitle) ) ) 
    fig = findobj('type','figure','name',figtitle) ; 
    fig = figure(fig); hold on ; 
else
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
    ylabel('Norm. Activity') 
    xlabel('I_{opto} (\mu A . cm^{-2})') 
end
R_PV_SOM = [] ;
I_PV_SOM = [] ;

for i=1:length(v_Iprtr) 

    Iext(prtrPop) = IextBL(prtrPop) + v_Iprtr(i) ; 
    Rates = linsolve(J,-Iext.') ; 

    fprintf('Iopto %.2f: ',v_Iprtr(i))
    if(all(Rates>0))
        fprintf('Full Balance: ') 
        fprintf('%.3f ', Rates) 
        fprintf('\n')
        
        for j=1:nbpop
            %plot(v_Iprtr(i), Rates(j)./Baseline(j),'color',cl{j},'marker','o','markerSize',mkSize)
            % plot(v_Iprtr(i), Rates(j),'color',cl{j},'marker','o','markerSize',mkSize)
        end
    else
        fprintf('No Full Balance \n ') 
    end

    for j=1:length(Conditions)
        % fprintf('Partial Balance (%d): ', j) 
        % fprintf('%d ', Conditions(j,:) ) 
        % fprintf('\n')

        Icond = Iext ;
        Jcond = J ;
        Rates = zeros(1,nbpop).' ; 
        Rcond = zeros(1,nbpop) ;            
        nonInputs = zeros(1,nbpop) ;            
        
        Idx = find(Conditions(j,:).') ;
        nonIdx = sort(find(Conditions(j,:).'==0),'descend') ;

        Rcond = Rcond(Idx) ;
        Icond = Iext(Idx) ;

        for k=1:length(nonIdx)
            Jcond(nonIdx(k),:) = [] ;
            Jcond(:,nonIdx(k)) = [] ;
        end

        if(length(Icond)>1)
            Rcond = linsolve(Jcond,-Icond.') ; 
        else
            if(Jcond~=0)
                Rcond = -Icond./Jcond ;
            else
                Rcond = 0 ;
            end
        end

        % fprintf('Rates ') 
        % fprintf('%.3f ',Rcond) 
        % fprintf('\n')
        
        % nonInputs = Iext(nonIdx).' + Jcond*Rcond ; 

        for k=1:length(Idx) 
            Rates(Idx(k)) = Rcond(k) ;
        end
        
        Inputs = Iext.' + J*Rates ; 
        nonInputs = Inputs(nonIdx) ;
        
        % for k=1:length(nonIdx)
        %     fprintf('IDX %d Input %.2f \n', nonIdx(k), nonInputs(k)) 
        % end

        if(all(Rcond>0) && all(nonInputs<=0) )
            fprintf('Partial Balance (%d): ', j) 
            fprintf('%d ', Conditions(j,:) ) 
            fprintf('\n') 
            fprintf('Rates ') 
            fprintf('%.3f ',Rcond) 
            fprintf('\n')

            if(length(Rcond)>2)
                R_PV_SOM = [R_PV_SOM max(Rcond(1),Rcond(2))] 
            else
                R_PV_SOM = [R_PV_SOM Rcond] 
            end
            I_PV_SOM = [I_PV_SOM  v_Iprtr(i)] ;

            for k=1:length(Rcond)
                %plot(v_Iprtr(i), Rcond(k)./Baseline(Idx(k)),'color',cl{Idx(k)},'marker','x','markersize',mkSize)
                %plot(v_Iprtr(i), Rcond(k),'color',cl{Idx(k)},'marker','x','markersize',mkSize)
            end
        end
    end
    fprintf('\n')

    if(nbpop==4)

        Icond = Iext ;
        Jcond = J ;
        Rcond = zeros(1,nbpop) ;            

        Jcond(4,:) = [] ;
        Jcond(:,4) = [] ;
        Icond(4) = [] ;

        Jcond(3,:) = [] ;
        Jcond(:,1) = [] ;
        Icond(3) = [] ;
        
        Rcond = linsolve(Jcond,-Icond.') ; 
        nonInputs = Iext(4) + J(4,2) * Rcond(1) + J(4,3) * Rcond(2) ;

        if(all(Rcond>0) && all(nonInputs<=0) ) 
            fprintf('me=1/sqrtK mi=O(1) ms=O(1) : ')
            fprintf('Rates ') 
            fprintf('%.3f ',Rcond) 
            % plot(v_Iprtr(i), Rcond(2),'color',cl{3},'marker','x','markersize',mkSize)
            % plot(v_Iprtr(i), Rcond(1),'color',cl{2},'marker','x','markersize',mkSize)

            
             % R_PV_SOM = [R_PV_SOM Rcond] 
             % I_PV_SOM = [I_PV_SOM  v_Iprtr(i)] ;

            % plot(v_Iprtr(i), Rcond(2)./Baseline(3),'color',cl{3},'marker','x','markersize',mkSize)
            % plot(v_Iprtr(i), Rcond(1)./Baseline(2),'color',cl{2},'marker','x','markersize',mkSize)
        end
        fprintf('\n')
    end

end

whos R_PV_SOM 
whos I_PV_SOM 

plot(I_PV_SOM, R_PV_SOM(1,:)./Baseline(2),'-','color',cl{2})

% plot(I_PV_SOM, R_PV_SOM(1,:)./Baseline(2),'--','color',cl{2})
% plot(I_PV_SOM, R_PV_SOM(2,:)./Baseline(3),'--','color',cl{3})
