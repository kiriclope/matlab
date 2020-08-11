clear all ; 
GlobalVars 

g = 1 ; 
IextBL = ExternalInput(model,nbpop,dir) ;
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ; 
Cpt = CptNeuron(nbpop,nbN) ; 

v_Iprtr = v_Iprtr(1):v_Iprtr(2):v_Iprtr(3) ; 
FIGPERPOP = 0 ; 

J = ImportJab(model,nbpop,dir) ; 
Iext = IextBL ; 

u=[]; 
b=[]; 

% A = [2 2] ; 
% B = 10*[IextBL(1) IextBL(2)] ; 
% C = [0 1] ;

for i=1:length(v_Iprtr) 

    Iext(prtrPop) = IextBL(prtrPop) + v_Iprtr(i) ; 
    % Iext(prtrPop+1) = IextBL(prtrPop+1) + v_Iprtr(i) ; 
    % Iext(prtrPop+2) = IextBL(prtrPop+2) + v_Iprtr(i) ; 
    
    % for j=1:nbpop-1
    %     Iext(j) = max(( C(j) - A(j) ) * v_Iprtr(i) + B(j) , 0) ; 
    % end 

    % if(Iext(1)==0)
    %     Iext(2) = v_Iprtr(i) ;
    % end
    % [u b] = RateInputDist(model,nbpop,dir,g.*Iext,K,g,[],false) ; 

    RatesMF = BalRatesMF(model,nbpop,dir,Iext*.01,J,0) ; 
    [u b] = RateInputDist(model,nbpop,dir,Iext*.01,K,g,J,false,u,b) ; 
    Rates = QchAvgTF(u,b)*1000 ; 
    
    fprintf('I_opto ') 
    fprintf('%.3f ', v_Iprtr(i)) 


    fprintf('P_opto ') 
    %    fprintf('%.3f ', P0 .* ( exp( v_Iprtr(i) ./ I0 ) - 1 ) ) 
    if(IF_POWER==1) 
        %fprintf('%.3f ',  P0 .* ( exp(  v_Iprtr(i) ./ I0 ) - 1 ) ) ;
        fprintf('%.3f ', Pinf ./ ( 1 + exp(-( v_Iprtr(i)*2 - I0 )/Iinf ) ) ) ;
    elseif(IF_POWER==2)  
        fprintf('%.3f ',  P0 .* ( exp( v_Iprtr(i)*20 ./ I0 ) - 1 ) ) ; 
    end

    fprintf(' Rates ') 
    for j=1:nbpop 
        m(j,i) = Rates(j) ; 
        mu(j,i) = u(j) ;
        fprintf('%.3f ', m(j,i) ) 
    end 
    fprintf('\n') 

    fprintf('MF Rates ') 
    for j=1:nbpop 
        MFrates(j,i) = RatesMF(j)*1000 ; 
        fprintf('%.3f ', MFrates(j,i) ) 
    end 
    fprintf('\n') 

end

if(~FIGPERPOP) 
    figtitle = sprintf('%s_RatesVsIopto%s_MF',dir,popList(prtrPop)) ; 

    if(IF_POWER~=0) 
        figtitle = sprintf('%s_%d',figtitle,IF_POWER) ; 
    end
    if( ishandle( findobj('type','figure','name',figtitle) ) )
        fig = findobj('type','figure','name',figtitle) ; 
        fig = figure(fig); hold on ; 
    else
        if(IF_NORM) 
            figtitle = sprintf('%s_Norm',figtitle) ; 
        end
        fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
        xlabel('I_{opto} (\mu A . cm^{-2})') 
        if(IF_NORM) 
            ylabel('Norm. Activity') 
        else
            ylabel('Activity (Hz)') 
        end
    end
end

if IF_POWER==1
    %v_Iprtr = P0 .* ( exp( v_Iprtr ./ I0 ) - 1 )  ; 
    v_Iprtr = Pinf ./ ( 1 + exp(-( v_Iprtr*2 - I0 )/Iinf ) ) ;
    %v_Iprtr = exp( ( sqrt(K)*v_Iprtr*m0 -.4 ) / 0.2253 ) ;
elseif IF_POWER==2
    v_Iprtr = P0 .* ( exp( v_Iprtr *20 ./ I0 ) - 1 ) ; 
    %v_Iprtr = -P0 * log( 1 - v_Iprtr ./ I0 ) ; % Profile David 
end 

for i=nbpop:-1:1    

    if(FIGPERPOP)
        if(i==1 || i==2) 
            figtitle = sprintf('%s_RatesVsIopto%s_MF_EI',dir,popList(prtrPop)) ; 
        else 
            figtitle = sprintf('%s_RatesVsIopto%s_MF_SV',dir,popList(prtrPop)) ; 
        end
        fig = popPerFig(i,dir,figtitle) ;
        xlabel('I_{opto}') 
        ylabel('Norm. Activity') 
    end

    if(IF_POWER~=0) 
        figtitle = sprintf('%s_%d',figtitle,IF_POWER) ; 
    end
    if(IF_POWER~=0) 
        xlabel('P_{opto} (mW/mm^2)') 
    end
    
    if(IF_NORM)
        NormRates = m(i,:) ./ m(i,1) ; 
    else 
        NormRates = m(i,:) ;  
        ylabel('Activity (Hz)') 
    end

    if( (i==2 || i==4 ) && IF_NORM) 
        plot(v_Iprtr(IDX:end),ones(1, length(v_Iprtr(IDX:end)) ),'--','Color','k')
    
        if(IF_MF_RATES)
            if(i==2)
                for j=1:2
                    MFpop = MFrates(j,:)./MFrates(j,1) ;
                    plot(v_Iprtr, MFpop, '--','Color',cl{j}) 
                end
            elseif(i==4)
                for j=3:4
                    MFpop = MFrates(j,:)./MFrates(j,1) ;
                    plot(v_Iprtr, MFpop, '--','Color',cl{j}) 
                end
            end

        end

    end

    if(i==1 && IF_NORM)
        plot(v_Iprtr(IDX:end)  , NormRates(IDX:end), '-','Color',cl{i})
        % if(FIGPERPOP || nbpop==2) 
        %     plot(v_Iprtr(IDX:end)  , NormRates(IDX:end), '-','Color',cl{i})
        % else
        %     plot(v_Iprtr(IDX:end), NormRates(IDX:end), '+', ...e
        % 'MarkerEdgeColor',cl{i},'MarkerSize',6,'MarkerFaceColor', ...
        %     cl{i},'LineWidth', 1) 
        % end
    else                                
        plot(v_Iprtr(IDX:end)  , NormRates(IDX:end), '-','Color',cl{i}) 
        plot(v_Iprtr(IDX:end), NormRates(IDX:end), 'o', ...
             'MarkerEdgeColor',cl{i},'MarkerSize',1,'MarkerFaceColor','none','LineWidth', 1) 
    end
    %ylim([0.2 1.1])

    % if(i==3)
    %     plot(v_Iprtr(IDX:end), NormRates(IDX:end), '+', ...e
    %     'MarkerEdgeColor',cl{i},'MarkerSize',6,'MarkerFaceColor', ...
    %         cl{i},'LineWidth', 1) 
    % end
    drawnow ;
    
    if(IF_LOGSCALE)
        xlim([.01 100])
        set(gca,'Xscale', 'log')
        ylim([.01 10])
        set(gca,'Yscale', 'log') 
    else
        if(IF_LOGSCALEX)
            xlim([.01 100]) 
            set(gca,'Xscale', 'log')
            ylim([0 4])
        end
    end

    if(IF_SAVE & i==1)
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
        fprintf('Writing %s \n',figdir)
        try
            mkdir(figdir) ;
        end
        ProcessFigure(fig, fullfile(figdir,figtitle)) ;
    end

end

% figtitle = sprintf('%s_InputVsIopto%s_MF',dir,popList(prtrPop)) ;
% fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
% for i =1:nbpop 
%     plot(v_Iprtr(IDX:end)  , mu(i,IDX:end), '-','Color',cl{i}) 
% end 
% xlabel('I_{opto}')
% ylabel('\mu')
% hold off ; 