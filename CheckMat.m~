function [] = CheckMat(model,nbpop,n,K,IF_SPACE,Sigma,IF_Nk,DisplayOn,IF_SAVE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [] = Space_CheckMat(model,nbpop,n,K,Sigma,IF_Nk,DisplayOn,IF_SAVE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L = pi ;
    SCALING = 1.5 ;

    N = n*10000;
    % N = n*2500 ;
    nbN = nbNeuron(nbpop,n,IF_Nk,[]) ;
    Cpt = CptNeuron(nbpop,nbN) ;
        
    if strcmpi(IF_SPACE,'Ring') | strcmpi(IF_SPACE,'Gauss')  | strcmpi(IF_SPACE,'Exp')
        if(length(Sigma)==1 & nbpop==1)
            str_Matrix = sprintf(['../Connectivity/%dpop/N%d/K%d/%s/CrecI%.4f/Cij_Matrix.dat'],nbpop,n,K,IF_SPACE,Sigma)
        else
            if(nbpop==1)
                str_Matrix = sprintf(['../Connectivity/%dpop/N%d/K%d/%s/CrecI%.4f/Cij_Matrix.dat'],nbpop,n,K,IF_SPACE,Sigma)
            elseif(nbpop==2)
                if(length(Sigma)==nbpop)
                    if(length(K)==1)
                        str_Matrix = sprintf(['../Connectivity/%dpop/N%d/' ...
                                            'K%d/%s/CrecE%.4fCrecI%.4f/Cij_Matrix.dat'],nbpop,n,K,IF_SPACE,Sigma(1),Sigma(2)) 
                    else
                        str_Matrix = sprintf(['../%s/Connectivity/%dpop/N%d/' ...
                                            'KE%dKI%d/%s/CrecE%.4fCrecI%.4f/Cij_Matrix.dat'],model,nbpop,n,K(1),K(2),IF_SPACE,Sigma(1),Sigma(2)) 
                    end
                elseif(length(Sigma)>=nbpop)
                    str_Matrix = sprintf(['../%s/Connectivity/%dpop/N%d/' ...
                                        'K%d/%s/CrecEE%.4fCrecEI%.4fCrecIE%.4fCrecII%.4f/Cij_Matrix.dat'],model,nbpop,n,K,IF_SPACE,Sigma(1),Sigma(2),Sigma(3),Sigma(4)) 
                end
            else
                str_Matrix = sprintf(['../%s/Connectivity/%dpop/N%d/' ...
                                    'K%d/%s/CrecE%.4fCrecI%.4fCrecS%.4fCrecV%.4f/Cij_Matrix.dat'],model,nbpop,n,K,IF_SPACE,Sigma(1),Sigma(2),Sigma(3),Sigma(4)) 
            end
        end
    else
        str_Matrix = sprintf(['../Connectivity/%dpop/N%d/K%d/Cij_Matrix.dat'],nbpop,n,K)
    end

    fMatrix = fopen(str_Matrix,'rt') ;
    M = fread(fMatrix,'*int32') ;
    size(M) 
    M = reshape(M,N,N) ;
    M = double(M) ;
    size(M) 

    popList = ['E','I','S','X'] ;
    cl = {[1 0 0] [0 0 1] [0 1 0] [0.7 0.7 0.7]} ;
    
    for i=1:nbpop
        for j=1:nbpop
            % meanX = 0 ;
            meanY = 0 ;
                         
            A = M(Cpt(i)+1:Cpt(i+1),Cpt(j)+1:Cpt(j+1)) ; % A(i,j) j Pres to i Post, 

            for k=1:nbN(i)
                % meanX = meanX + sum(A(:,k)) ;
                meanY = meanY + sum(A(k,:)) ;
            end
            
            fprintf('Average nbPreS %s%s: %.3f \n',popList(j),popList(i),meanY./nbN(j))
            
            for k=1:nbN(j)
                P(i,j,k) = sum(diag(A,-k)) + sum(diag(A,nbN(j)-k)) ; % P(i,j) j to i
            end
        end
    end
    

    for i=1:nbpop        
        for j=1:nbpop
            
            Np = nbN(j) ;
            
            X = linspace(-pi,pi,Np)./L*SCALING ;                

            Pij = squeeze(P(i,j,1:Np)) ;
            % length(Pij)
            % length(X(1:Np/2)) 
            % length(X(Np/2+1:end))
            
            %figname=sprintf('Connection Probability %s%s Sigma %.4f',popList(j),popList(i),Sigma) ;
            if(length(Sigma)>1)
                figname=sprintf('Connection Probability %s%s Sigma %.4f',popList(j),popList(i),Sigma(i)) ;
                fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
            else
                figname=sprintf('Connection Probability %s%s Sigma %.4f',popList(j),popList(i),Sigma) ;
                fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
            end
            % plot(X,squeeze(P(i,j,:))./nbN(j),'color',cl{i})

            plot(X(1:Np/2),Pij(Np/2+1:end)./Np,'color',cl{i})
            plot(X(Np/2+1:end),Pij(1:Np/2)./Np,'color',cl{i})
             
            fit = smooth(X,Pij./Np,.1,'sgolay') ;

            % plot(X,fit,'linewidth',1,'color','k')

            plot(X(1:Np/2),fit(Np/2+1:end),'linewidth',1,'color','k')
            plot(X(Np/2+1:end),fit(1:Np/2),'linewidth',1,'color','k')

            xlabel('x')
            ylabel('Connection Probability')

            yLimits = get(gca,'YLim') ;
            
            if(length(Sigma)>1)
                SigScale = Sigma(i)./L*SCALING ;
            elseif(length(Sigma)>nbpop)                
                SigScale = Sigma(i+(j-1)*nbpop)./L*SCALING ;
            else
                SigScale = Sigma./L*SCALING ;
            end

            % plot([SigScale SigScale],[yLimits(1) yLimits(2)],'--','color',cl{i},'linewidth',1) 
            % plot([3*SigScale 3*SigScale],[yLimits(1) yLimits(2)],'--','color',cl{i},'linewidth',1) 
            % plot([-SigScale -SigScale],[yLimits(1) yLimits(2)],'--','color',cl{i},'linewidth',1) 
            % plot([-3*SigScale -3*SigScale],[yLimits(1) yLimits(2)],'--','color',cl{i},'linewidth',1) 
             
            if(strfind(model,'Binary'))
                xlabel('\phi')
                set(gca,'xtick',[-pi/2:pi/4:pi/2],'xticklabel',{'-\pi/2','-\pi/4', '0', '\pi/4', '\pi/2'})
                xlim([-pi/2,pi/2])
            end

            if(strfind(model,'IF'))
                xlabel('x')
                % set(gca,'xtick',[-1 0 1],'xticklabel',{'-L/2', '0', 'L/2'})
                % set(gca,'xtick',sort([ round(3*SigScale,2) round(-3*SigScale,2) -SCALING 0 SCALING] ))
                set(gca, 'FontSize', 10) 
                
                xlim([-SCALING,SCALING]) 
            end

            drawnow ;
            hold off ; 

            if(IF_SAVE)
                set(gca,'FontSize',2)
                figdir = sprintf(['./Figs/IF/%dpop/Connectivity/'],nbpop);
                try
                    mkdir(figdir)
                end
                ProcessFigure(fig, fullfile(figdir,figname),1.25) ;
            end

        end
    end
    
end