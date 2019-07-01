function [] = CheckMat2D(nbpop,n,K,IF_SPACE,Sigma,IF_Nk,DisplayOn,IF_SAVE)

    L = 1 ;
    SCALING = 1 ;

    N = n*10000 ;
    nbN = nbNeuron(nbpop,n,IF_Nk,[]) 
    
    N = 76800 ; 
    % nbN = 10000 ; 
    Cpt = CptNeuron(nbpop,nbN) 

    if(nbpop==1)
        str_Matrix = sprintf(['../Connectivity/%dpop/N%d/K%d/%s2D/CrecI%.4f/Cij_Matrix.dat'],nbpop,n,K,IF_SPACE,Sigma)
    else
        %str_Matrix = sprintf(['../Connectivity/%dpop/N%d/K%d/%s2D/CrecE%.4fCrecI%.4f/Cij_Matrix.dat'],nbpop,n,K,IF_SPACE,Sigma(1),Sigma(2)) 
        str_Matrix = sprintf(['../Connectivity/%dpop/N%d/K%d/Cij_Matrix.dat'],nbpop,n,K) 

    end

    fMatrix = fopen(str_Matrix,'rt') ;
    M = fread(fMatrix,'*float') ;
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

            % figname=sprintf('Connection Matrix %s%s ',popList(j),popList(i) );
            % fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
            % imagesc(A) ; 
            % xlabel('Post') 
            % ylabel('Pre') 
            
            for k=1:nbN(i)
                % meanX = meanX + sum(A(:,k)) ;
                meanY = meanY + sum(A(k,:)) ; 
            end
            
            %fprintf('Average nbPreS %s%s: %.3f %.3f \n',popList(i),popList(j),meanX./nbN(i),meanY./nbN(j))
            fprintf('Average nbPreS %s%s: %.3f \n',popList(j),popList(i),meanY./nbN(i))
            
            fprintf('Pb diag %.2f\n', sum(diag(A))/nbN(j) ) 
            
            Pb = [] ; 
            Pk = [] ; 
            for k=1:nbN(i)
                % Pb(k) = sum(diag(A,-k)) + sum(diag(A,nbN(j)-k)) ; % P(i,j) j to i
                Pk(k) = sum(A(k,:)) ;
            end 
            
            fprintf('Pk size %d\n', length(Pk) ) 
            
            % figname=sprintf('Connection Probability %s%s Sigma %.2f',popList(i),popList(j),Sigma(j)) ;
            % fig = figure('Name',figname,'NumberTitle','off') ; hold on ;

            % xlabel('Neuron')
            % ylabel('K_i')
            
            % plot(1:length(Pk),Pk) 

            figname=sprintf('In-degree_%s%s',popList(i),popList(j)) ;
            fig = figure('Name',figname,'NumberTitle','off') ; hold on ;

            xlabel('Neuron')
            ylabel('K_i')
            
            plot(1:length(Pk),Pk,'color',cl{i}) 

            if(IF_SAVE) 
                figdir = sprintf(['./Paper/%dpop/Connections/N%dK%dg%.2f'],nbpop,n,K,1) ;
                fprintf('Writing %s \n',figdir) 
                try
                    mkdir(figdir) ;
                end
                ProcessFigure(fig, fullfile(figdir,figname)) ;
            end

            figname=sprintf('In-degree_Dist_%s%s',popList(i),popList(j)) ;
            fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
            
            xlabel('K_i')
            ylabel('pdf')
            
            h = histogram( Pk , 27, 'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i},'EdgeAlpha',1,'Linewidth',2) ;

            if(IF_SAVE) 
                figdir = sprintf(['./Paper/%dpop/Connections/N%dK%dg%.2f'],nbpop,n,K,1) ;

                fprintf('Writing %s \n',figdir) 
                try
                    mkdir(figdir) ;
                end
                ProcessFigure(fig, fullfile(figdir,figname)) ;
            end
            
            % Pij = squeeze(Pb(i,j,:))./nbN(j) ; 
            % whos Pij
            % fprintf('Pij size %d\n', length(Pij) )

            % Pij = reshape(Pij, sqrt(length(Pij)), sqrt(length(Pij)) ) ; 
            
            % M = Pij ;
            
            % % M(1:sqrt(length(Pij))/2,sqrt(length(Pij))/2+1:end) = Pij(sqrt(length(Pij))/2+1:end,1:sqrt(length(Pij))/2) ;
            % % M(sqrt(length(Pij))/2:end,sqrt(length(Pij))/2:end) = Pij(1:sqrt(length(Pij))/2+1,1:sqrt(length(Pij))/2+1) ;
            % % M(sqrt(length(Pij))/2:end,1:sqrt(length(Pij))/2) = Pij(1:sqrt(length(Pij))/2+1,sqrt(length(Pij))/2+1:end) ;
            % % M(1:sqrt(length(Pij))/2,1:sqrt(length(Pij))/2+1) = Pij(sqrt(length(Pij))/2+1:end,sqrt(length(Pij))/2:end) ;
            
            % %M = interp2(M,'cubic') ;

            % figname=sprintf('Connection Probability %s%s Sigma %.2f',popList(j),popList(i),Sigma(i)) ;
            % fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
            
            % imagesc(M) ;
            % circle(100,100,3*Sigma(i)*100) ; 
            
            % % surf(M)
            % % shading interp
            % % view(90,0)

            % h = colorbar ; 
            % xlabel('X (mm)')
            % ylabel('Y (mm)')

            % set(gca,'xtick',linspace(0,length(M(:,1)),3),'xticklabel',{ linspace(0,1,3) })
            % set(gca,'ytick',linspace(0,length(M(:,1)),3),'yticklabel',{ linspace(0,1,3) })
            % xlim([0 length(M(:,1))])
            % ylim([0 length(M(:,1))])
 
            % set(gca,'xtick',[0 length(M(:,1))/4 length(M(:,1))/2 3*length(M(:,1))/4 length(M(:,1))], ...
            %         'xticklabel',{-L/2, -L/4, 0, L/4, L/2})
            
            % set(gca,'ytick',[0 length(M(:,1))/4 length(M(:,1))/2 3*length(M(:,1))/4 length(M(:,1))], ...
            %         'yticklabel',{-L/2, -L/4, 0, L/4, L/2})
           
        end
    end
    
    function h = circle(x,y,r)
        hold on ;
        th = 0:pi/50:2*pi ;
        xunit = r * cos(th) + x ;
        yunit = r * sin(th) + y ;
        h = plot(xunit, yunit, '--k','linewidth',1) ;
    end
    
end