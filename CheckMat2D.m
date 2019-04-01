function [] = CheckMat2D(nbpop,n,K,IF_SPACE,Sigma,IF_Nk,DisplayOn,IF_SAVE)
    L = 2 ;
    SCALING = 1 ;

    N = n*10000 ;
    nbN = nbNeuron(nbpop,n,IF_Nk,[]) 
    Cpt = CptNeuron(nbpop,nbN) 

    if(nbpop>=2)
        sqrtNp = [ sqrt(nbN(1)) sqrt(nbN(2)) ] ;
    else 
        sqrtNp = sqrt(nbN(1)) ;
    end

    if(nbpop==1)
        str_Matrix = sprintf(['../Connectivity/%dpop/N%d/K%d/%s2D/CrecI%.4f/Cij_Matrix.dat'],nbpop,n,K,IF_SPACE,Sigma)
    else
        str_Matrix = sprintf(['../Connectivity/%dpop/N%d/K%d/%s2D/CrecE%.4fCrecI%.4f/Cij_Matrix.dat'],nbpop,n,K,IF_SPACE,Sigma(1),Sigma(2)) 
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

            for k=1:nbN(i)
                % meanX = meanX + sum(A(:,k)) ;
                meanY = meanY + sum(A(k,:)) ;
            end
            
            %fprintf('Average nbPreS %s%s: %.3f %.3f \n',popList(i),popList(j),meanX./nbN(i),meanY./nbN(j))
            fprintf('Average nbPreS %s%s: %.3f \n',popList(j),popList(i),meanY./nbN(j))
            
            for k=1:nbN(i)
                P(i,j,k) = sum(diag(A,-k)) + sum(diag(A,nbN(j)-k)) ; % P(i,j) j to i
            end
        end
    end

    for i=1:nbpop
        for j=1:nbpop            
            
            Pij = squeeze(P(i,j,:))./nbN(j) ;
            Pij = reshape(Pij,sqrtNp(j),sqrtNp(j)) ;
            
            M = Pij ;
            M(1:sqrtNp(j)/2,sqrtNp(j)/2+1:end) = Pij(sqrtNp(j)/2+1:end,1:sqrtNp(j)/2) ;
            M(sqrtNp(j)/2:end,sqrtNp(j)/2:end) = Pij(1:sqrtNp(j)/2+1,1:sqrtNp(j)/2+1) ;
            M(sqrtNp(j)/2:end,1:sqrtNp(j)/2) = Pij(1:sqrtNp(j)/2+1,sqrtNp(j)/2+1:end) ;
            M(1:sqrtNp(j)/2,1:sqrtNp(j)/2+1) = Pij(sqrtNp(j)/2+1:end,sqrtNp(j)/2:end) ;

            M = interp2(M,'cubic') ;

            figname=sprintf('Connection Probability %s%s Sigma %.2f',popList(j),popList(i),Sigma(i)) ;
            fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
            
            imagesc(M) ;
            circle(100,100,3*Sigma(i)*100) ; 
            h = colorbar ;

            xlabel('X (mm)')
            ylabel('Y (mm)')

            set(gca,'xtick',linspace(0,length(M(:,1)),3),'xticklabel',{ linspace(0,1,3) })
            set(gca,'ytick',linspace(0,length(M(:,1)),3),'yticklabel',{ linspace(0,1,3) })
            xlim([0 length(M(:,1))])
            ylim([0 length(M(:,1))])

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