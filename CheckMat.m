clear all ; 
GlobalVars

str_Matrix = ImportData('Connectivity', nbpop, dir, 'Cij_Matrix', N, K, g, IF_RING, Crec, [], '', 1, 1) 
n = N ;
[nbN N] = nbNeuron(nbpop,N,IF_Nk,[]) ; 
Cpt = CptNeuron(nbpop,nbN) ; 

fMatrix = fopen(str_Matrix,'rt') ; 
M = fread(fMatrix,'*float') ; 
size(M) ; 
M = reshape(M,N,N) ; 
M = double(M)' ; 
size(M) ;

popList = ['E','I','S','X'] ;
cl = {[1 0 0] [0 0 1] [0 1 0] [0.7 0.7 0.7]} ;
if(nbpop==1) 
    popList = 'I' ; 
    cl = {[0 0 1]} ; 
end

if(~FIGPERPOP)
    figname=sprintf('Connection_Probability_CrecE_%.4fCrecI%.4f',Crec(1),Crec(2));
    fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
    xlabel('X (mm)')
    ylabel('Pb')
end

for i=1:nbpop
    for j=1:nbpop
        meanX = 0 ; 
        meanY = 0 ; 
        
        A = M(Cpt(i)+1:Cpt(i+1),Cpt(j)+1:Cpt(j+1)) ; % A(i,j) j Pres to i Post, 

        % figname=sprintf('Connection Matrix %s%s ',popList(i),popList(j) );
        % fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
        % xlabel('Post') 
        % ylabel('Pre') 
        % imagesc(A) ; 
        
        for k=1:nbN(i)
            % meanX = meanX + sum(A(:,k)) ; 
            meanY = meanY + sum(A(k,:)) ; 
        end
        
        % fprintf('Average nbPreS %s%s: %.3f %.3f \n',popList(i),popList(j),meanX./nbN(i),meanY./nbN(j))
        fprintf('Average nbPreS %s%s: %.3f \n',popList(i),popList(j),meanY./nbN(i))         
        fprintf('Pb diag %.2f\n', sum(diag(A))/nbN(j) ) 
        
        Pb = [] ; 
        Pk = [] ; 
        IdvCorr = [] ;

        for k=1:nbN(i)
            % Pb(k) = (mean(diag(A,k)) + mean(diag(A,nbN(i)-k)) )/2 ; 
            Pb(k) = (sum(diag(A,k)) + sum(diag(A,nbN(i)-k)) ) ./ nbN(i) ; 
            % IdvCorr = [ IdvCorr; xcorr( A(k,:) ) ];
            Pk(k) = sum(A(k,:)) ; 
        end 
        
        fprintf('Pk size %d\n', length(Pk) ) 

        % Pb = mean(IdvCorr) ; 


        if(FIGPERPOP)
            figname=sprintf('Connection_Probability_%s%s_Crec_%.4f',popList(i),popList(j),Crec(j));
            fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
            xlabel('X (mm)')
            ylabel('Pb')
        end

        X = linspace(-L/2,L/2,length(Pb)) ; 
        
        plot(X(1:length(Pb)/2),Pb(length(Pb)/2+1:end),'color',cl{j}) 
        plot(X(length(Pb)/2+1:length(Pb)),Pb(1:length(Pb)/2),'color',cl{j}) 

        plot(4*Crec(j) *[1 1],[0 1], '--','Color',cl{i}) 
        plot(-4*Crec(j) *[1 1],[0 1], '--','Color',cl{i}) 
        
        % figname=sprintf('In-degree_%s%s',popList(i),popList(j)) ;
        % fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
        % xlabel('Neuron')
        % ylabel('K_i')
        
        % plot(1:length(Pk),Pk,'color',cl{j}) 

        % if(IF_SAVE) 
        %     figdir = sprintf(['./Paper/%dpop/Connections/N%dK%dg%.2f'],nbpop,n,K,1) ; 
        %     fprintf('Writing %s \n',figdir) 
        %     try
        %         mkdir(figdir) ;
        %     end
        %     ProcessFigure(fig, fullfile(figdir,figname)) ;
        % end
        
        % figname=sprintf('In-degree_Dist_%s%s',popList(i),popList(j)) ;
        % fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
        
        % xlabel('K_i') 
        % ylabel('pdf') 
        
        % h = histogram( Pk , 27, 'Normalization', 'pdf'
        % ,'DisplayStyle','stairs','EdgeColor',cl{j},'EdgeAlpha',1,'Linewidth',2)
        % ;        
    end
end

ylim([0 1])
xlim([0 L/2])

if(IF_SAVE) 
    figdir = sprintf(['./Paper/%dpop/Connections/N%dK%dg%.2f/%s/CrecE%.4fCrecI%.4f'],nbpop,n,K,1,IF_RING,Crec(1),Crec(2)) ; 
    
    fprintf('Writing %s \n',figdir) 
    try
        mkdir(figdir) ;
    end
    ProcessFigure(fig, fullfile(figdir,figname)) ;
end
