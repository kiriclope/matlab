function [Relbd Imlbd] = BalStabTsyn(model,nbpop,dir,Iext,n,K,g,J,IF_Nk,IF_IDV,DisplayOn,IF_SAVE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function Rate_BalStab(nbpop,dir,Iext,K)
%% Slow synapses with single presynaptic time constant
%% Ta d Sa/dt = -Sa + TF(Iext_a + sum_b Jab Sb)
%% Ta d dSa/dt = -dSa + dTF(Iext_a + sum_b Jab Sb) sum_b Jab dSb
%% T d dS/dt = ( -I + Da . J) dS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N = n*1000 ;

    if(isempty(Iext)) 
        Iext = ExternalInput(model,nbpop,dir) ;
    end
    if( nargin<8 | isempty(J) )
        J = ImportJab(model,nbpop,dir) ;
    end    
    if(nargin<9)
        IF_Nk = false ;
    end
    if(nargin<10)
        IF_IDV = false ;
    end
    if(nargin<11)
        DisplayOn = true ;
    end
    if(nargin<12)
        IF_SAVE = false ;
    end

    Tr = zeros(1,nbpop) ;
    Tr(1) = 3 ;
    Tr(2:nbpop) = 2 ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MF limit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J = g .* J ;

    Id = eye(nbpop) ;

    for i = 1:nbpop
        for j = 1:nbpop
            A(i,j) = J(i,j) / Tr(i) ; 
            Id(i,j) = Id(i,j) / Tr(i) ;
        end
    end

    try
        lbd = eig(-Id + sqrt(K) .* A ) ;
        Relbd = real(lbd) ;
        Imlbd = imag(lbd) ;
    catch
        lbd = zeros(1,nbpop) ;
        Relbd = real(lbd) ;
        Imlbd = imag(lbd) ;
    end
    
    
    if(DisplayOn)
        
        fprintf('MF\n Rates ')
        Rates = linsolve(J,-Iext.') ;
        for i=1:nbpop
            fprintf('%.3f | ', Rates(i)) ;
        end
        fprintf('\n') ;

        fprintf(' lbd ') ;
        for i=1:nbpop
            fprintf('%.3f + %.3fi | ', Relbd(i), Imlbd(i)) ;
        end
        fprintf('\n') ;

    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Self Consistant sol
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nbN = nbNeuron(nbpop,n,IF_Nk,[]) ;
    Cpt = CptNeuron(nbpop,nbN) ;
    
    count = 0 ;
    while count<1
        count = count+1 ;
        try
            [u b] = RateInputDist(model,nbpop,dir,Iext,K,1,J,false) ;
        catch
            fprintf('Error computing inputs\n') 
            u = zeros(1,nbpop) ;
            b = zeros(1,nbpop) ;
        end
    end

    if(DisplayOn)
        % fprintf('Inputs u ')
        % fprintf('%.3f ',u)
        % fprintf('b ')
        % fprintf('%.3f ',b)
        % fprintf('\n')
        fprintf('Finite K\n Rates ') 
        fprintf('%.3f | ',QchAvgTF(u,b))
        fprintf('\n')
    end

    z = normrnd(0,1,1,N) ;
    sympref('HeavisideAtOrigin',0) ;

    for i=1:nbpop
        for j=1:nbpop
            h(i,Cpt(j)+1:Cpt(j+1)) = heaviside( u(i) + sqrt(b(i)) .* z(Cpt(j)+1:Cpt(j+1) ) ) ;
        end
        Gain(i) = mean( h(i,:) ) ;
    end
    
    % if(DisplayOn)
    %     fprintf('Gain ')
    %     for i=1:nbpop
    %         fprintf('%.3f ',Gain(i)) ;
    %     end
    %     fprintf('\n')
    % end
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Rates dynamics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Id = eye(nbpop) ;

    for i = 1:nbpop
        for j = 1:nbpop
            G(i,j) = Gain(i) * J(i,j) / Tsyn(i,j) ;
            Id(i,j) = Id(i,j) / Tsyn(i,j) ;
        end
    end
    
    M = ( -Id + sqrt(K) .* G ) ;
    
    try
        lbd = eig(M) ;
        Relbd = real(lbd) ;
        Imlbd = imag(lbd) ;
    catch
        lbd = zeros(1,nbpop) ;
        Relbd = real(lbd) ;
        Imlbd = imag(lbd) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Synaptic dynamics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Id = eye(nbpop) ;
    % M = ( -Id./sqrt(K) + J.*mean(dh) )./Tsyn ;

    % Mab = -hab + K * Jab * mean( dF(htotb) ) * ( hbe - hbi - hbs - hbv) ;
    % for i=1:nbpop
    %     for j=1:nbpop
    %         M(i,j) = -1 + K * J(i,j) * mean( dh( Cpt(j)+1:Cpt(j+1) ) ) *
    % % M = J.*mean(dh)./Tsyn 
    
    % try
    %     lbd = eig(M) ;
    %     Relbd = real(lbd) ;
    %     Imlbd = imag(lbd) ;
    % catch
    %     lbd = zeros(1,nbpop) ;
    %     Relbd = real(lbd) ;
    %     Imlbd = imag(lbd) ;
    % end

    if(DisplayOn)
        fprintf(' lbd ') ;
        for i=1:nbpop
            fprintf('%.3f + %.3fi | ', Relbd(i), Imlbd(i)) ;
        end
        fprintf('\n') ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Individual Rates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(IF_IDV)
        fprintf('From Cij Matrix \n') ;

        J = J./sqrt(K) ;
        
        str_Matrix = sprintf(['../Connectivity/%dpop/N%d/K%d/Cij_Matrix.dat'],nbpop,n,K) ;
        disp(str_Matrix) ;
        fMatrix = fopen(str_Matrix,'rt') ;
        M = fread(fMatrix,'*float') ;
        size(M) ;
        M = reshape(M,N,N) ;
        M = double(M) ;
        size(M) ;

        sum(M(:,1)) ;
        sum(M(1,:)) ;

        Id = eye(N) ;
        
        for i=1:nbpop
            for j=1:nbpop
                M(Cpt(i)+1:Cpt(i+1),Cpt(j)+1:Cpt(j+1)) = J(i,j) ./ Tr(i) .* M(Cpt(i)+1:Cpt(i+1),Cpt(j)+1:Cpt(j+1)) ;
                Id(Cpt(i)+1:Cpt(i+1),Cpt(j)+1:Cpt(j+1)) = Id(Cpt(i)+1:Cpt(i+1),Cpt(j)+1:Cpt(j+1)) ./ Tr(i) ;
            end
        end

        for i=1:nbpop
            for j=1:N
                M(Cpt(i)+1:Cpt(i+1),j) = h(i,j) .* M(Cpt(i)+1:Cpt(i+1),j) ;
            end
        end

        lbd = eig(-Id + M) ;
        for i=1:nbpop
            idx = find( lbd==max(lbd) ) ;
            lbdMax(i) = lbd(idx) ;
            lbd(idx) = [] ;
        end

        figname=sprintf('Eig_%s_N%d_K%.0f',dir,n,K) ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        plot(real(lbd),imag(lbd),'.')

        fprintf('lbdMax ')
        for i=1:nbpop
            fprintf('%.3f + %.3fi | ',real(lbdMax(i)),imag(lbdMax(i))) ;
            plot(real(lbdMax(i)),imag(lbdMax(i)),'r*','markersize',2)
        end
        fprintf('\n') ;
        
        xlabel('Re(\lambda)')
        ylabel('Im(\lambda)')
        
        % xlim([-2 2])
        % ylim([-2 2])
        % set(gca,'XTick', [-2:2])
        % set(gca,'XTick', [-2:2])
                     
        if(IF_SAVE)
            
            figdir = sprintf(['./Figs/%dpop/%s/stab/'],nbpop,dir);
            try
                mkdir(figdir)
            end

            ProcessFigure(fig, fullfile(figdir,figname)) ;
        end
            
    end     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = phi(x)
        out = exp(-x.^2./2)./sqrt(2.*pi);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = Phi(x)
        out = .5.*(1+erf( x./sqrt(2) ) ) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = QchAvgTF(u,b)
        if(b>0)
            out = u.*Phi( u./sqrt(b) ) + sqrt(b).*phi( u./sqrt(b) ) ;
        else
            out = u ;
        end
    end
   
end