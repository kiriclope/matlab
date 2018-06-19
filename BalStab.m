function [Relbd Imlbd] = BalStab(model,nbpop,dir,Iext,n,K,g,J,IF_Nk,IF_IDV,DisplayOn,IF_SAVE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function Rate_BalStab(nbpop,dir,Iext,K)
%% Slow synapses with single presynaptic time constant
%% Ta d Sa/dt = -Sa + TF(Iext_a + sum_b Jab Sb)
%% Ta d dSa/dt = -dSa + dTF(Iext_a + sum_b Jab Sb) sum_b Jab dSb
%% T d dS/dt = ( -I + Da . J) dS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    J = g.*J ;
    Tsyn = ImportTsyn(model,nbpop,dir) ;

    try
        lbd = eig(J./Tsyn) ;
        Relbd = real(lbd) ;
        Imlbd = imag(lbd) ;
    catch
        lbd = zeros(1,nbpop) ;
        Relbd = real(lbd) ;
        Imlbd = imag(lbd) ;
    end
    
    if(DisplayOn)
        fprintf('lbd J \n')
        for i=1:nbpop
            fprintf('%.3f + %.3fi\n',Relbd(i),Imlbd(i))
        end
    end

    N = n*10000 ;
    
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
        fprintf('Inputs u ')
        fprintf('%.3f ',u)
        fprintf('b ')
        fprintf('%.3f ',b)
        fprintf('\n')
        fprintf('Rates ')
        fprintf('%.3f ',QchAvgTF(u,b))
        fprintf('\n')
    end

    z = normrnd(0,1,1,N) ;
    sympref('HeavisideAtOrigin',0) ;

    for i=1:nbpop
        for j=1:nbpop
            h(i,Cpt(j)+1:Cpt(j+1)) = heaviside( u(i) + sqrt(b(i)).*z(Cpt(j)+1:Cpt(j+1) ) ) ;
        end
        Gain(i) = mean( h(i,:) ) ;
    end
    
    if(DisplayOn)
        fprintf('Gain ')
        for i=1:nbpop
            fprintf('%.3f ',Gain(i)) ;
        end
        fprintf('\n')
    end
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Rates dynamics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Tr = zeros(1,nbpop) ;
    Tr(1) = 2 ;
    Tr(2:nbpop) = 2 ;
    Id = eye(nbpop) ;

    for i = 1:nbpop
        for j = 1:nbpop
            G(i,j) = Gain(i) * J(i,j) / Tr(i) ;
            Id(i,j) = Id(i,j)/Tr(i) ;
        end
    end
    
    M = ( -Id./sqrt(K) + G ) ;
    % M = J.*mean(dh)./Tsyn 
    
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
        fprintf('lbd \n')
        for i=1:nbpop
            fprintf('%.3f + %.3fi\n ',Relbd(i),Imlbd(i))
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Individual Rates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(IF_IDV)
        J = J./sqrt(K) ;
        
        str_Matrix = sprintf(['../Rate/%dpop/Connectivity/N%d/K%d/Cij_Matrix.dat'],nbpop,n,K) 
        fMatrix = fopen(str_Matrix,'rt') ;
        M = fread(fMatrix,'*int32') ;
        size(M) ;
        M = reshape(M,N,N) ;
        M = double(M) ;
        size(M) ;

        sum(M(:,1)) ;
        sum(M(1,:)) ;
        
        for i=1:nbpop
            for j=1:nbpop
                M(Cpt(i)+1:Cpt(i+1),Cpt(j)+1:Cpt(j+1)) = J(i,j).*M(Cpt(i)+1:Cpt(i+1),Cpt(j)+1:Cpt(j+1)) ;
            end
        end
        for i=1:N
            M(:,i) = dh(i).*M(:,i) ;
        end
        
        whos M
        lbd = eig(M) ;
        figname=sprintf('Eig_%s_N%d_K%.0f',dir,n,K) ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        plot(real(lbd),imag(lbd),'.')
        
        xlabel('Re(\lambda)')
        ylabel('Im(\lambda)')

        if(IF_SAVE)
            xlim([-2 2])
            ylim([-2 2])
            set(gca,'XTick', [-2:2])
            set(gca,'XTick', [-2:2])

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