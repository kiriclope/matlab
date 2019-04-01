clear all ;
GlobalVars

IextBL = ExternalInput(model,nbpop,dir) ; 
v_Jab = v_Jab(1):v_Jab(2):v_Jab(3) ; 

J = ImportJab(model,nbpop,dir) ; 
G = abs(J) ;
Iext = IextBL ; 
Iext(2) = IextBL(2) + Iprtr ; 

for i=1:length(v_Jab) 
    J(1,1) = v_Jab(i) ; 

    Rates = linsolve(J,-IextBL.') ; 
    RatesPrtr = linsolve(J,-Iext.') ; 
    
    fprintf('Jab ') 
    fprintf('%.3f ', v_Jab(i)) 

    fprintf(' Rates ') 
    for j=1:nbpop 
        m(j,i) = Rates(j) ; 
        fprintf('%.3f ', m(j,i) ) 
    end 
    fprintf('\n') 

    fprintf(' Khi ') 
    for j=1:nbpop 
        khi(j,i) = ( RatesPrtr(j) - Rates(j) ) ./ Iprtr ./ Rates(j);
        fprintf('%.3f ', khi(j,i) ) 
    end 
    fprintf('\n') 

    % khi(1,i) = - G(1,2) ./ ( G(2,2) * Iext(1) - G(1,2) * Iext(2) ) ;
    % khi(2,i) = - J(1,1) ./ ( G(2,1) * Iext(1) - J(1,1) * Iext(2) ) ;

end

figtitle = sprintf('%s_RatesVsJab_MF',dir) ; 
fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
xlabel('J_{EE}')
ylabel('Rates ')

% xlim([.01 100])
% set(gca,'Xscale', 'log')
% ylim([.01 10])
% set(gca,'Yscale', 'log') 

for i=1:nbpop   
    NormRates = m(i,:) ;
    plot(v_Jab, NormRates, '-','Color',cl{i})    
end

if(IF_SAVE)
    figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
    fprintf('Writing %s \n',figdir)
    try
        mkdir(figdir) ;
    end
    ProcessFigure(fig, fullfile(figdir,figtitle)) ;
end


figtitle = sprintf('%s_KhiVsJab',dir) ; 
if( ishandle( findobj('type','figure','name',figtitle) ) )
    fig = findobj('type','figure','name',figtitle) ; 
    fig = figure(fig); hold on ; 
else
    fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
    xlabel('J_{EE}')
    ylabel('\chi ')
end
plot(v_Jab, zeros(1,length(v_Jab)), '--','Color','k')

% xlim([.01 100])
% set(gca,'Xscale', 'log')
% ylim([.01 10])
% set(gca,'Yscale', 'log') 

for i=1:nbpop
    NormKhi = khi(i,:) ; 
    plot(v_Jab, NormKhi, '--','Color',cl{i}) 
end

% plot([0.05 0.05],[-2 2], '--','Color','k')
% plot([0.4 0.4],[-2 2], '--','Color','k')

    if(IF_SAVE)
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_DATA) ;
        fprintf('Writing %s \n',figdir)
        try
            mkdir(figdir) ;
        end
        ProcessFigure(fig, fullfile(figdir,figtitle)) ;
    end

