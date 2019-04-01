clear all ;
GlobalVars

Iext = ExternalInput(model,nbpop,dir) ;    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

Iext(prtrPop) = Iext(prtrPop)+Iprtr ;

try
    data = ImportData(model,nbpop,dir,'Voltage',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ; 
    rates = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ; 
    for i=1:length(rates(1,:))-1
        IdvRates(i) = mean(rates(:,i+1)) ;
    end
end

tps = data(:,1)./1000 ;
% Volt = 3 .* data(:,2:end) - 40 ; 
Volt = data(:,2:end) ;

whos Volt
whos IdvRates

for i=1:nbpop
    
    % nId = randi([1+(i-1)*10 i*10]) ;
    % figname=sprintf('VoltTrace_%s_id%d_m_i=%.3f',popList(i), nId, IdvRates(nId+Cpt(i)) ) ;
    % fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
    
    % plot(tps,Volt(:,nId),'LineWidth',1,'color',cl{i})
    % xlabel('t (s)')
    % ylabel('Vm (mV)')
    % xlim([0 10])
    % drawnow;
    % if(IF_SAVE)
    %     figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
    %     fprintf('Writing %s \n',figdir)
    %     try
    %         mkdir(figdir)
    %     end
        
    %     ProcessFigure(fig, fullfile(figdir, figname), 3, [5*1.33*3, ...
    %                         3] );
    %     pause(.25);
    % end        
    % hold off ;

    for j=1:10
        figname=sprintf('VoltTrace_%s_id%d_m_i=%.3f',popList(i), j, IdvRates(j+Cpt(i)) ) ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        
        plot(tps,Volt(:,j+(i-1)*10),'LineWidth',1,'color',cl{i})
        xlabel('t (s)')
        ylabel('Vm (mV)')
        xlim([0 10])
        drawnow;
        if(IF_SAVE)
            figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
            fprintf('Writing %s \n',figdir)
            try
                mkdir(figdir)
        end
        
        ProcessFigure(fig, fullfile(figdir, figname), 3, [5*1.33*3,3] );
        end        
        hold off ;
    end
end