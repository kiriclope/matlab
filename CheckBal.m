function CheckBal(model,nbpop,dir,Iext,K,file,n,g,IF_Nk,IF_RING,Crec,Cff,IF_IEXT,nPrtr,IF_SAVE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Rate_CheckBal(nbpop,dir,Iext,K,file,n,g,IF_Nk,IF_PE,nPrtr,Iprtr,IF_RING,Cff,Crec,IF_IEXT)
% Compares simulation to analytics of the balance state
% if : - file=IdvInputs, plots input distributions
%      - file=IdvRates, plots rate distributions
% Utils : - ImportData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(isempty(Iext))
        Iext = ExternalInput(model,nbpop,dir) ;
    end

    if(strfind(dir,'2pop'))
        cl = {[1 0 0] [0 0 1] [0 1 1] } ;
    else
        cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ;
    end

    p = [] ;
    nbN = nbNeuron(nbpop,n,IF_Nk,p) ;
    Cpt = CptNeuron(nbpop,nbN) ;
    THRESHOLD = .1 ;

    nb = 5 ;
    

    if(nargin<8)
        IF_Nk = false ;
    end

    if(nargin<10)
        IF_RING = false ;
        Crec = 0 ;
        Cff = 0 ;
    end

    if(nargin<13)
        IF_IEXT=false ;
        nPrtr = 2 ;
    end

    if(IF_IEXT==false)
        nPrtr = 2 ;
    end

    Iprtr=Iext(nPrtr) ;

    if(nargin<15)
        IF_SAVE=false ;
    end

    if strcmpi(file,'rpop')
        try
            data = ImportData(model,nbpop,dir,'Population_Rates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ;
            for i=2:length(data(1,:))
                PopRates(i-1,:) = data(:,i) ;
                meanRate(i-1) = mean(data(:,i)) ;
            end
            tps = data(:,1)./1000 ;
        catch
            fprintf('DATA NOT FOUND\n')
            for i=2:length(data(1,:))
                PopRates(i-1) = nan ;
                meanRate(i-1) = nan ;
            end
            tps = [] ;
        end
                
        figname=sprintf('PopulationRates') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

        for i=1:nbpop
            if(i==1 | i==4)
                plot(tps,PopRates(i,:),'--','color',cl{i})
            else
                plot(tps,PopRates(i,:),'color',cl{i})
            end
        end    
        xlabel('t (a.u.)')
        ylabel('m_\alpha')
        drawnow ;

        fprintf(' Mean Rates \n')
        fprintf(' %.3f ', meanRate)
        fprintf('\n')
    end

    if strcmpi(file,'Movie')
        set(0, 'DefaultFigureRenderer', 'zbuffer') ;

        try
            data = ImportData(model,nbpop,dir,'Raster',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ; 
            Spikes = sortrows(data,2) ;
            Spikes(:,2) = Spikes(:,2)./1000 ;
        catch
            fprintf('FILE NOT FOUND\n')
        end

        [nbSpk tIdx]= hist(Spikes(:,2),unique( Spikes(:,2) ) ) ;
        CumSumSpk = [0 cumsum(nbSpk)] ;
        
        for i=1:nbpop
            for j=1:nbN(i)

                ix(i,j) = mod( j+1, sqrt(nbN(i)) ) ;
                X(i,j) = ix(i,j)/sqrt(nbN(i)) ;
                
                iy(i,j) = ceil( j/sqrt(nbN(i)) ) ;
                Y(i,j) = iy(i,j)/sqrt(nbN(i)) ;

                fprintf('ix %d iy %d X %.3f Y %.3f \r',ix(i,j),iy(i,j),X(i,j),Y(i,j))
            end
        end
        fprintf('\n')

        tw = 10 ;
        tl = 1 ;

        nIdx = {} ;
        for i=1:(length(CumSumSpk)-1)/tl
            fprintf('# %d nbSpk %.3f \r',tIdx(i),nbSpk(i)) 
            nIdx = [nIdx ; Spikes(CumSumSpk(i)+1:CumSumSpk(i)+nbSpk(i)-1,1).'] ;
        end
        fprintf('\n') 
        
        % nIdx{1}
        % [nIdx{1:2}]

        popList = ['E' 'I' 'S' 'X'] ;

        writerObjE = VideoWriter('Eraster.avi');
        open(writerObjE);

        writerObjI = VideoWriter('Iraster.avi');
        open(writerObjI);

        nframe = length(1:tw:(length(CumSumSpk)-1)/tl-tw) ;
        % vectMov = [] ;
        % mov(1:nbpop,1:nframe) = struct('nbpop',[],'cdata',[],'colormap',[]) ;
        
        % for i=1:nbpop
        %     vectMov = [vectMov ; mov(1:nbpop,1:nframe)]
        % end

        figname = 'XY raster' ;
        fig = figure('Name',figname,'NumberTitle','off','Visible','off') ;
        plot(1,1,'.')
        set(gca,'nextplot','replacechildren','Visible','off')
        xlabel('X (mm)')
        ylabel('Y (mm)')
        xlim([0 1])
        ylim([0 1])
        drawnow ;

        % figI = figure(2) ;
        % plot(1,1,'.')
        % set(gca,'nextplot','replacechildren','Visible','off')
        % xlabel('X (mm)')
        % ylabel('Y (mm)')
        % xlim([0 1])
        % ylim([0 1])
        % drawnow ;
        % hAxI = gca;

        % figure('visible','off'); 
        for i=1:tw:(length(CumSumSpk)-1)/tl-tw
            Xax = cell(nbpop) ;
            Yax = cell(nbpop) ;
            nId = [nIdx{i:i+tw}] ;

            for j=1:length(nId)
                for k=1:nbpop-1
                    if(Cpt(k)<=nId(j) & nId(j)<Cpt(k+1))
                        fprintf('t %f %s nId %d',tIdx(i),popList(k),nId(j))
                        fprintf(' ix %d iy %d \r',ix(k,nId(j)-Cpt(k)+1),iy(k,nId(j)-Cpt(k)+1))
                        Xax{k} = [Xax{k} X(k, nId(j)-Cpt(k)+1 )] ;
                        Yax{k} = [Yax{k} Y(k, nId(j)-Cpt(k)+1 )] ;
                    end
                end
            end
            
            for k=1:nbpop
                if(k==1)
                    plot(Xax{k},Yax{k},'.','color',cl{k},'MarkerSize',5)
                    drawnow                     
                    writeVideo(writerObjE,getframe(gca)) ;
                else
                    plot(Xax{k},Yax{k},'.','color',cl{k},'MarkerSize',5)
                    drawnow
                    writeVideo(writerObjI,getframe(gca)) ;
                end
            end
        end

        close(writerObjE);
        close(writerObjI);
        fprintf('\n')
        % movie2avi(mov, 'raster.avi', 'compression', 'None') ;
    end


    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%
    % Spike Times      %
    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%

    if strcmpi(file,'Raster')
        try
            data = ImportData(model,nbpop,dir,file,n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr);
            Spikes = sortrows(data) ; 
            Spikes(:,2) = Spikes(:,2)./1000. ;
        catch
            fprintf('FILE NOT FOUND\n')
        end

        for i=1:nbpop
            Spikes( Spikes(:,1)<Cpt(i+1) & Spikes(:,1)>Cpt(i+1)-0.75*nbN(i), : ) = [] ;
        end

        reSizeIdx = 1:length( unique( Spikes ) ) ;

        [nbSpk nidx]= hist(Spikes(:,1),unique( Spikes(:,1) ) ) ;
        CumSumSpk = [0 cumsum(nbSpk)] ;
        
        SpkTimes = {} ;
        for i=1:length(CumSumSpk)-1
            fprintf('# %d nIdx %d nbSpk %.3f ',nidx(i),reSizeIdx(i),nbSpk(i)) 
            try
                SpkTimes = [ SpkTimes;{Spikes(CumSumSpk(i)+1:CumSumSpk(i)+nbSpk(i)-1,2).'} ] ; 
                ISI = SpkTimes{i}(2:nbSpk(i)-1)-SpkTimes{i}(1:nbSpk(i)-2) ; 
                CV(i) = CVfunc(ISI) ;
                CV2(i) = CV2func(ISI) ;
            catch
                CV(i) = nan ;
                CV2(i) = nan ;
            end
            fprintf('CV %.3f CV2 %.3f',CV(i),CV2(i))
            fprintf('\r')
        end
        fprintf('\n')
        
        popList = ['E' 'I' 'S' 'X'] ;
        figname=sprintf('Raster_Iext%.2f',Iext(nPrtr)) ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

        for j=1:nbpop
            for i=1:length(CumSumSpk)-1
                if(Cpt(j)<=nidx(i) & nidx(i)<Cpt(j+1))
                    % fprintf('%s # %d Cpt %d \r',popList(j),nidx(i),Cpt(j+1))
                    % plot(SpkTimes{i},nidx(i).*ones(nbSpk(i)-1,1),'.','Color',cl{j},'MarkerSize',.001) 
                    plot(SpkTimes{i},reSizeIdx(i).*ones(nbSpk(i)-1,1),'.','Color',cl{j},'MarkerSize',.001) 
                end
            end
        end
        % fprintf('\n')

        xlabel('t (s)')
        ylabel('#')
        xlim([0 .5])
        drawnow ;

        hold off ;

        if(IF_SAVE)
            figdir = sprintf(['./Figs/IF/%dpop/%s/'],nbpop,dir);
            if(IF_RING)
                figdir = sprintf(['./Figs/IF/%dpop/%s/Crec%.2fCff%.2f'],nbpop,dir,Crec,Cff);
            end
            try
                mkdir(figdir)
            end
            
            ProcessFigure(fig, fullfile(figdir,figname)) ;
        end        
        
        % figname=sprintf('CVDist_Iext%.2f',Iext(nPrtr)) ;
        % fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        % for i=1:nbpop
        %     h = histogram(CV(Cpt(i)<nidx & nidx<=Cpt(i+1)),100,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i}) ;
        % end
        % xlabel('CV')
        % ylabel('pdf')
        % hold off ;
        
        % if(IF_SAVE)
        %     figdir = sprintf(['./Figs/IF/%dpop/%s/'],nbpop,dir);
        %     if(IF_RING)
        %         figdir = sprintf(['./Figs/IF/%dpop/%s/Crec%.2fCff%.2f'],nbpop,dir,Crec,Cff);
        %     end
        %     try
        %         mkdir(figdir)
        %     end
            
        %     ProcessFigure(fig, fullfile(figdir,figname),1.25) ;
        % end        

        % figname=sprintf('CV2Dist_Iext%.2f',Iext(nPrtr)) ;
        % fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        % for i=1:nbpop
        %     h = histogram(CV2(Cpt(i)<nidx & nidx<=Cpt(i+1)),100,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i}) ;
        % end
        % xlabel('CV2')
        % ylabel('pdf')
        % hold off ;

        % if(IF_SAVE)
        %     figdir = sprintf(['./Figs/IF/%dpop/%s/'],nbpop,dir);
        %     if(IF_RING)
        %         figdir = sprintf(['./Figs/IF/%dpop/%s/Crec%.2fCff%.2f'],nbpop,dir,Crec,Cff);
        %     end
        %     try
        %         mkdir(figdir)
        %     end
            
        %     ProcessFigure(fig, fullfile(figdir,figname),1.25) ;
        % end        

    end
    
    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%
    % Membrane Voltage %
    %%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%

    if strcmp(file,'Voltage')
        try
            data = ImportData(model,nbpop,dir,'Voltage',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ; 
            rates = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ; 
            for i=1:length(data(1,:))-1
                IdvRates(i) = mean(rates(:,i+1)) ;
            end
        end
        tps = data(:,1)./1000 ;
        Volt = data(:,2:end) ;
        nId = randi([1 10]) ;
        figname=sprintf('#%d m_i=%.3f',nId, IdvRates(nId) ) ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        plot(tps,Volt(:,nId),'LineWidth',.25)
        xlabel('t (s)')
        ylabel('Vm')
        hold off ;
    end

    %%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%
    % Input Distribution %
    %%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%

    if strcmp(file,'IdvInputs')
        try
            data = ImportData(model,nbpop,dir,file,n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ; 
            for i=1:length(data(1,:))-1
                IdvInputs(i) = mean(data(:,i+1)) ;
            end
            tps = data(:,1)./1000 ;
            for i=1:nbpop
                for j=1:length(data(:,1))
                    PopInput(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1))) ;
                end
            end
        catch
            if(strfind(dir,'2pop'))
                for i=1:n*2/3*10000
                    IdvInputs(i) = nan ;
                end
            else
                for i=1:n*10000
                    IdvInputs(i) = nan ;
                end
            end
            tps = [] ;
            PopInput = [] ;
        end
        
        figname=sprintf('InputDist') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

        for i=1:nbpop
            MeanInput(i) = mean( IdvInputs( Cpt(i)+1:Cpt(i+1) ) ) ;
            VarInput(i) = var( IdvInputs( Cpt(i)+1:Cpt(i+1) ) ) ;            
            h = histogram(log(IdvInputs( Cpt(i)+1:Cpt(i+1) )),100, ...
                          'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i}) ;
            % set(gca, 'XScale', 'log')
        end

        xlabel('log(I_{syn})')
        ylabel('pdf')

        hold off ;

        fprintf('Simuls : ')
        fprintf('u ')
        fprintf('%.3f | ', MeanInput)
        fprintf('b ')
        fprintf('%.3f | ', VarInput)
        fprintf('\n')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inset PopAvgInput evolution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % size(data)
        % length(tps)
        % size(PopInput)
        
        axes('Position',[.25 .7 .2 .2]) ; hold on ;
        box off
        for i=1:nbpop
            plot(tps,PopInput(i,:),'color',cl{i})
        end    
        xlabel('t (s)')
        ylabel('<I_{syn}>_i')
        set(gca,'FontSize',8)
        hold off ;

        axes('Position',[.25 .35 .2 .2]) ; hold on ;
        box off
        for i=1:nbpop
            for k=1:nb
                nId = randi([Cpt(i)+1 Cpt(i+1)]) ;
                patchline(tps,data(:,nId),'linestyle','-','edgecolor',cl{i},'edgealpha',.25) 
            end
        end
        xlabel('t (s)')
        ylabel('I_{syn_i}')
        set(gca,'FontSize',8)

        if(IF_SAVE)
            figdir = sprintf(['./Figs/IF/%dpop/%s/'],nbpop,dir);
            if(IF_RING)
                figdir = sprintf(['./Figs/IF/%dpop/%s/Crec%.2fCff%.2f'],nbpop,dir,Crec,Cff);
            end

            try
                mkdir(figdir)
            end
            
            ProcessFigure(fig, fullfile(figdir,figname),1.25) ;
        end        

        hold off ;

    end

    %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%
    % Rate Distribution %
    %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%

    if strcmp(file,'IdvRates')
        try
            data = ImportData(model,nbpop,dir,file,n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ; 
            for i=1:length(data(1,:))-1
                IdvRates(i) = mean(data(:,i+1)) ;
            end
            tps = data(:,1)./1000 ;
            for i=1:nbpop
                for j=1:length(data(:,1))
                    PopRate(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1))) ;
                end
            end
        catch
            if(strfind(dir,'2pop'))
                for i=1:n*2/3*10000
                    IdvRates(i) = nan ;
                end
            else
                for i=1:n*10000
                    IdvRates(i) = nan ;
                end
            end
            tps = [] ;
            PopRate = [] ;
        end
        
        figname=sprintf('RateDist') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

        for i=1:nbpop
            MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
            VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ; 
            m = IdvRates( Cpt(i)+1:Cpt(i+1) ) ; 
            m = m(m>THRESHOLD) ; 
            h = histogram(log(m),100,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i}) ; 
            % set(gca, 'XScale', 'log')
        end

        xlabel('log(m)')
        ylabel('\rho(m)')

        hold off ;

        MF = BalRatesMF(model,nbpop,dir,Iext,[]) ;
        fprintf('MF : ')
        fprintf('%.3f | ', MF)
        fprintf('\n')

        fprintf('Simuls : ')
        fprintf('%.3f | ', MeanRate)
        fprintf('\n')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inset PopAvgRate evolution %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % size(data)
        % length(tps)
        % size(PopInput)

        axes('Position',[.7 .7 .2 .2]) ; hold on ;
        box off
        for i=1:nbpop
            plot(tps,PopRate(i,:),'color',cl{i})
        end    
        xlabel('t (s)')
        ylabel('<m>_i')
        set(gca,'FontSize',8)

        for i=1:nbpop
            for k=1:nb
                nId = randi([Cpt(i)+1 Cpt(i+1)]) ;
                patchline(tps,data(:,nId),'linestyle','-','edgecolor',cl{i},'edgealpha',.25) 
            end
        end
        xlabel('t (s)')
        ylabel('m_i')
        set(gca,'FontSize',8)
        hold off ;

        axes('Position',[.3 .7 .2 .2]) ; hold on ;
        box off
         
        for i=1:nbpop        
            m = IdvRates(1,Cpt(i)+1:Cpt(i+1)) ;
            idx = find(m>THRESHOLD) ;
            m = m(idx) ;
            
            patchline(idx,m,'linestyle','-','edgecolor',cl{i},'edgealpha',.25) 
            fit = smooth(idx,m,.01,'sgolay') ;
            plot(idx,fit,'linewidth',.2,'color','k')
        end
        
        xlabel('#')
        ylabel('[m_i]_t')
        set(gca,'FontSize',8)
        hold off ;

        if(IF_SAVE)
            figdir = sprintf(['./Figs/IF/%dpop/%s/'],nbpop,dir);
            if(IF_RING)
                figdir = sprintf(['./Figs/IF/%dpop/%s/Crec%.2fCff%.2f'],nbpop,dir,Crec,Cff);
            end
            try
                mkdir(figdir)
            end
            
            ProcessFigure(fig, fullfile(figdir,figname),1.25) ;
        end

    end    

    if strcmp(file,'Scatter')
        try
            data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ; 
            for i=1:length(data(1,:))-1
                IdvRates(i) = mean(data(:,i+1)) ;
            end
        catch
            for i=1:n*10000
                IdvRates(i) = nan ;
            end
        end

        try
            Prtr = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr+.4) ; 
            for i=1:length(data(1,:))-1
                IdvRatesPrtr(i) = mean(Prtr(:,i+1)) ;
            end
        catch
            for i=1:n*10000
                IdvRatesPrtr(i) = nan ;
            end
        end

        popList = ['E' 'I' 'S' 'X'] ;
        
        if( strcmp(IF_RING,'Ring') | strcmp(IF_RING,'Gauss') | strcmp(IF_RING,'Exp') )
            Cth = .125 ;
        else
            Cth = 100 ;
        end

        for i=1:nbpop
            
            figname=sprintf('Scatter%s',popList(i)) ;
            fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

            MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
            VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ; 

            m = IdvRates( Cpt(i)+1:Cpt(i+1) ) ;
            idx = find(m>=THRESHOLD) ;
            m = m(idx) ; 

            mPrtr = IdvRatesPrtr( Cpt(i)+1:Cpt(i+1) ) ; 
            mPrtr = mPrtr(idx) ;

            X = linspace(-pi,pi,length( IdvRates(Cpt(i)+1:Cpt(i+1) ) ) ) ;
            X = X(idx) ;

            idx = find(abs(X)<=Cth) ; % Cth
            X = X(idx) ;
            m = m(idx) ; 
            mPrtr = mPrtr(idx) ;
            
            sc = scatter(m,mPrtr,.1,cl{i},'MarkerEdgeAlpha',.5) ;
            set(sc,'SizeData',.1);

            xLimits = get(gca,'XLim');
            yLimits = get(gca,'YLim');
            plot([0 xLimits(2)],[0 xLimits(2)],'--k','Linewidth',.5) 
            ylim(yLimits)

            xlabel('Baseline')
            ylabel('Perturbartion')


            if(IF_SAVE)
                figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff) ;
                fprintf('Writing %s \n',figdir)
                try
                    mkdir(figdir)
                end
                
                ProcessFigure(fig, fullfile(figdir,figname),1.25) ;
            end


            hold off ;
        end


        MF = BalRatesMF(model,nbpop,dir,Iext,[]) ;
        fprintf('MF : ')
        fprintf('%.3f | ', MF)
        fprintf('\n')

        fprintf('Simuls : ')
        fprintf('%.3f | ', MeanRate)
        fprintf('\n')

        hold off ;

    end    

    %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%
    % Spatial Profile %
    %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%

    if strcmp(file,'Space')

        ErrFlag = 0 ;

        try
            data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ;
            DataLength = length(data(1,:))-1 ;

            if(IF_IEXT)
                IextBL = ExternalInput(model,nbpop,dir) ;
                
                Baseline = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,IextBL(nPrtr)) ;
                for i=1:DataLength
                    IdvRatesBL(i) = mean(Baseline(:,i+1)) ;
                end
            else 
                for i=1:DataLength
                    IdvRatesBL(i) = 1 ;
                end
            end

            for i=1:DataLength
                IdvRates(i) = mean(data(:,i+1)) ;
            end

        catch
            fprintf('ERROR FILE NOT FOUND\n')
            for i=1:DataLength
                IdvRates(i) = nan ;
            end
            ErrFlag = 1 ;
        end
        
        if ErrFlag == 0
            for i=1:nbpop
                if(i==1 | i==2)
                    figname=sprintf('%s_Spatial_EI',dir) ;
                else
                    figname=sprintf('%s_Spatial_SV',dir) ;
                end
                if( ishandle( findobj('type','figure','name',figname) ) )
                    fig = findobj('type','figure','name',figname) ; 
                    figure(fig); hold on ; 
                else
                    fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
                    xlabel('x (mm)')
                    ylabel('Norm. Rates')
                end
                
                popLength = length(Cpt(i)+1:Cpt(i+1)) ;

                if(IF_IEXT)
                    BL = IdvRatesBL( Cpt(i)+1:Cpt(i+1) ) ;
                    leftBL = fliplr( BL(1:popLength/2) ) ;
                    rightBL = BL(popLength/2+1:end) ;

                    mBL = leftBL + rightBL ;

                    idx = find(mBL>=THRESHOLD) ;
                    mBL = mBL(idx) ;
                end
                
                m = IdvRates( Cpt(i)+1:Cpt(i+1) ) ;
                MeanRate(i) = mean(m) ;
                VarRate(i) = var(m) ; 
                
                mLeft = fliplr( m(1:popLength/2) );
                mRight = m(popLength/2+1:end) ;
                mTot = mLeft + mRight ;

                MeanRate(i) = mean(mTot) ;
                VarRate(i) = var(mTot) ; 

                X = linspace(0,pi,popLength/2) ; 

                % length(mLeft) 
                % length(mRight)
                % length(X) 

                % mTot = mTot(idx)./mBL ;
                X = X(idx) ;            

                idxMax = max( find(X<=Cff/2) ) ;

                fit = smooth(idx(1:idxMax),mTot(1:idxMax),.001,'lowess') ;
                patchline(X(1:idxMax),fit,'linestyle','-','edgecolor',cl{i},'edgealpha',.1)

                fit = smooth(idx(idxMax+1:end),mTot(idxMax+1:end),.001,'lowess') ;
                patchline(X(idxMax+1:end),fit,'linestyle','-','edgecolor',cl{i},'edgealpha',.1)
                
                fit = smooth(idx(1:idxMax),mTot(1:idxMax),.05,'lowess') ;
                plot(X(1:idxMax),fit,'linewidth',2,'color',cl{i})

                fit = smooth(idx(idxMax+1:end),mTot(idxMax+1:end),.05,'lowess') ;
                plot(X(idxMax+1:end),fit,'linewidth',2,'color',cl{i})
                
                % % X = linspace(-pi,pi,length( IdvRates(Cpt(i)+1:Cpt(i+1) ) ) ) ; 
                % % idx = find( abs(X)<=pi ) ;
                % % X = X(idx) ;
                % % fprintf(' %d %d \n',idx(1),idx(length(idx)))
                
                % % m = IdvRates(idx(1)+Cpt(i):idx(end)+Cpt(i) ) ;
                % % fprintf('Rates %.3f \n',mean(m))
                            
                % X = linspace(-pi,pi,nbN(i)) ;
                % X = X(idx) ;
                
                % fit = smooth(idx,m,.001,'lowess') ;
                % patchline(X,fit,'linestyle','-','edgecolor',cl{i},'edgealpha',.1) 
                
                % fit = smooth(idx,m,.05,'lowess') ;
                % plot(X,fit,'linewidth',2,'color',cl{i})
                % b = fitsine(X,fit.') ;
                % Sol = @(x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);
                % plot(X,Sol(X),'--k')
                
                % patchline(X(idx),m,'linestyle','-','edgecolor',cl{i},'edgealpha',.25) 
                % fit = smooth(idx,m,.1,'sgolay') ;
                % plot(X(idx),fit,'linewidth',.2,'color',cl{i})
                xlim([0 2])
                ylim([0 2.5])
                
                plot([4*Crec(1) 4*Crec(1)],[0 3],'--r', 'linewidth', 1)
                plot([4*Crec(2) 4*Crec(2)],[0 3],'--b', 'linewidth', 1)
                plot([Cff/2 Cff/2],[0 3],'--k', 'linewidth', 1)
                
                if(IF_SAVE & (i==2 | i==4))
                    
                    figdir = FigDir(model,nbpop,dir,n,K,g,IF_RING,Crec,Cff) ;
                    fprintf('Writing %s \n',figdir)
                    
                    try
                        mkdir(figdir) ;
                    end
                    
                    ProcessFigure(fig, fullfile(figdir,figname)) ;
                end
            end
              
            MF = BalRatesMF(model,nbpop,dir,Iext,[]) ;
            fprintf('MF : ')
            fprintf('%.3f | ', MF)
            fprintf('\n')
            
            fprintf('Simuls : ')
            fprintf('%.3f | ', MeanRate)
            fprintf('\n')
            
        end    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(file,'MeanRates')
        try
            data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ; 
            for i=2:length(data(1,:))
                IdvRates(i) = mean(data(:,i)) ;
            end
            tps = data(:,1) ;
            for i=1:nbpop
                for j=1:length(data(:,1))
                    PopRate(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1)+1)) ;
                end
            end
        catch
            fprintf(' FILE NOT FOUND \n')
            if(strfind(dir,'2pop'))
                for i=1:n*2/3*10000
                    IdvRates(i) = nan ;
                end
            else
                for i=1:n*10000
                    IdvRates(i) = nan ;
                end
            end
            tps = [] ;
            PopRate = [] ;
        end
        
        figname=sprintf('Mean Rates') ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 

        for i=1:nbpop
            MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1)+1 ) ) ;
            VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1)+1 ) ) ;
            plot(tps,PopRate(i,:),'color',cl{i})
        end   
        xlabel('t (ms)')
        ylabel('<m>_i')

        drawnow ;
        hold off ;

        fprintf('Simuls : ')
        fprintf('%.3f | ', MeanRate)
        fprintf('\n')
    end

    %%%%%%%%%%%%%%%%%%%%%
    if strcmpi(file,'Bump')

        % set(0, 'DefaultFigureRenderer', 'zbuffer') ;
        figname = 'Bump' ;
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ;
        xlabel('x')
        ylabel('[m_i]_x') 
        drawnow ;
        axis tight manual 
        set(gca,'nextplot','replacechildren') ; 
        
        %writerObj = VideoWriter('Bump.avi');
        writerObj = VideoWriter('CosBumpEI.avi');
        open(writerObj);

        try
            data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ;
        end

        binLength = length(data(:,1)) ;
        neuronLength = length(data(1,:))-1 ;
        phsList  = [] ;
        fprintf('binLength %d, nbNeuron %d \n',binLength,neuronLength) 

        for k=1:binLength

            fprintf('bin %d ',k)
            fprintf('Rates ')
            for i=1:nbpop
                fprintf('%.3f ', mean( data(k,Cpt(i)+2:Cpt(i+1)+1) ) )
            end

            fprintf('Phase ')
            for i=1:nbpop

                m = data(k,Cpt(i)+2:Cpt(i+1)+1) ;
                idx = find(m>=THRESHOLD) ;
                m = m(idx) ;
                
                X = linspace(-pi/2,pi/2,nbN(i)) ;
                Z = exp(2*1j*X) ;
                X = X(idx) ;
                
                % fit = smooth(idx,m,.001,'lowess') ;
                % patchline(X,fit,'linestyle','-','edgecolor',cl{i},'edgealpha',.1) 

                fit = smooth(idx,m,.05,'lowess') ;
                plot(X,fit,'linewidth',2,'color',cl{i}) ; hold on ;
                xlim([-pi/2 pi/2])
                ylim([0 15])
                xlabel('x')
                ylabel('[m_i]_x') 

                args = fitcos(X,fit.') ;
                Sol = @(x)  args(1).*(cos(2*pi*x./args(2) + 2*pi/args(3))) + args(4) ;
                plot(X,Sol(X),'--','color',cl{i})

                phs0 =  wrapToPi( angle( dot(Z(idx),m) ) ) ;
                phs = wrapToPi(pi/args(3)) ;
                % plot([phs phs],[0 20],'--k')
                fprintf('%.3f %.3f ',phs,phs0)
                phsList = [phsList phs0] ;
                drawnow ;
            end
            
            hold off ;
            fprintf('\r')

            frame=getframe(gcf);
            writeVideo(writerObj,frame) ;
            % clf('reset') ;

        end

        close(writerObj);
        
        figure()
        plot(data(:,1)/1000,phsList)
        xlabel('x')
        ylabel('\phi_0')

        fprintf('\n')

    end
    %%%%%%%%%%%%%%%%%%%%%

    if strcmpi(file,'Phase')

        % set(0, 'DefaultFigureRenderer', 'zbuffer') ;
        figname = 'Phase' ;
        fig = figure('Name',figname,'NumberTitle','off') ;

        try
            data = ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr) ;
        end

        binLength = length(data(:,1)) ;
        neuronLength = length(data(1,:))-1 ;
        phsList  = [] ;
        fprintf('binLength %d, nbNeuron %d \n',binLength,neuronLength) 

        IF_PHI0 = 1 ;
        for k=1:binLength

            fprintf('bin %d ',k)
            fprintf('Rates ')
            for i=1:nbpop
                fprintf('%.3f ', mean( data(k,Cpt(i)+2:Cpt(i+1)+1) ) )
            end
            
            fprintf('Phase ')
            for i=1:1

                m = data(k,Cpt(i)+2:Cpt(i+1)+1) ;
                idx = find(m>=THRESHOLD) ;
                m = m(idx) ;
                
                X = linspace(-pi/2,pi/2,nbN(i)) ;
                Z = exp(2*1j*X) ;
                X = X(idx) ;
                
                fit = smooth(idx,m,.05,'lowess') ;

                args = fitcos(X,fit.') ;
                Sol = @(x)  args(1).*(cos(2*pi*x./args(2) + 2*pi/args(3))) + args(4) ;

                phs0 = abs( wrapToPi( 2 * angle(dot(Z(idx),m)) ) ) -pi/2 ;

                if(IF_PHI0)
                    IF_PHI0 = 0 ;
                    PHI0 = phs0 ;
                end

                % phs0 = phs0 - PHI0 ;
                phs = wrapToPi(2*pi/args(3)) ;

                fprintf('%.3f %.3f ',phs,phs0)
                phsList = [phsList phs0] ;
            end
            
            fprintf('\r')
        end
        
        plot(data(:,1)./1000,phsList)
        xlabel('t (s)')
        ylabel('\phi_0')

        fprintf('\n')

    end

    %%%%%%%%%%%%%%%%%%%%%
    if strcmpi(file,'Bump2D')
        
        try
            data=ImportData(model,nbpop,dir,'IdvRates',n,K,g,IF_RING,Crec,Cff,IF_IEXT,nPrtr,Iprtr);
        end
        
        binLength = length(data(:,1)) ;
        neuronLength = length(data(1,:))-1 ;
        
        figname = 'Bump2D' ;
        fig = figure('Name',figname,'NumberTitle','off') ;
        set(gca,'nextplot','replacechildren') ; 
        
        for k=1:binLength        
            for i=1:1
                
                Rates = data(k,Cpt(i)+2:Cpt(i+1)+1) ;
                Rates = reshape( Rates, sqrt( length( Rates ) ), sqrt( length( Rates ) ) );

                M = interp2(Rates,'cubic') ;                                    
                imagesc(M) ;

                h = colorbar ;
                caxis([0 10])
                xlabel('\theta_x ')
                ylabel('\theta_y ')
                
                xlim([0 length(M(:,1))])
                ylim([0 length(M(:,1))])

                drawnow ;
            end
            
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%
    %% Utils             %
    %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%

    function out = CVfunc(ISI)
        out = sqrt(var(ISI))./mean(ISI) ;
    end

    function out = CV2func(ISI)
        dISI = abs(ISI(1:end-1)-ISI(2:end)) ;
        sumISI = ISI(1:end-1)+ISI(2:end) ;
        out = 2.* mean( dISI )./mean( sumISI ) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = SmoothSignal(x)
        windowSize = 5 ; 
        Tw = (1./windowSize).*ones(1,windowSize) ; 
        % Tw = ones(1,length(x))/length(x) ;
        out = filter(Tw,.08,x) ;
    end

    function out = M1Component(x)
        dPhi = pi./length(x) ;
        y = 1:length(x) ;
        out = 2.*abs(x*exp(-2j.*y.*dPhi).' )./length(x) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    function out = NthFourierMoment(x, n, phaseShift)
        out = [] ;
        nPoints = length(x) ;
        if(n==0)
            preFactor = 1 ;
        else
            preFactor = 2 * n ;
        end
        y = linspace(0, pi, nPoints) - phaseShift ;
        out = preFactor.*x*cos(2.*n.*y).'./nPoints ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%

end
