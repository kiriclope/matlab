clear all ;
GlobalVars

Iext = ExternalInput(model,nbpop,dir) ;    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

Iext(prtrPop) = Iext(prtrPop) + Iprtr ; 

try
    data = ImportData(model,nbpop,dir,'IdvRates',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iext(prtrPop)) ; 
catch
    return ;
end

for i=1:length(data(1,:))-1
    IdvRates(i) = mean(data(:,i+1)) ;
end
tps = data(:,1)./1000 ;
for i=1:nbpop
    for j=1:length(data(:,1))
        PopRate(i,j) = mean(data(j,Cpt(i)+2:Cpt(i+1))) ; 
    end
end

for i=1:nbpop
    
    switch i
      case 1
        figname='RatesDistE' ;
      case 2
        figname='RatesDistI' ;
      case 3
        figname='RatesDistS' ;
      case 4
        figname='RatesDistV' ;
    end
            
    if( ishandle( findobj('type','figure','name',figname) ) )
        fig = findobj('type','figure','name',figname) ; 
        figure(fig); hold on ; 
    else
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        xlabel('log(Rates)')
        ylabel('pdf')
    end
    
    alp = 1 ; 
    
    MeanRate(i) = mean( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ;
    VarRate(i) = var( IdvRates( Cpt(i)+1:Cpt(i+1) ) ) ; 
    m = IdvRates( Cpt(i)+1:Cpt(i+1) ) ; 

    m = m(m>THRESHOLD) ; 
    
    if(strfind(model,'Binary')) 
        h = histogram(m,27,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i},'EdgeAlpha',.2,'Linewidth',2) ;
        if(i==1)
            [u a b] = Bin_InputDist(nbpop,dir,Iext*m0,K,1,0) ; 
            x =.01:.01:1;
            rho = @(m,u,a,b) sqrt((a-b)./b).*exp( -(1-u-sqrt(a-b).*sqrt(2).*erfcinv(2.*m) ).^2./(2.*b) + erfcinv(2.*m).^2 ) ; 
        end 
        patchline(x,rho(x,u(i),a(i),b(i)),'linestyle','-','edgecolor',cl{i},'edgealpha',1) 
        xlabel('m')
        ylabel('\rho(m)') 
        % set(gca,'xscale','log') 
    else
        h = histogram(log(m)/log(10),27,'Normalization', 'pdf' ,'DisplayStyle','stairs','EdgeColor',cl{i},'EdgeAlpha',alp,'Linewidth',2) ;
        xlim([-2 2])
        set(gca,'xtick',[-2 -1 0 1 2],'xticklabel',{'10^{-2}','10^{-1}','10^0','10^1','10^2'})
    end

    drawnow;
    if(IF_SAVE)
        figdir = FigDir(model,nbpop,dir,N,K,g,IF_RING,Crec,Cff,IF_IEXT) ;
        fprintf('Writing %s \n',figdir)
        try 
            mkdir(figdir) 
        end 
        ProcessFigure(fig, fullfile(figdir,figname), 2.2, [1.33*2.2, 2.2]) ;
    end        
    hold off ;
end


MF = BalRatesMF(model,nbpop,dir,Iext,[]) ;
fprintf('MF : ')
fprintf('%.3f | ', MF)
fprintf('\n')

if(strfind(model,'Binary')) 
    QchAvgTF = @(u,a) .5*erfc( ( 1-u )./sqrt(2.*a ) ) ; 
    fprintf('Numerics : ')
    fprintf('%.3f | ', QchAvgTF(u,a) )
    fprintf('\n')
end

fprintf('Simuls : ')
fprintf('%.3f | ', MeanRate)
fprintf('\n')
