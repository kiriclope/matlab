clear all ; 
GlobalVars 

Iext = ExternalInput(model,nbpop,dir) ;    
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

% Iext(prtrPop) = Iext(prtrPop)+Iprtr ;

try
    data = ImportData(model,nbpop,dir,'Raster',N,K,g,IF_RING,Crec,Cff,IF_IEXT,prtrPop,Iprtr) ; 
    Spikes = sortrows(data) ; 
    Spikes(:,2) = Spikes(:,2)./1000. ;
catch
    fprintf('FILE NOT FOUND\n') ;
    return ;
end

Tw = .5 ;
WinSize = ceil( Spikes(end,2) ) / Tw ;
fprintf('Duration %ds Tw %ds \n', WinSize*Tw, Tw) 

reSizeIdx = 1:length( unique( Spikes ) ) ;
[nbSpk nidx]= hist(Spikes(:,1),unique( Spikes(:,1) ) ) ; 
CumSumSpk = [0 cumsum(nbSpk)] ; 

Fano = NaN(Cpt(nbpop+1),1) ;

fprintf('nbN ')
for j=1:nbpop
    fprintf('%d ', nbN(j))
end 

fprintf('\nlength ROI ')
for j=1:nbpop
    X=[] ; 
    Y=[] ; 
    for i=1:nbN(j)
        X(i) = L * mod( double(i), sqrt( double( nbN(j) ) ) ) / sqrt( double( nbN(j) ) ) ; 
        Y(i) = L * floor( double(i) / sqrt( double( nbN(j) ) ) ) / sqrt( double( nbN(j) ) ) ; 
    end    
    ROI{j} = find( (X-L/2).^2 + (Y-L/2).^2 <= Cth.^2 / 4 ) ; 
    fprintf('%d ',length(ROI{j}))
end
fprintf('\n')

SpkTimes = {} ; 
for i=1:length(CumSumSpk)-1
    fprintf('# %d nIdx %d nbSpk %d ', nidx(i), reSizeIdx(i), nbSpk(i)) 
       
    SpkTimes = Spikes(CumSumSpk(i)+1:CumSumSpk(i+1),2).' ; 
    
    for j=1:WinSize
        SpkCount(i,j) = length( intersect( find( SpkTimes<j*Tw ), find( SpkTimes>=(j-1)*Tw ) ) ) ;
    end    
    
    for j=1:nbpop
        if( ismember(nidx(i)-Cpt(j),ROI{j}) ) 
            Fano( nidx(i)+1 ) = var( SpkCount(i,:) ) / mean( SpkCount(i,:) )  ;
        end
    end

    fprintf('Fano %.3f ', Fano(nidx(i)+1) ) 
    fprintf('\r')

end


fprintf('\n')
fprintf('mean Fano ') 

for i=1:nbpop
    switch i 
      case 1
        figname='FanoDistE' ;
      case 2
        figname='FanoDistI' ;
      case 3
        figname='FanoDistS' ;
      case 4
        figname='FanoDistV' ;
    end
    
    if( ishandle( findobj('type','figure','name',figname) ) ) 
        fig = findobj('type','figure','name',figname) ; 
        figure(fig); hold on ; 
    else
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        xlabel('Fano Factor') 
        ylabel('Probability') 
    end
    
    FF = Fano( Cpt(i)+1:Cpt(i+1) ) ; 
    FF(isnan(FF)) = [] ;

    fprintf('%.2f ', mean( FF ) ) 
    h = histogram( FF , 30, 'Normalization', 'probability' ,'DisplayStyle','stairs','EdgeColor',cl{i},'EdgeAlpha',alp,'Linewidth',2) ; 
    xlim([0 3]) 
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
fprintf('\n')