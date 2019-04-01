function fig = popPerFig(pop,dir,figtitle) 
    if(pop==1 || pop ==2)
        figname=sprintf('%s_%s_EI',dir,figtitle) ;
    else
        figname=sprintf('%s_%s_SV',dir,figtitle) ;
    end

    if( ishandle( findobj('type','figure','name',figname) ) )
        fig = findobj('type','figure','name',figname) ; 
        fig = figure(fig); hold on ; 
    else
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
    end
end
