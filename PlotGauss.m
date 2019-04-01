function PlotGauss(Cff)    
    
    L = 2 ;
    x = -L/2:.01:L/2 ;

    figname = 'DeltaGauss';
    if( ishandle( findobj('type','figure','name',figname) ) )
        fig = findobj('type','figure','name',figname) ; 
        figure(fig); hold on ; 
    else
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        xlabel('x (mm)')
        ylabel('I_{opto}')
    end

    plot(x, Gauss(CircDist(x), Cff),'b--')  

    figname = 'DeltaGauss';
    if( ishandle( findobj('type','figure','name',figname) ) )
        fig = findobj('type','figure','name',figname) ; 
        figure(fig); hold on ; 
    else
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        xlabel('x (mm)')
        ylabel('I_{opto}')
    end

    figname = 'DeltaGauss';
    if( ishandle( findobj('type','figure','name',figname) ) )
        fig = findobj('type','figure','name',figname) ; 
        figure(fig); hold on ; 
    else
        fig = figure('Name',figname,'NumberTitle','off') ; hold on ; 
        xlabel('x (mm)')
        ylabel('I_{opto}')
    end

    plot(x, DeltaGauss(CircDist(x), Cff),'b') 
    plot([Cff Cff], [0 exp(-1/2)/sqrt(2*pi)/Cff], 'k--')  
    plot([-Cff -Cff], [0 exp(-1/2)/sqrt(2*pi)/Cff], 'k--') 
    plot([4*Cff 4*Cff], [0 1/sqrt(2*pi)/Cff], 'k--') 
    plot([-4*Cff -4*Cff], [0 1/sqrt(2*pi)/Cff], 'k--') 
    xlim([-L/2 L/2])

    function out = DeltaGauss(X,Sigma)
        out = [] ;
        
        Norm = exp( - 1 / 2 ) / sqrt(2 * pi ) / Sigma ;
        for i=1:length(X)
            if( abs(X(i)) <= Sigma )
                out = [out exp( - 1 / 2 ) / sqrt(2 * pi ) / Sigma / Norm ] ; 
            else
                if(abs(X(i)) <= 4 * Sigma )
                    out = [out exp( - X(i) * X(i) / 2 /Sigma /Sigma ) / sqrt(2 * pi ) / Sigma /Norm ] ;
                else
                    out = [out nan] ;
                end 
            end
        end
    end
    
    function out = Gauss(X,Sigma)
        out = [] ;
        Norm = exp( - 1 / 2 ) / sqrt(2 * pi ) / Sigma ;
        for i=1:length(X)
            if(abs(X(i)) <= 4 * Sigma )            
                out = [out exp( - X(i) * X(i) / 2 / Sigma /Sigma ) / sqrt(2 * pi ) / Sigma /Norm ] ; 
            else
                out = [out nan] ;
            end 
        end
    end

    function out = CircDist(X) 
        out = [] ; 
        for i=1:length(X) 
            if( mod(abs(X(i)),L) > .5*L ) 
                out = [ out, mod(abs(X(i)),L)-L ] ;
            else  
                out = [ out, L/2*mod(abs(X(i)),L) ] ;
            end 
        end                                 
    end

end