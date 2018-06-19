function ProcessFigure(hFigure, fileName,h,paperSize)
% ProcessFigure(hFigure, fileName, varargin)
% [fontSize, paperSize]
% History:
% Shrisha - Created
   
    fontSize = 12 ;
    labelFntSize = 10 ;
    lw=1.5 ;
    msz = 8 ;

    if(nargin<3)
        h = 3 ;
    else
        fontSize = 12 ;
        labelFntSize = 10 ;
        lw=.0001/3 ;
        msz = .05 ;
    end
    
    if(nargin<4)
        paperSize = [1.33*h, h] ;
    end

    try
        allAxesInFigure = findall(hFigure, 'type', 'axes') ; 
        
        for hAxis = allAxesInFigure'
            set(hAxis, 'box', 'off');
            set(hAxis, 'TickDir', 'out');
            set(hAxis,'FontSize', fontSize);
            set(get(hAxis,'XLabel'), 'FontSize', fontSize);
            set(get(hAxis,'YLabel'), 'FontSize', fontSize);
            txtHdl = findall(hAxis, 'type', 'text');
            for kTxtHdl = txtHdl
                set(kTxtHdl, 'fontSize', fontSize);
            end
        end
    
        set(hFigure, 'PaperPosition', [0, 0, paperSize]);
        set(hFigure, 'PaperSize', paperSize);
        set(hFigure, 'MarkerSize', msz);
        %    set(hFigure, 'Renderer', 'Painters');
    end    
    
    saveas(hFigure, [fileName '.pdf'], 'pdf');
    saveas(hFigure, [fileName '.fig'], 'fig');
    saveas(hFigure, [fileName '.svg'], 'svg');
    % saveas(hFigure, [fileName '.png'], 'png');

end