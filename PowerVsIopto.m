clear all ;
GlobalVars

v_Iprtr = v_Iprtr(1):v_Iprtr(2):v_Iprtr(3) ;

IF_POWER=2 ;

if IF_POWER==1
    Popto = @(X) Pinf ./ ( 1 + exp(-(X - I0 )/Iinf ) ) ;
    Iopto = @(Y) I0./20 .* ( log( 1 + Y / P0 ) ) ; 
elseif IF_POWER==2
    Popto = @(X) -P0 * log( 1 - X ./ I0 ) ;
    Iopto = @(Y)  I0 .* (1 - exp( - Y ./ P0 ) ) ;  
    Popto = @(X) P0 .* ( exp( 20*X ./ I0 ) - 1 ) ; 
    Iopto = @(Y) I0/5 .* ( log( 1 + Y / P0 ) ) ; 
end

fprintf('%s Iprtr %.3f Pprtr %.3f\n', dir, prtrAmp, Popto( prtrAmp ) )

% figtitle = sprintf('PowerVsIopto') ;
% if( ishandle( findobj('type','figure','name',figtitle) ) )
%     fig = findobj('type','figure','name',figtitle) ; 
%     fig = figure(fig); hold on ; 
% else
%     fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
%     xlabel('\Gamma_{opto} (mW/mm^2)') 
%     ylabel('I_{opto} (nA)') 
% end

% X = 0:.01:3 ; 
% plot( X, Iopto(X), '-') ; 
% %plot( Popto(X), X)
% % Y = 0:5:Iopto(X(end)) ; 
% % plot( Popto(Y), Y, 'x','markersize',4) ; 

% drawnow ;
% hold off ;