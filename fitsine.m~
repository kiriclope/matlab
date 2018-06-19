function out = fitsine(x,y) 

    warning off ;
    yu = max(y);
    yl = min(y);
    yr = (yu-yl); % Range of ‘y’
    yz = y-yu+(yr/2);
    % zx = x(yz .* circshift(yz,[0 1]) <= 0) % Find zero-crossings
    % per = 2*mean(diff(zx)) % Estimate period
    ym = mean(y);  % Estimate offset
    
    fitfcn = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4); % Function to fit
    fcn = @(b) sum((fitfcn(b,x) - y).^2);   % Least-Squares cost function
    
    options = optimset('MaxFunEvals',10000,'MaxIter',10000);

    per =  pi+2*pi*rand ;
    
    exitflag = -1 ;
    TOLERANCE = 10E-5 ;
    out = zeros(1,4) ;

    while exitflag<=0 | fcn(out)<TOLERANCE
        [out,fval,exitflag,output] = fminsearch(fcn, [yr;  per;  -1;  ym],options) ; 
        per = 2*pi*rand ;
    end

end