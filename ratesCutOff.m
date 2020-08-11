function [m idx RatesOut ROI] = ratesCutOff(Rates1, Rates2, THRESHOLD, Cth, DIM, L)
    
    nbN = length(Rates2) ; 
    idx = find(Rates2>=THRESHOLD) ; 
    RatesCO = Rates1(idx) ; 
    RatesOut = [] ;

    if( length(idx)==0 ) 
        m = .01 ; 
        return ; 
    else
        if(DIM==1)
                    
            X = linspace(-L/2, L/2, nbN) ;            
            % sizeROI = length( find(abs(X) <= Cth/2 ) );
            
            X = X(idx) ;
            ROI = find(abs(X)<=Cth) ;
            
            if( length(ROI)==0 )
                m = .01 ;
            else                
                RatesOut = RatesCO(ROI) ; 
                % fprintf(' %d %d ', idx(1), idx(length(idx)))                
                % m = sum(Rates) ./ sizeROI ; 
                m = mean(RatesOut) ; 
            end
            
        else

            X = [] ;
            Y = [] ;

            for i=1:nbN
                X(i) = L * mod( double(i), sqrt( double( nbN ) ) ) / sqrt( double( nbN ) ) ; 
                Y(i) = L * floor( double(i) / sqrt( double( nbN ) ) ) / sqrt( double( nbN ) ) ; 
            end
            
            % sizeROI = length( find( (X-L/2).^2 + (Y-L/2).^2 <= Cth .^2 / 4 ) ) ;
            
            X = X(idx) ;
            Y = Y(idx) ;
             
            ROI = find( (X-L/2).^2 + (Y-L/2).^2 <= Cth.^2 / 4 ) ; 
            
            if( length(ROI)==0 ) 
                m = .01 ;
            else                
                RatesOut = RatesCO(ROI) ;
                % RatesOut = RatesOut(RatesOut>=THRESHOLD) ;
                % fprintf(' %d %d ', idx(1), idx(length(idx))) 
                % m = sum(Rates) ./ sizeROI ; 
                m = max(mean(RatesOut),.01) ; 
            end
            
        end
    end

end