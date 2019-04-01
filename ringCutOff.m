function Rates = ringCutOff(Cth, DIM, L)
    
    if(DIM==1)

        X = linspace(-L/2,L/2,length( IdvRates(i, Cpt(j)+1:Cpt(j+1) ) ) ) ;
        Idv = IdvRates(i,Cpt(j)+1:Cpt(j+1) ) ;
        
        ROI = length(find(abs(X)<=Cth)) ;
        idx = find(Idv>=THRESHOLD) ;
        
        if( length(idx)==0 )
            m(j,i) = 0 ;
        else
            
            X = X(idx) ;
            Idv = Idv(idx) ;
            idx = find(abs(X)<=Cth) ; % Cff
            
            if( length(idx)==0 )
                m(j,i) = 0 ;
            else
                
                X = X(idx) ;
                Idv = Idv(idx) ;
                fprintf(' %d %d ',idx(1),idx(length(idx)))
                
                m(j,i) = sum( Idv )./nbROI ;
            end
        end

end