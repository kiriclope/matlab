clear all ;

IF_PAPERMAIN = 1 ;

l_dir2pop={'L5C' 'L5B'} ;
l_dir4pop={'S1L5'} ;

l_nbpop=[2 4] ; 
l_Iprtr=[.2 .4] ;

for idxPop=1:length(l_nbpop)    
    nbpop = l_nbpop(idxPop) ;
    if(nbpop==2) 
        for idx_dir=1:length(l_dir2pop)
            dir = l_dir2pop{idx_dir} ;
            %RatesVsIopto 
            % RatesDist
            % CVDist
        end
    else
        for idx_dir=1:length(l_dir4pop)
            dir = l_dir4pop{idx_dir} ;
            %RatesVsIopto 
        end

        for idx_dir=1:length(l_dir4pop)
            dir = l_dir4pop{idx_dir} ;
            for idx_prtr=1:length(l_Iprtr) 
                Iprtr = l_Iprtr(idx_prtr) ; 
                PieChart 
            end
        end

        % for idx_dir=1:length(l_dir4pop)
        %     dir = l_dir4pop{idx_dir} ;
        %     % RatesDist
        %     CVDist
        % end
        
    end
end


IF_PAPERMAIN = 0 ;
