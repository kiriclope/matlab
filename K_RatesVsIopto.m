clear all ; 

IF_ROBUST = 1 ; 

l_K = [500 1000 2000] ; 
l_mk = ['^' '+' 'o'] ; 
l_alp = [.3 .6 .9] ; 

mkSize = 5 ; 

l_K = [50 100 500] ; 
l_mk = ['+' '^' 'o'] ; 
l_alp = [.3 .6 .9] ; 

for i=1:length(l_K) 
    K = l_K(i) ; 
    mk = l_mk(i) ; 
    alp = l_alp(i) ; 
    if(i==1) 
        IF_MF_RATES = 1 ; 
    else
        IF_MF_RATES = 0 ; 
    end
    RatesVsIopto 
end 
