x = [.03 .1 .2 .5 1] ;
y = [.011 .042 .093 .21 .28] ;
err = [0.273 0.143 .194 .095 .179] ;
logerr = 0.434 .* err ;

figure() ;
errorbar(x,y,err)

figure() ;
errorbar(x,y,logerr)
