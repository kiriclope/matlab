clear all ;
GlobalVars

J = ImportJab(model,nbpop,dir,0) 

X=0:.01:1 ;

for i=1:length(X)
    for j=1:length(X)
        Det(i,j) = det(J) ;
        lbd(i,j) = sign( max( real(eig(J)) ) );
        if(det(J)>0)
            Slope(i,j) = sign( -X(i)*J(4,3) - X(j)*J(4,1) );
        else
            Slope(i,j) = 0 ;
        end
    end 
end

imagesc(X,X,Slope) ; hold on ;
plot(J(1,1),-J(1,3),'k*')
hold off ;
h = colorbar ;
caxis([-1 1])