R = [] ;
Y = [] ;
X=50:50:1000 ;
cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7] [0.7 0.7 0.2]} ;

for i=1:length(X)
    [lbd r] = RateStabilityK(4,X(i)) ;
    Y = [Y lbd] ;
    R = [R r] ;
end

dim = length(lbd) ;

subplot(3,1,1)
for i=1:length(r)
    plot(X,R(i,:),'-o','color',cl{i}) ; hold on ;
end
ylabel('Rates')
hold off ;

subplot(3,1,2)
for i=1:dim
    plot(X,real(Y(i,:)),'-o','color',cl{i}) ; hold on ;
end
ylabel('Re \lambda')
hold off ;

subplot(3,1,3)
for i=1:dim
    plot(X,imag(Y(i,:)),'-o','color',cl{i}) ; hold on ;
end
ylabel('Im \lambda')
hold off ;

xlabel('K')
