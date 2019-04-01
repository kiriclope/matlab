clear all ;


% y = [.1 .2 .4 .9 1.8] ;
% x = [.8 1 1.2 1.25 1.5] *2 ;
% %myfit = fittype('a + b*log(x)','dependent',{'y'},'independent', {'x'},'coefficients',{'a','b'}) ;
% myfit = fittype('a*log(1+b*x)','dependent',{'y'},'independent', {'x'},'coefficients',{'a','b'}) ;
% % myfit = fittype('a*( 1 - exp(-x/b) )','dependent',{'y'},'independent', {'x'},'coefficients',{'a','b'}) ;


% y = [0 0.0955    0.1592    0.3183    0.4775    0.6366    1.0504 1.5915    2.5465    4.7746] ;
% x = linspace(0,1.5,length(y)) ;

% myfit = fittype('a/(1+exp(-(x-b)/c))','dependent',{'y'},'independent', {'x'},'coefficients',{'a','b','c'}) ;
% f = fit(x',y',myfit) 
% plot(f) ; 
% view([90 -90]) ;
% ylabel('\Gamma_{opto}')
% xlabel('I_{opto}')
% xlim([0 3])
% ylim([0 3])

x = [0 0.0955    0.1592    0.3183    0.4775    0.6366    1.0504 1.5915    2.5465    4.7746] ; 
y = linspace(0,1.5,length(x)) ; 
myfit = fittype('a*log(1+b*x)','dependent',{'y'},'independent', {'x'},'coefficients',{'a','b'}) ; 
f = fit(x',y',myfit) 
plot(f)
xlabel('\Gamma_{opto}') 
ylabel('I_{opto}') 
xlim([0 3])
ylim([0 3])
