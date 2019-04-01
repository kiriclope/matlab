Vr=0 ;
Vth = 1 ;
nbSteps = 100 ;
Vk = linspace(Vr,Vth,nbSteps) ;
Delta = abs(Vk(1) - Vk(2)) ;

for i=1:nbSteps
    pk(i) = pk(i+1) * exp(Delta * G(i+1) ) + Delta * H(i+1)