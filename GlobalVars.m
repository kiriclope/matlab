warning off ; 
clear all ;

%%%%%%%%%%%%%%%%%%%%%

popList = ['E' 'I' 'S' 'V'] ;
cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ; 

%%%%%%%%%%%%%%%%%%%%%

model = 'LIF' ; 
nbpop = 2 ; 
dir = 'L5' ; 
file = 'MeanRates' ; 

Vl = -55. ; 
Vth = -40. ; 
Vr=[0.0 -80.] ; 
m0 = .001 ; 

%%%%%%%%%%%%%%%%%%%%% 

N = 2 ; 
g = 1 ; 
IF_Nk = 0 ; 

J = ImportJab(model,nbpop,dir) ; 
Iext = ExternalInput(model,nbpop,dir) ; 

for i=1:nbpop
    for j=1:nbpop
        G(i,j) = -(Vl-Vr(j))* abs(J(i,j)) ; 
    end
end

G0 = -Iext.*(Vl-Vr(1)) * m0 * 1000. ;
meanfield_rates = linsolve(G,-G0.') ;

fprintf('Mean field rates: ')
fprintf('%.3f ', meanfield_rates)
fprintf('\n')

if(nbpop==2)
    fprintf('Ie/Ii>Je/Ji')
    fprintf(' %.2f>', G0(1)/G(1,1) / (G0(2)/G(2,1) ) )
    fprintf('%.2f \n', G(1,2)/G(1,1) / (G(2,2)/G(2,1) ) )
end

nbN = nbNeuron(nbpop,N,IF_Nk,[]) ; 
Cpt = CptNeuron(nbpop,nbN) ; 

%%%%%%%%%%%%%%%%%%%%%
IF_CONDUCTANCES = 1 ;
RHO = 0.0 ;
%%%%%%%%%%%%%%%%%%%%%

nbIdv = [10 10 10 10] ; 

%%%%%%%%%%%%%%%%%%%%%

IF_LOOP = 0 ;
if(~IF_LOOP) 
    K = 400 ; 
    mkSize = 5 ;     
    IF_ROBUST = 0 ; 
    mk = 'o' ; 
    alp = 1 ; 
end

mk = 'o' ; 

%%%%%%%%%%%%%%%%%%%%%

THRESHOLD = .1 ; 
CV_THRESHOLD = .1 ; 

%%%%%%%%%%%%%%%%%%%%%

IF_NORM = 1 ; 
FIGPERPOP = 0 ; 
IF_IDVTraces = 1 ; 
IF_COUNTPIF = 0 ; 

IF_POWER = 2 ; 
I0 = 8. ; 
P0 = .5 ; 

IDX = 2 ; 
IF_CORRECTION = 0 ; 
IF_LOGSCALE = 1 ; 
IF_LOGSCALEX = 0 ; 

IF_MF_RATES = 0 ; 
IF_SAVE = 1 ; 

%%%%%%%%%%%%%%%%%%%%%

IF_PROP = 0 ; 
IF_PROPWEAK = 0 ; 

%%%%%%%%%%%%%%%%%%%%%

IF_IEXT = '' ; 
prtrPop = 2 ; 
Iprtr = Iext ; 

prtrAmp = 0 ; %.5 ; % .28 et .45

if(prtrPop>0)
    Iprtr(prtrPop) = Iprtr(prtrPop) + prtrAmp ; 
else
    Iprtr = prtrAmp * ones(1,nbpop) ;
end

v_Iprtr = 0:.1:1. ; 
%v_Iprtr = [.1:.1:.9,1:1:10 ]; 

%%%%%%%%%%%%%%%%%%%%%

IF_RING = '' ; 
L = 3 ; 
DIM = 1 ; 

IF_SPACELOOP = exist('IF_SPACELOOP') ;
if(~IF_SPACELOOP)
    Cff = [.1] ; 
end

Cff=.075;
Crec = [.125 .075 0.125 .075] ; 

IF_Dij = 0 ; 
Dij = [1.0 1.0 1.6667 1.0] ; 

Cth = 100 ; 

% Cff = [] ; 
% Crec = [.25 .25 .25 .075] ; 
% Cth = 100 ; 

v_Cff = .075:.025:.375 ;  
%v_Cff = .2:.025:.375 ; 

%%%%%%%%%%%%%%%%%%%%%
