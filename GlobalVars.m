warning off ; 
model = 'Auta' ; 
nbpop = 1 ; 
dir = 'Test' ; 
popList = ['E' 'I' 'S' 'V'] ; 

N = 1 ; 
K = 500 ; 
g = 10 ; 

IF_Nk = 0 ; 

IF_IEXT = 'AUTA' ; 
IF_DATA = 'AUTA' ; 
IF_ROBUST = 0 ; 
mk = '^' ; 
alp = .75 ; 

THRESHOLD = 0 ; 
CV_THRESHOLD = .1 ; 

prtrPop = 1 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% model = 'Binary' ; 
% nbpop = 2 ; 
% dir = 'Test' ; 

% N = 1 ; 
% K = 500 ; 
% g = 1 ; 

% IF_Nk = 0 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

L = 1.5 ; 

if(N==7 || N==5)
    IF_Nk = 1 ; 
    IF_RING = 'Gauss2D' ; 
    DIM = 2 ; 
    Crec = [.375 .125 .375 .25] ; 
    if(nbpop==2) 
        Crec = [.375 .125] ; 
    end
    Cff = 2 ; 
    Cth = 2 ; 
else
    IF_RING = '' ; 
    DIM = 1 ; 
    Cff = 2 ; 
    Crec = [] ; 
    Cth = 100 ; 
end

m0 = .20 ; 
r_Neuron = .003 ; % cm
A_Neuron = 4*pi*r_Neuron*r_Neuron * 4 ; 

Iprtr = 0 ; 
v_Iprtr = [0 .05 1.5] ; 
ITR = 2 ; 
%v_Iprtr = [0 .001 1] ;  

v_Jab = [0 .05 1.5] ; 

IF_NORM = 1 ; 
IF_IDVTraces = 1 ; 
IF_COUNTPIF = 0 ; 
IF_LOGSCALE = 0 ; 
IF_LOGSCALEX = 0 ; 
FIGPERPOP = 1 ; 
IF_POWER = 0 ; 

if K==Inf
    v_Iprtr = [0 .0025 .2] ; 
    IF_IDVTraces = 0 ; 
    IF_COUNTPIF = 0 ; 
    IF_LOGSCALE = 0 ; 
    IF_LOGSCALEX = 0 ; 
    FIGPERPOP = 0 ; 
    IF_POWER = 0 ; 
end
    
nbIdv = 10 ; 
IDX = 1 ; 

if IF_POWER==1 
    if strfind(dir,'L23') 

        P0 = 1.25*.25 ; 
        I0 = 5*.25 ; 
        IDX = 1 ;

        Pinf = 4641 ;
        I0 = 9.221 ;
        Iinf= 0.8546 ;

        Pinf = 1490 ;
        I0 = 14.66 ;
        Iinf= 1.625 ; 

        Pinf = 3082 ;
        I0 = 3.792 ;
        Iinf= 0.3512 ;
       
    elseif strfind(dir,'S1L5')
        P0 = 3.75 ; 
        I0 = 2.5 ;
        IDX = 1 ;

        Pinf = 4641 ;
        I0 = 9.221 ;
        Iinf= 0.8546 ;
    elseif( strfind(dir,'L5') ) 
        P0 = 3.75 ; 
        I0 = 2.5 ;
        IDX = 1 ;
        
        Pinf = 4641 ;
        I0 = 9.221 ;
        Iinf= 0.8546 ;

        Pinf = 2432 ;
        I0 = 4.94 ;
        Iinf= 0.4699 ;
        
        Pinf = 3082 ;
        I0 = 3.792 ;
        Iinf= 0.3512 ;
        
    end
elseif IF_POWER==2
    if strfind(dir,'L23') 

        P0 = .1 ; 
        I0 = 1 ;
        IDX = 1 ;

        P0 = 1.25*.25 ; 
        I0 = 5*.25 ; 
        IDX = 1 ;
        
        I0 = 8.006 ;
        P0 = 0.4281 ;

    elseif strfind(dir,'S1L5')
        P0 = 3.75 ; 
        I0 = 2.5 ;
        IDX = 1 ;

        P0 = .5 ; 
        I0 = 10 ;
        IDX = 1 ;

        I0 = 24.0 ;
        P0 = .8 ;
    else
        IDX = 1 ;
        I0 = 2.64 ;
        P0 = 0.13 ;

        I0 = 8.006 ;
        P0 = 0.4281 ;
    end
end

cl = {[1 0 0] [0 0 1] [0 1 0]  [0.7 0.7 0.7]} ; 

IF_SAVE = 0 ; 
