clear all ;
GlobalVars

IextBL = ExternalInput(model,nbpop,dir) ; 
nbN = nbNeuron(nbpop,N,IF_Nk,[]) ;
Cpt = CptNeuron(nbpop,nbN) ;

v_Jab = v_Jab(1):v_Jab(2):v_Jab(3) ; 
nbIdv = 10 ;

J = ImportJab(model,nbpop,dir) ; 
Iext = IextBL ; 
Iext(prtrPop) = IextBL(prtrPop) + Iprtr ; 

dir0 = dir ;

v_Ia = .8:.2:2 ;

%v_Ia = 1 ;
for l=1:length(v_Ia)

    dir = sprintf('%s_Ie%.1f',dir0,v_Ia(l)) ;
    
    for i=1:length(v_Jab) 
        
        BL = ImportData(model, nbpop, dir, 'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, IextBL(prtrPop), 0, 0 , 1 , v_Jab(i) ) ;
        
        try
            for j=1:length(BL(1,:))-1
                IdvRatesBL(i,j) = mean( BL(:,j+1) ) ;
            end
        catch
            for j=1:nbN(nbpop)
                IdvRatesBL(i,j) = nan ; 
            end
        end

        for j=1:nbpop 
            mBL(l,j,i) = 0 ;
            RatesBL = IdvRatesBL(i, Cpt(j)+1:Cpt(j+1)) ;
            % mBL(j,i) = mean(RatesBL) ; 
            [mBL(l,j,i) idx] = ratesCutOff(RatesBL, RatesBL, THRESHOLD, Cth, DIM, L) ;
        end                                 

        fprintf('Jab ') 
        fprintf('%.3f ', v_Jab(i)) 

        fprintf(' RatesBL ') 
        for j=1:nbpop         
            fprintf('%.3f ', mBL(l,j,i) )
        end    
        fprintf('\n')

        Prtr = ImportData(model, nbpop, dir, 'IdvRates', N, K, g, IF_RING, Crec, Cff, IF_DATA, prtrPop, Iext(prtrPop), 0, 0 , 1 , v_Jab(i)) ;
        try
            for j=1:length(Prtr(1,:))-1
                IdvRatesPrtr(i,j) = mean(Prtr(:,j+1)) ;
            end
        catch
            for j=1:nbN(nbpop)
                IdvRatesPrtr(i,j) = nan ; 
            end
        end

        for j=1:nbpop 
            mPrtr(j,i) = 0 ;
            RatesPrtr = IdvRatesPrtr(i, Cpt(j)+1:Cpt(j+1)) ;
            % mPrtr(j,i) = mean(RatesPrtr) ;
            [mPrtr(j,i) idx] = ratesCutOff(RatesPrtr, RatesPrtr, THRESHOLD, Cth, DIM, L) ;
        end    

        fprintf(' RatesPrtr ') 
        for j=1:nbpop         
            fprintf('%.3f ', mPrtr(j,i) )
        end 
        fprintf('\n')
        

        fprintf(' Khi ') 
        for j=1:nbpop 
            khi(l,j,i) = ( mPrtr(j,i) - mBL(l,j,i) ) ./ Iprtr ; 
            fprintf('%.3f ', khi(l,j,i) ) 
        end 
        fprintf('\n') 

        NormKhi(l,i) = khi(l,2,i) ./ khi(l,1,i) ; 
        mE(l,i) = mBL(l,1,i) ;
        mI(l,i) = mBL(l,2,i) ;

    end
end

% figtitle = sprintf('L5_KhiVsJab') ; 
% fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 
% xlabel('J_{EE}')
% ylabel('Norm. \chi ')
% plot(v_Jab, zeros(1,length(v_Jab)), '--','Color','k')

% for i=1:nbpop
%     Y = khi(1,i,:) ./ mBL(1,i,:) ; 
%     plot(v_Jab, squeeze(Y), 'd','Color',cl{i},'markersize',1)     
% end

figtitle = sprintf('%s_KhiVsJabVsIe',dir) ; 
fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 

imagesc([0 3.0],[.4 1],NormKhi)
xlabel('J_E')
ylabel('I_E')
h = colorbar ;
ylabel(h, ' \chi_I / \chi_E ')
%caxis([0 1])

figtitle = sprintf('%s_mEVsJabVsIe',dir) ; 
fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 

imagesc([0 3.0],[.4 1],mE)
xlabel('J_E')
ylabel('I_E')
h = colorbar ;
ylabel(h, ' m_E ')
%caxis([0 1])

figtitle = sprintf('%s_mIVsJabVsIe',dir) ; 
fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 

imagesc([0 3.0],[.4 1],mI)
xlabel('J_E')
ylabel('I_E')
h = colorbar ;
ylabel(h, 'm_I ')
