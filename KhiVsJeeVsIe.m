clear all ; 
GlobalVars 

v_Jab = 0:.2:2 ; 
v_Ia = 0.5:.2:1.5 ; 
v_Ia = 2:.2:3 ; 

% v_Jab = .5:.1:2.5 ; 
% v_Ia = 2:.2:3. ; 

for l=1:length(v_Ia)

    % dir = sprintf('%s_Ie%.1f',dir0,v_Ia(l)) ;

    for i=1:length(v_Jab) 
        
        % BL = ImportData(model, nbpop, dir, 'IdvRates', N, K, g, IF_RING, Crec, Cff, 'Jab_Loop', prtrPop, [Iext(prtrPop) v_Jab(i)] ) ;
        
        % file = sprintf('Jee%.4f/Mean',v_Jab(i)) ; 
        %BL = ImportData(model, nbpop, dir, file, N, K, g, IF_RING, Crec, Cff, '', prtrPop, [] ) ; 

        try 
            BL = ImportData(model, nbpop, dir, 'Mean', N, K, g, IF_RING, ...
                            Crec, Cff, 'JabLoop', prtrPop, [Iext(prtrPop) ...
                                v_Ia(l) v_Jab(i)] ) ;
            BL(1,:) = [] ; 
        catch
            BL = NaN(1,nbpop+1)
        end

        for j=1:nbpop 
            mBL(l,j,i) = 0 ;
            mBL(l,j,i) = mean(BL(:,j+1)) ;
        end                                 

        fprintf('Jab ') 
        fprintf('%.3f ', v_Jab(i)) 

        fprintf(' RatesBL ') 
        for j=1:nbpop  
            fprintf('%.3f ', mBL(l,j,i) ) 
        end 
        fprintf('\n') 

        try
            Prtr = ImportData(model, nbpop, dir, 'Mean', N, K, g, IF_RING, ...
                              Crec, Cff, 'JabLoop', prtrPop, [Iext(prtrPop)+prtrAmp ...
                                v_Ia(l) v_Jab(i)] ) ;
            Prtr(1,:) = [] ;
        catch
            Prtr = NaN(1,nbpop+1)
        end

        for j=1:nbpop 
            mPrtr(j,i) = mean(Prtr(:,j+1)) ;
        end    

        fprintf(' RatesPrtr ') 
        for j=1:nbpop         
            fprintf('%.3f ', mPrtr(j,i) ) 
        end 
        fprintf('\n')
        
        fprintf(' Khi ') 
        for j=1:nbpop 
            khi(l,j,i) = ( mPrtr(j,i) - mBL(l,j,i) ) ./ Iprtr(prtrPop) ; 
            fprintf('%.3f ', khi(l,j,i) ) 
        end 
        fprintf('\n') 

        NormKhi(l,i) = khi(l,2,i) ./ khi(l,1,i) ; 
        mE(l,i) = mBL(l,1,i) ; 
        mI(l,i) = mBL(l,2,i) ; 
        % mS(l,i) = mBL(l,3,i) ; 
        % mV(l,i) = mBL(l,4,i) ; 

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

imagesc([20*v_Jab(1) 20*v_Jab(end)],[20*v_Ia(1) 20*v_Ia(end)],NormKhi)
%imagesc([0 30],[.4 2],NormKhi)
xlabel('J_{EE} (\mu A.ms.cm^-2)')
ylabel('J_{E0} (\mu A.ms.cm^-2)') 
h = colorbar ;
ylabel(h, ' \chi_I / \chi_E ')
caxis([-.25 1.5])
xlim([20*v_Jab(1) 20*v_Jab(end)])
ylim([20*v_Ia(1) 20*v_Ia(end)])

figtitle = sprintf('%s_mEVsJabVsIe',dir) ; 
fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 

imagesc([20*v_Jab(1) 20*v_Jab(end)],[20*v_Ia(1) 20*v_Ia(end)],mE)
xlabel('J_{EE} (\mu A.ms.cm^-2)')
ylabel('J_{E0} (\mu A.ms.cm^-2)') 
h = colorbar ;
caxis([0 20])
ylabel(h, ' r_E (Hz)')
xlim([20*v_Jab(1) 20*v_Jab(end)])
ylim([20*v_Ia(1) 20*v_Ia(end)])

figtitle = sprintf('%s_mIVsJabVsIe',dir) ; 
fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 

imagesc([20*v_Jab(1) 20*v_Jab(end)],[20*v_Ia(1) 20*v_Ia(end)],mI)
xlabel('J_{EE} (\mu A.ms.cm^-2)')
ylabel('J_{E0} (\mu A.ms.cm^-2)') 
h = colorbar ;
caxis([0 20])
ylabel(h, 'r_I (Hz)')
xlim([20*v_Jab(1) 20*v_Jab(end)])
ylim([20*v_Ia(1) 20*v_Ia(end)])

% figtitle = sprintf('%s_mEVsJabVsIe',dir) ; 
% fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 

% imagesc([20*v_Jab(1) 20*v_Jab(end)],[20*v_Ia(1) 20*v_Ia(end)],mS)
% xlabel('J_{EE} (\mu A.ms.cm^-2)')
% ylabel('J_{E0} (\mu A.ms.cm^-2)') 
% h = colorbar ;
% caxis([0 20])
% ylabel(h, ' r_S (Hz)')
% xlim([20*v_Jab(1) 20*v_Jab(end)])
% ylim([20*v_Ia(1) 20*v_Ia(end)])

% figtitle = sprintf('%s_mIVsJabVsIe',dir) ; 
% fig = figure('Name',figtitle,'NumberTitle','off') ; hold on ; 

% imagesc([20*v_Jab(1) 20*v_Jab(end)],[20*v_Ia(1) 20*v_Ia(end)],mV)
% xlabel('J_{EE} (\mu A.ms.cm^-2)')
% ylabel('J_{E0} (\mu A.ms.cm^-2)') 
% h = colorbar ;
% caxis([0 20])
% ylabel(h, 'r_V (Hz)')
% xlim([20*v_Jab(1) 20*v_Jab(end)])
% ylim([20*v_Ia(1) 20*v_Ia(end)])