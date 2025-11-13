%% use pre processed energy balane and crujra to run zeppetello model over most land pixels.
% check whether it all fits in RAM or needs to be done in chunks

clear


% test = ncinfo("crujra.v2.4.5d.spfh.1990.365d.noc.nc") ; 


cd('F:\CRUJRA\processed')

load('crujra_datetime_vector.mat')
load('datetime_1990_2014.mat')
load('row_cru_v.mat')
load('col_cru_v.mat')


SRB_file = matfile("Clara_SRB_full_array_linear.mat") ; 
tmp_file = matfile("tmp_full_array.mat") ;
pre_file = matfile("pre_full_array.mat") ;
spfh_file = matfile("spfh_full_array.mat") ;

% get VIP monthly NDVI for f_E
NDVI_file = matfile('F:\VIP_NDVI\VIP_NDVI_array_linear.mat') ; 


%tmp from ERA 5 monthly
cd('F:\ERA_5_monthly')
tmp_ERA_file = matfile("T2M_ERA_daily_interp.mat") ;
                    

load('F:\CRUJRA\processed\datetime_1990_2014.mat')

% srb_test = SRB_file.Clara_SRB_full_array(15000:16000,:)    ;
% tmp_test = tmp_file.tmp_full_array(15000:16000,:)          ;
% pre_test = pre_file.pre_full_array(15000:16000,:)          ;
% spfh_test = spfh_file.spfh_full_array(15000:16000,:)       ;
% NDVI_test = NDVI_file.VIP_NDVI_array (15000:16000,:)       ;


srb_test = SRB_file.Clara_SRB_full_array    ;
tmp_test = tmp_file.tmp_full_array          ;
pre_test = pre_file.pre_full_array          ;
spfh_test = spfh_file.spfh_full_array       ;
NDVI_test = NDVI_file.VIP_NDVI_array        ;

% 
% 
% Ts              = NaN(size(tmp_test,1),size(tmp_test,2)) ; 
% TC              = NaN(size(tmp_test,1),size(tmp_test,2)) ; 
% ms              = NaN(size(tmp_test,1),size(tmp_test,2)) ; 
% mD              = NaN(size(tmp_test,1),size(tmp_test,2)) ; 
% Tra_s_vec              = NaN(size(tmp_test,1),size(tmp_test,2)) ; 
% Tra_r_vec              = NaN(size(tmp_test,1),size(tmp_test,2)) ; 
% Es_vec              = NaN(size(tmp_test,1),size(tmp_test,2)) ; 
% Cap_flux_vec              = NaN(size(tmp_test,1),size(tmp_test,2)) ; 
% f_E_vec              = NaN(size(tmp_test,1),size(tmp_test,2)) ; 
% 
% 
% cd('F:\Matlab_Code\Zeppe_model')
% 
% for i = 1:1000
% 
% 
% steps_per_day = 2 ; 
% [Ts(i,:),TC(i,:),ms(i,:),mD(i,:),Tra_s_vec(i,:),Tra_r_vec(i,:),Es_vec(i,:),Cap_flux_vec(i,:), f_E_vec(i,:)] = ...
%         The_Model_f_E_vec(srb_test(i,:),tmp_test(i,:),spfh_test(i,:),pre_test(i,:)./ 86400 ./ 2,NDVI_test(i,:),steps_per_day) ; 
% 
% 
% i
% end
% 
% 
% 
% for i = 1:size(tmp_test,1)
% 
% 
% steps_per_day = 2 ; 
% [Ts(i,:),TC(i,:),ms(i,:),mD(i,:),Tra_s_vec(i,:),Tra_r_vec(i,:),Es_vec(i,:),Cap_flux_vec(i,:), f_E_vec(i,:)] = ...
%         The_Model(srb_test(i,:),tmp_test(i,:),spfh_test(i,:),pre_test(i,:)./86400./2,steps_per_day) ; 
% 
% 
% i
% end


% 
% 
% cd('F:\Matlab_Code\Zeppe_model')
% 
% 
% ms            = NaN(67420,9131) ; 
% f_E           = NaN(67420,9131) ; 
% 
% steps_per_day = 2 ; 
% 
% parfor i = 1:67420
% 
% srb_test = SRB_file.Clara_SRB_full_array(i,:)    ;
% tmp_test = tmp_file.tmp_full_array(i,:)          ;
% pre_test = pre_file.pre_full_array(i,:)         ;
% spfh_test = spfh_file.spfh_full_array(i,:)      ;
% NDVI_test = NDVI_file.VIP_NDVI_array(i,:)        ;
% 
% 
% [~,~,ms(i,:),~,~,~,~,~, f_E(i,:)] = ...
%         The_Model(srb_test,tmp_test,spfh_test,pre_test./86400./5,steps_per_day) ; 
% 
% i
% end
% 
% 
% 
% 
% test = sum(pre_file.pre_full_array,2,'omitnan') ; 
% test = test ./ 25 ; 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% for n = 1:200
%          % 
%          % Fig_Panel = figure ;
%          % plot(Ts(n,:) - 273,'g.-') ; hold on ; plot(TC(n,:) - 273,'k.-')  ;
%          % legend('Tsoil','Tcanopy')  
%          % xlabel('time')
%          % ylabel('Temperature [°Celsius]')
%          % set(gca,'FontSize',15)
%          % % saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model\Ts_Tc_1y_time_series','svg')
%          % % close 
%          % % yyaxis right
% 
%           Fig_Panel = figure ; 
%          plot(ms(n,:),'b.-') ; hold on ; %plot(mD(n,:),'r.-')  ;
%          % legend('SM surf','SM deep')  
%          xlabel('time')
%          ylabel('Soil moisture [m³/m³]')
%          set(gca,'FontSize',15)      
%          title(n)
%          pause(1)
%          close
% 
% end
% 
% 
% 


%% put mean back on map
% cd('F:\CRUJRA\processed')
% load('row_cru_v.mat')
% load('col_cru_v.mat')
% 
% 
% mean_ms = NaN(360,720) ; 
% 
% for i = 1:67420
% 
%     mean_ms(row_cru_v(i),col_cru_v(i)) = mean(ms(i,:),2,'omitnan')  ; 
% 
%     i
% end
% 
% 
% imagesc(mean_ms)
% [xs , ys] = getpts() ; xs = round(xs) ; ys= round(ys) ; 
% close
% 
% 
% for i = 1:length(xs)
% 
%     index = find( row_cru_v == ys(i) & col_cru_v == xs(i) ) ; 
% 
%          Fig_Panel = figure ; 
%          plot(ms(index,:),'b.-') ; hold on ; %plot(mD(n,:),'r.-')  ;
%          % legend('SM surf','SM deep')  
%          xlabel('time')
%          ylabel('Soil moisture [m³/m³]')
%          ylim([0 1])
%          set(gca,'FontSize',15)      
%          title(i)
%          % pause(2)
%          % close
% 
% end
% 
% 
% 


%% run model with f_e vector

%  spfh_testinfo = ncinfo('F:\CRUJRA\spfh\crujra.v2.4.5d.spfh.1990.365d.noc.nc') ; 


cd('F:\Matlab_Code\Zeppe_model')

% shorter with 1000 steps spinup
 Ts            = NaN(67420,9131+1000) ; 
 % TC            = NaN(67420,9131) ;
 ms             = NaN(67420,9131+1000) ;
 % mD            = NaN(67420,9131+1000) ;
 % f_E           = NaN(67420,9131) ;





% load('F:\projects\SM_long_term_DDs\Zeppe_model_output_03\means_ms_md.mat')
%  imagesc(mean_ms)  ; % clim([0 400])
% [xs , ys] = getpts() ; xs = round(xs) ; ys= round(ys) ; 
% close



% spatial extract on subset
% spatial_extract = row_cru_v == ys(1) &  col_cru_v == xs(1) ; 
% find(spatial_extract)

steps_per_day = 1 ;

cd('F:\Matlab_Code\Zeppe_model')


%i = 31785 % rainforest
% i = 61472 ; % = high coastal pixel
for i = 1:67420


%  i = spatial_extract    
srb_test = SRB_file.Clara_SRB_full_array_linear(i,:)        ;
tmp_test = tmp_file.tmp_full_array(i,:)              ;
pre_test = pre_file.pre_full_array(i,:)              ;
spfh_test = spfh_file.spfh_full_array(i,:)           ; % kg/kg
NDVI_test = NDVI_file.VIP_NDVI_array(i,:)            ;
  tmp_ERA_test = tmp_ERA_file.T2M_ERA_daily_interp(i,:) ;           
  tmp_ERA_test = tmp_ERA_test + 273.15 ; 

 if mean(NDVI_test,'omitnan') < 0  || range(NDVI_test) < 0.1 || all(isnan(NDVI_test))
     continue
 end

% figure
% plot(NDVI_test(100,:))


 % NDVI_test = 0.3 .* sin(linspace(0,10*pi,9131)) + 0.45 ; 
% movmean to make all cilmatic inputs basically monthly
 % spfh_test = movmean(spfh_test,32,2,'omitnan') ; 

% plot(tmp_test(2000:4000)-273.15) ; hold on  ; yyaxis right
% plot(spfh_test(2000:4000))
% 
% figure  ; plot(tmp_ERA_test(2000:4000)) ; hold on ; plot(tmp_test(2000:4000))  ; yyaxis right
% plot(srb_test(2000:4000))
%  figure ; plot(NDVI_test)

% srb_mean = repmat(mean(srb_test),[1 9131]) ; 


% daily CRUJRA temp
% [~,~,ms(i,:),mD(i,:),~,~,~,~, f_E(i,:)] = ...
        % The_Model_f_E_vec(srb_test,tmp_test ,spfh_test,pre_test./86400./2,NDVI_test,steps_per_day) ; 
 

% get 1000 length as spinup

srb_test =      cat(2,srb_test(1:1000),srb_test) ; 
tmp_test =     cat(2,tmp_test(1:1000),tmp_test) ; 
tmp_ERA_test =  cat(2,tmp_ERA_test(1:1000),tmp_ERA_test) ; 
spfh_test =     cat(2,spfh_test(1:1000),spfh_test) ;  
pre_test =      cat(2,pre_test(1:1000),pre_test) ; 
NDVI_test =     cat(2,NDVI_test(1:1000),NDVI_test) ; 


% monthly interpolated to daily ERA
 % [Ts(i,:),~,ms(i,:),~,~,~,~,~, ~] = ...
 %         The_Model_f_E_vec(srb_test,tmp_ERA_test ,spfh_test,pre_test./86400,NDVI_test,steps_per_day) ; 

 [Ts(i,:),~,ms(i,:),~,~,~,~,~, ~] = ...
         The_Model_f_E_vec(srb_test,tmp_ERA_test ,spfh_test,pre_test./86400,NDVI_test,steps_per_day) ;
 
i
end


test2 = mean(ms,1,'omitnan') ; 


Ts = Ts(:,1001:end) ;
ms = ms(:,1001:end) ;
mD = mD(:,1001:end) ;




% rescale to same volumetric units this is critical and should rescale all
% SM values to roughly the same scale as SMAP. We can also it to 0-0.6?
% This would make it basically the same as SMAP might be useful
theta_max = 0.6 ; % pore space
ms = ms .* theta_max ;
mD = mD .* theta_max ;


figure 
histogram(ms)
histogram(f_E)
histogram(Ts - 273.15)


test = min(Ts,[],2,'omitnan') ; 

test2 = sum(Ts < 0 ,2,'omitnan') ; 




% plot(srb_test, tmp_test,'.')
 f_E = (NDVI_test - (-0.1)) ./ (1 - (-0.1)) ; 


% plot(NDVI_test,'.')
% plot((NDVI_test - (-0.1)) ./ (1 - (-0.1)),'.')



% maybe short extra seciton here. Set Pixels that are not in ESA CCI SM NAN
% here. This might prevent wrong sampling cause of rainforest and desert
% areas that are not in ESA CCI
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_02')
% load('SM_CCI_Combined_mean_CRU.mat')

load('F:\CRUJRA\Ancil_data\CRUJRA_land_mask.mat')

for i = 1:67420
    dummy_ESA_CCI = CRUJRA_land_mask(row_cru_v(i),col_cru_v(i)) ; 
    if isnan(dummy_ESA_CCI)
     % Ts(i,:) = NaN ; 
     % TC(i,:) = NaN ; 
    ms(i,:) = NaN ; 
    % mD(i,:) = NaN ; 
    end
    i
end


test = sum(all(isnan(ms),2)) ; 




save('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\Ts','Ts','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\ms','ms','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\mD','mD','-v7.3')


% 
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output_04\mD','mD','-v7.3')
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output_04\TC','TC','-v7.3')
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output_04\f_E','f_E','-v7.3')
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\ms_ERA','ms','-v7.3')
% 


% save('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_nu_20\ms','ms','-v7.3')
% 
% 
% mean_ms = NaN(360,720) ; 
% mean_mD = NaN(360,720) ;
% mean_Ts = NaN(360,720) ;
% mean_TC = NaN(360,720) ;
% mean_f_E = NaN(360,720) ;
% 
% mean_ms = NaN(360,720) ;
% mode_ms = NaN(360,720) ;
% std_fe_anomalies =  NaN(360,720) ;
% 
% NDVI_mean =  NaN(360,720) ;
% 
% for i = 1:67420
% 
%     mean_ms(row_cru_v(i),col_cru_v(i)) = mean(ms(i,:),2,'omitnan')  ; 
%     mode_ms(row_cru_v(i),col_cru_v(i)) = mode(ms(i,:),2)  ;    
%      % mean_mD(row_cru_v(i),col_cru_v(i)) = mean(mD(i,:),2,'omitnan')  ;    
%     mean_Ts(row_cru_v(i),col_cru_v(i)) = mean(Ts(i,:),2,'omitnan')  ; 
%      % mean_TC(row_cru_v(i),col_cru_v(i)) = mean(TC(i,:),2,'omitnan')  ; 
% 
%      % mean_f_E(row_cru_v(i),col_cru_v(i)) = mean(f_E(i,:),2,'omitnan')  ; 
%       % std_fe_anomalies(row_cru_v(i),col_cru_v(i)) = std(f_E(i,:),1,'omitnan')  ;
%       % NDVI_test = NDVI_file.VIP_NDVI_array(i,:)            ;
%       % NDVI_mean(row_cru_v(i),col_cru_v(i)) = mean(NDVI_test,'omitnan')  ;
% 
%       % std_fe_anomalies(row_cru_v(i),col_cru_v(i)) = range(f_E(i,:),'omitnan')  ;
%     i
% end
% 
% 
% test = NDVI_mean ; 
% test(test < 0) = NaN ; 
% imagesc(test) ; figure ; imagesc(NDVI_mean)
% 
% 
% 
% 
% % h1 = pcolor(flipud(std_fe_anomalies)) ; set(h1,'LineStyle','none')
% test = std_fe_anomalies ; 
% test(test < 0.15) = NaN ; 
% figure ; h1 = pcolor(flipud(test)) ; set(h1,'LineStyle','none')
% 
% % h1 = pcolor(flipud(mean_ms)) ; set(h1,'LineStyle','none')
% % h1 = pcolor(flipud(mean_f_E)) ; set(h1,'LineStyle','none')
% % h1 = pcolor(flipud(std_fe_anomalies)) ; set(h1,'LineStyle','none')
% % h1 = pcolor(flipud(mean_Ts-273.15)) ; set(h1,'LineStyle','none')  ; clim([-0 50])
% 
% imagesc(mean_Ts-273.15) ;   clim([-10 40])
% [xs , ys] = getpts() ; xs = round(xs) ; ys= round(ys) ; 
% close
% 
% for i = 1:length(ys)
% 
%     index = find(row_cru_v == ys(i) & col_cru_v == xs(i)) ; 
% 
%     figure
%     plot(ms(index,3000:4000),'o-') ; ylim([0 1])
%     hold on 
%     yyaxis right 
%     plot(Ts(index,3000:4000)-273.15,'o-') ; ylim([-50 50])
% 
% end
% 
% 
% 
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output\mean_TC','mean_TC','-v7.3')
% 
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_03')
% 
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\ms','ms','-v7.3')
% 
% 
% figure
% imagesc(mean_Ts -273) ; clim([-40 40]) ;
% imagesc(mean_mD)
% 
% mean_ms(mean_ms > 0.9) = NaN ; 
% 
% figure ; imagesc(mean_ms)
% % imagesc(mode_ms)
% % imagesc(std_fe_anomalies)
% 
% [xs , ys] = getpts() ; xs = round(xs) ; ys= round(ys) ; 
% close
% 
% 
% for i = 1:length(xs)
% 
%     index = find( row_cru_v == ys(i) & col_cru_v == xs(i) ) ; 
% 
%          Fig_Panel(i) = figure('units','centimeters','position',[10 3 40 17])  ;
%          plot(datetime_1990_2014(5941:6941+365*6),ms(index,5941:6941+365*6),'b.-') ; 
%          hold on ;
%           % plot(datetime_1990_2014(5941:6941+365*6),mD(index,5941:6941+365*6),'r.-')  ;
%          ylim([0 1])
%          ylabel('Soil moisture [m³/m³]')
%          yyaxis right
%           % plot(datetime_1990_2014(5941:6941+365*6),T2M_CRU_daily_anomaly(index,5941:6941+365*6),'k.-')  ;
%           legend('SM surf','SM deep')  
%          xlabel('time')
%          % ylim([0 1])
%          set(gca,'FontSize',15)     
%          legend('SM_s','SM_d')
%          title(strcat('Sahel_',num2str(i)))
%          % pause(2)
%          % close
%          % yyaxis right
%          % plot(f_E(index,:),'g.-') 
% end
% 
% 
% %    plot(ms(index,:),T2M_CRU_daily_anomaly(index,:),'.')
% 
% 
% DDLength = 4 ; 
% cd('E:\Daten Baur\Matlab code\Project IGARSS multi frq tau\noodles_L_C_X')
% [NDry,timev,SMv] = DryDowns_SM(1:length(ms(index,5941:6941+365*6)),ms(index,5941:6941+365*6),DDLength)  ; 
% 
% hold on
% yyaxis left
% for i = 1:length(SMv)
%     SM_dummy = SMv{i} ; 
%     tt_dummy = timev{i} ; 
%      plot(datetime_1990_2014(5940 + tt_dummy),SM_dummy,'go-') ; 
%      hold on
% 
% end
% 
% xlim([datetime_1990_2014(5941), datetime_1990_2014(6941+365*6)])
% legend off
% 
% 
%  index = 102
%  plot(datetime_1990_2014(1:6941+365*6),T2M_CRU_daily_anomaly(index,1:6941+365*6),'k.-')  ;
%  hold on
%  plot(datetime_1990_2014(1:6941+365*6),movmean(T2M_CRU_daily_anomaly(index,1:6941+365*6),365),'g.-')  ;
% 
% 
% 
% for i = 1:size(Fig_Panel,2)
% 
% saveas(Fig_Panel(i),strcat('F:\projects\SM_long_term_DDs\figures\zeppe_model\SM_shallow_deep_time_series_Sahel_dewpoint',num2str(i)),'svg')
% 
% end
% close all
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% imagesc(mean_Ts)  ; clim([0 400])
% [xs , ys] = getpts() ; xs = round(xs) ; ys= round(ys) ; 
% close
% 
% 
% for i = 1:length(xs)
% 
%     index = find( row_cru_v == ys(i) & col_cru_v == xs(i) ) ; 
% 
%          Fig_Panel = figure ; 
%          plot(Ts(index,:)-273 ,'g.-') ; hold on ;
%          % plot(TC(index,:) - 273,'k.-')  ;
%          % legend('SM surf','SM deep')  
%          xlabel('time')
%          ylabel('Soil moisture [m³/m³]')
%          %ylim([0 1])
%          set(gca,'FontSize',15)      
%          title(i)
%          legend('Ts','Tc')
%          % pause(2)
%          % close
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% compute NDVI f_E anomalies. Maybe keep them as NDVI so it is easier to understand ..reverse minmax scalin
% 
% [y_NDVI, m_NDVI, d_NDVI] = ymd(datetime_1990_2014) ;
% 
% f_E_NDVI_daily_anomaly = NaN(67420,9131) ; 
% 
% for i = 1:size(f_E,1)
% 
% 
%     f_E_dummy = f_E(i,:) ; 
% 
% 
% 
% 
%         % remove very long term dynamics 5 year mean
%         VIP_NDVI_5y_movmean = movmean(f_E_dummy,1825) ; 
%         % detrend
%         % VIP_NDVI_day_interp = VIP_NDVI_day_interp - VIP_NDVI_5y_movmean ; 
% 
% 
%         % % remove from mean year
%         NDVI_extract_day_mean = NaN(365,1) ; 
%         for j = 1:365
%             [y, m, d] = ymd(datetime_1990_2014(j)) ; 
%             extract_index = find(d == d_NDVI & m == m_NDVI) ; 
%             NDVI_extract_day_mean(j) = mean(f_E_dummy(extract_index),'omitnan') ; 
% 
%         end
%         NDVI_extract_day_mean = repmat(NDVI_extract_day_mean,[34,1]) ; 
%         NDVI_extract_day_mean = NDVI_extract_day_mean(1:length(f_E_dummy)) ; 
% 
% 
%         f_E_NDVI_daily_anomaly(i,:) = f_E_dummy - NDVI_extract_day_mean' ; 
% 
% 
%        % plot(f_E_dummy,'-')    
%        % hold on
%        % plot(VIP_NDVI_5y_movmean,'-')
%        % plot(NDVI_extract_day_mean,'-')       
% 
%        % plot( f_E_NDVI_daily_anomaly(i,:))
% plot(f_E_dummy' - NDVI_extract_day_mean,'-')    
% 
% 
% 
% i
% end
% 
% 
% 
% 
% 
% 
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output\f_E_NDVI_daily_anomaly','f_E_NDVI_daily_anomaly','-v7.3')
% 
% 
% 
% %% get mean figures of forcing
% 
% 
% mean_tmp = NaN(360,720,1) ; 
% mean_srb = NaN(360,720,1) ; 
% mean_pre = NaN(360,720,1) ; 
% mean_spfh = NaN(360,720,1) ; 
% mean_dewtmp = NaN(360,720,1) ; 
% 
% 
% 
% for i = 1:67420
% 
%     % mean_tmp(row_cru_v(i),col_cru_v(i)) = mean(tmp_file.tmp_full_array(i,:),2,'omitnan')  ; 
%     % mean_srb(row_cru_v(i),col_cru_v(i)) = mean(SRB_file.Clara_SRB_full_array(i,:),2,'omitnan')  ; 
%     % mean_pre(row_cru_v(i),col_cru_v(i)) = mean(pre_file.pre_full_array(i,:),2,'omitnan')  ;     
%     % mean_spfh(row_cru_v(i),col_cru_v(i)) = mean(spfh_file.spfh_full_array(i,:),2,'omitnan')  ;     
% 
%     % mean_dewtmp(row_cru_v(i),col_cru_v(i)) = mean(Tdew(i,:),2,'omitnan')  ;     
% 
% 
%     i
% end
% 
% 
% 
% 
% 
% 
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output\mean_tmp','mean_tmp','-v7.3')
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output\mean_srb','mean_srb','-v7.3')
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output\mean_pre','mean_pre','-v7.3')
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output\mean_spfh','mean_spfh','-v7.3')
% 




%% compute T2M anomalies. 
cd('F:\CRUJRA\processed')

load('datetime_1990_2014.mat')
load('tmp_full_array.mat')

[y_NDVI, m_NDVI, d_NDVI] = ymd(datetime_1990_2014) ;

T2M_CRU_daily_anomaly = NaN(67420,9131) ; 


for i = 1:size(tmp_full_array,1)


    f_E_dummy = tmp_full_array(i,:) ; 




        % remove very long term dynamics 5 year mean
        VIP_NDVI_5y_movmean = movmean(f_E_dummy,1825) ; 
        % detrend
        % VIP_NDVI_day_interp = VIP_NDVI_day_interp - VIP_NDVI_5y_movmean ; 
   

        % % remove from mean year
        NDVI_extract_day_mean = NaN(365,1) ; 
        for j = 1:365
            [y, m, d] = ymd(datetime_1990_2014(j)) ; 
            extract_index = find(d == d_NDVI & m == m_NDVI) ; 
            NDVI_extract_day_mean(j) = mean(f_E_dummy(extract_index),'omitnan') ; 

        end
        NDVI_extract_day_mean = repmat(NDVI_extract_day_mean,[34,1]) ; 
        NDVI_extract_day_mean = NDVI_extract_day_mean(1:length(f_E_dummy)) ; 


        T2M_CRU_daily_anomaly(i,:) = f_E_dummy - NDVI_extract_day_mean' ; 


       % plot(f_E_dummy,'-')    
       % hold on
       % plot(VIP_NDVI_5y_movmean,'-')
       % plot(NDVI_extract_day_mean,'-')       

       % plot( T2M_CRU_daily_anomaly(i,:))
      % plot(f_E_dummy' - NDVI_extract_day_mean,'-')    



i
end



save('F:\projects\SM_long_term_DDs\Zeppe_model_output\T2M_CRU_daily_anomaly','T2M_CRU_daily_anomaly','-v7.3')




%% process monthly era .. so we use same anomalies and drivers for model and observaiton. Although daily Cruncep might be better
% just to check whether temp anomalies are different if we use era
clear


% monthly LAI will be only until 2014
cd('F:\ERA_5_monthly')
load('T2M_monthly_all.mat')
T2M_datetime = datetime('Jan-1979') ; 
load('T2m_datetime_ERA.mat') ; 
T2m_datetime_ERA = dateshift(T2m_datetime_ERA, 'start', 'day');
load('F:\CRUJRA\processed\datetime_1990_2014.mat') ; 

time_remove = T2m_datetime_ERA > max(datetime_1990_2014) | T2m_datetime_ERA < min(datetime_1990_2014) ; 
% time_remove(132) = 0 ; 
% time_remove(433) = 0 ; 


T2M_monthly_all(:,:,time_remove) = [] ; 
T2m_datetime_ERA(time_remove) = [] ; 


cd('F:\CRUJRA')
lon_cru = ncread('crujra.v2.4.5d.pre.2022.365d.noc.nc','lon')    ; 
lat_cru = ncread('crujra.v2.4.5d.pre.2022.365d.noc.nc','lat')    ; 
time_cru = ncread('crujra.v2.4.5d.pre.2022.365d.noc.nc','time')    ;

lat_cru = repmat(lat_cru,[1 720]) ; 
lon_cru = repmat(lon_cru',[360 1]) ; 


cd('F:\CLARA_surface_energy_budget\dataset-satellite-surface-radiation-budget-6bc9c8e5-baa4-4ec3-a219-683414681996')
fileslist_clara = string(ls('SRB*2014*')); 
lat_clara = ncread(fileslist_clara(1),'lat') ; 
lon_clara = ncread(fileslist_clara(1),'lon') ; 
lat_clara = repmat(lat_clara,[1 1440]) ; 
lon_clara = repmat(lon_clara',[720 1]) ; 


% pre allocate
T2M_monthly_all_interp = NaN(360,720,300) ; 

for i = 1:300

    dummy = T2M_monthly_all(:,:,i) ; 
    dummy = interp2(lon_clara,lat_clara,dummy,lon_cru,lat_cru) ; 
    T2M_monthly_all_interp(:,:,i) = dummy ; 
i
end




clear

% now put into 2D and inteprolate to daily like for the first analysis
cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_02')
load('T2M_monthly_all_interp.mat')
load('T2m_datetime_ERA.mat') ; 

cd('F:\CRUJRA\Ancil_data')
load('CRUJRA_land_mask.mat')
load('CRUJRA_land_mask_vec.mat')


T2m_datetime_ERA = dateshift(T2m_datetime_ERA, 'start', 'day');
load('F:\CRUJRA\processed\datetime_1990_2014.mat') ; 
time_remove = T2m_datetime_ERA > max(datetime_1990_2014) | T2m_datetime_ERA < min(datetime_1990_2014) ; 

T2m_datetime_ERA(time_remove) = [] ; 



    tmp_dummy = reshape(T2M_monthly_all_interp,[720*360,300]) ; 
    % remove ocean rows
    tmp_dummy(~CRUJRA_land_mask_vec,:) = [] ;
   [Lia, day_positions] = ismember(T2m_datetime_ERA, datetime_1990_2014') ; 



T2M_ERA_daily_interp = NaN(67420, 9131) ;

for i = 1:size(tmp_dummy,1)

T2M_ERA_daily_interp(i,:) = interp1(day_positions,tmp_dummy(i,:),1:9131,'linear')  ; 
T2M_ERA_daily_interp(i,9102:end) = T2M_ERA_daily_interp(i,9102-365:9131-365) ; 

i
end


save('F:\ERA_5_monthly\T2M_ERA_daily_interp','T2M_ERA_daily_interp','-v7.3') ; 




%% compute anomalies from T2M Era daily interpolated 



[y_NDVI, m_NDVI, d_NDVI] = ymd(datetime_1990_2014) ;
T2M_ERA_daily_anomaly = NaN(67420,9131) ; 


for i = 1:size(T2M_ERA_daily_interp,1)


    f_E_dummy = T2M_ERA_daily_interp(i,:) ; 



        % remove very long term dynamics 5 year mean
        VIP_NDVI_5y_movmean = movmean(f_E_dummy,1825) ; 
        % detrend
        % VIP_NDVI_day_interp = VIP_NDVI_day_interp - VIP_NDVI_5y_movmean ; 
   

        % % remove from mean year
        NDVI_extract_day_mean = NaN(365,1) ; 
        for j = 1:365
            [y, m, d] = ymd(datetime_1990_2014(j)) ; 
            extract_index = find(d == d_NDVI & m == m_NDVI) ; 
            NDVI_extract_day_mean(j) = mean(f_E_dummy(extract_index),'omitnan') ; 

        end
        NDVI_extract_day_mean = repmat(NDVI_extract_day_mean,[34,1]) ; 
        NDVI_extract_day_mean = NDVI_extract_day_mean(1:length(f_E_dummy)) ; 


        T2M_ERA_daily_anomaly(i,:) = f_E_dummy - NDVI_extract_day_mean' ; 


       % plot(f_E_dummy,'-')    
       % hold on
       % plot(VIP_NDVI_5y_movmean,'-')
       % plot(NDVI_extract_day_mean,'-')       

       % plot( f_E_NDVI_daily_anomaly(i,:))
% plot(f_E_dummy' - NDVI_extract_day_mean,'-')    



i
end



save('F:\ERA_5_monthly\T2M_ERA_daily_anomaly','T2M_ERA_daily_anomaly','-v7.3') ; 





%% compute NDVI  anomalies. Maybe keep them as NDVI so it is easier to understand ..reverse minmax scalin
clear



load('F:\CRUJRA\processed\datetime_1990_2014.mat')
load('F:\projects\SM_long_term_DDs\Zeppe_model_output\f_E.mat')


cd('F:\VIP_NDVI')
load('VIP_NDVI_array_linear.mat')



[y_NDVI, m_NDVI, d_NDVI] = ymd(datetime_1990_2014) ;

NDVI_daily_anomaly_5y_365seas = NaN(67420,9131) ; 

parfor i = 1:size(VIP_NDVI_array,1)


    f_E_dummy = VIP_NDVI_array(i,:) ; 




        % remove very long term dynamics 5 year mean
        VIP_NDVI_5y_movmean = movmean(f_E_dummy,1825) ; 
        % detrend
        % VIP_NDVI_day_interp = VIP_NDVI_day_interp - VIP_NDVI_5y_movmean ; 
   

        % % remove from mean year
        NDVI_extract_day_mean = NaN(365,1) ; 
        for j = 1:365
            [y, m, d] = ymd(datetime_1990_2014(j)) ; 
            extract_index = find(d == d_NDVI & m == m_NDVI) ; 
            NDVI_extract_day_mean(j) = mean(f_E_dummy(extract_index),'omitnan') ; 

        end
        NDVI_extract_day_mean = repmat(NDVI_extract_day_mean,[34,1]) ; 
        NDVI_extract_day_mean = NDVI_extract_day_mean(1:length(f_E_dummy)) ; 


        NDVI_daily_anomaly_5y_365seas(i,:) = f_E_dummy - NDVI_extract_day_mean' ; 


       % plot(f_E_dummy,'-')    
       % hold on
       % plot(VIP_NDVI_5y_movmean,'-')
       % plot(NDVI_extract_day_mean,'-')       

       % plot( f_E_NDVI_daily_anomaly(i,:))
% plot(f_E_dummy' - NDVI_extract_day_mean,'-')    



i
end




cd('F:\VIP_NDVI')

load('VIP_NDVI_array_linear.mat')

NDVI_daily_anomaly_02 = NaN(67420,9131) ; 


parfor i = 1:size(VIP_NDVI_array,1)


    VIP_NDVI_array_dummy = VIP_NDVI_array(i,:) ; 
    VIP_NDVI_movmean = movmean(VIP_NDVI_array_dummy,121) ; 
    NDVI_daily_anomaly_02(i,:) = VIP_NDVI_array_dummy - VIP_NDVI_movmean ; 

    % plot(f_E_dummy) ; hold on 
    % plot(VIP_NDVI_movmean) ; plot(f_E_dummy - VIP_NDVI_movmean)


i
end



histogram(randsample(NDVI_daily_anomaly_02(:),100000))
hold on
histogram(randsample(NDVI_daily_anomaly(:),100000))






save('F:\projects\SM_long_term_DDs\Zeppe_model_output\f_E_NDVI_daily_anomaly','f_E_NDVI_daily_anomaly','-v7.3')
save('F:\VIP_NDVI\NDVI_daily_anomaly','NDVI_daily_anomaly','-v7.3') ; 

save('F:\VIP_NDVI\NDVI_daily_anomaly_5y_365seas','NDVI_daily_anomaly_5y_365seas','-v7.3') ; 


% here saves for tests runs with 
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\ms','ms','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\f_E_NDVI_daily_anomaly','f_E_NDVI_daily_anomaly','-v7.3')












