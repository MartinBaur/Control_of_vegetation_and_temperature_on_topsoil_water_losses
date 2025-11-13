%% MJB 11-5-2024 use correct albedo runs and redo analyses.
clear


%% extract drydowns for 02 run


cd('F:\CRUJRA\processed')

load('row_cru_v.mat')
load('col_cru_v.mat')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output')
% load('ms.mat')



 cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests')
 load('ms_ERA.mat')
 % load('f_E.mat')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08')
%  load('ms.mat')


% load('means_global.mat')


% model is every 12 h .. so maybe do length of 8? Not gonna make any
% difference but fpr consistency
DDLength = 4  ;  % keep consistent with SMAP. 4 day min DD
sminterp = linspace(0,0.6,60) ; 
NDry_count = 0 ; 
count = 1 ;

% load('md')
% load('DD_cru_cols.mat')
% load('DD_cru_rows.mat')


% load fe anomaly from 01 run. Does not change
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output\')
% NDVI anomaly
% load('f_E_NDVI_daily_anomaly.mat')


% T2M CRUJRA anomaly
 % load('T2M_CRU_daily_anomaly.mat')


% T2M CRUJRA 
 % load('F:\CRUJRA\processed\tmp_full_array.mat')

% NDVI anomaly
% load('F:\VIP_NDVI\NDVI_daily_anomaly.mat')
load('F:\VIP_NDVI\NDVI_daily_anomaly_5y_365seas.mat')


% NDVI 
% load('F:\VIP_NDVI\VIP_NDVI_array_linear.mat')




% alternative ERA5 monthly interpolated anomalies
 cd('F:\ERA_5_monthly')
  load('T2M_ERA_daily_anomaly.mat')


% %% load tests
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests')
% load('f_E_NDVI_daily_anomaly.mat')
% load('ms.mat')


%%

DD_dSM_dt_shallow_interpsm = NaN(12900000,60) ; 
% DD_dSM_dt_deep_interpsm = NaN(12900000,60) ; 


% DD_NDVI_interpsm = NaN(11000000,50) ; 
 % DD_dSM_dt_deep_interpsm = NaN(4000000,50) ; 
% DD_Ts_interpsm = NaN(4000000,50) ; 
 % DD_f_E_anomaly_interpsm = NaN(12000000,50) ; 
  DD_T2M_ERA_anomaly_interpsm = NaN(12900000,60) ; 
  DD_NDVI_anomaly_interpsm = NaN(12900000,60) ; 

 % DD_tmp_CRUJRA_interpsm = NaN(11000000,50) ; 
% DD_T2M_ERA_anomaly_interpsm = NaN(10000000,50) ; 


DD_cru_rows = NaN(12900000,1) ; 
DD_cru_cols = NaN(12900000,1) ; 
DD_NDry = NaN(12900000,1) ; 
DD_timev =  NaN(12900000,2) ; 
load('F:\CRUJRA\processed\datetime_1990_2014.mat')



cd('E:\Daten Baur\Matlab code\Project IGARSS multi frq tau\noodles_L_C_X')
count = 1 ; 
 

for i = 1:length(ms)

time = 1:9131 ; 
ms_dummy = ms(i,:) ; 

 row_dummy = row_cru_v(i,:) ; 
 col_dummy = col_cru_v(i,:) ; 
 % mD_dummy = mD(i,:) ;
 % Ts_dummy = Ts(i,:) ; 
 % f_E_anomaly_dummy = f_E_NDVI_daily_anomaly(i,:) ; 
 % tmp_dummy = tmp(i,:) ;
  NDVI_anomaly_dummy = NDVI_daily_anomaly_5y_365seas(i,:) ; 


  % NDVI_dummy = VIP_NDVI_array(i,:) ; 
  % tmp_dummy = tmp_full_array(i,:) ;

   % tmp_anomaly_dummy = T2M_CRU_daily_anomaly(i,:) ; 
   tmp_anomaly_dummy = T2M_ERA_daily_anomaly(i,:) ; 


% maybe filter out weird coastal pixels and ponded conditions. They all
% reach fully ponded conditions quickly check the implication of this
 % if (mode(ms_dummy) == 1 || all(isnan(ms_dummy)))
 if (all(isnan(ms_dummy)))
     continue
 end


 [NDry,timev,SMv] = DryDowns_SM(time,ms_dummy,DDLength)  ; 
 NDry_count = NDry_count + NDry ;

% process them right after detection
 for j = 1:NDry

     DD_SM = SMv{j,1} ; 
     % remove DDs with little change?
      if range(DD_SM) < 0.05 
          continue
      end

     % rather remove drydowns with full saturation?
        if any(DD_SM > 0.59)
            continue
        end

      % if any(DD_SM == 1)
          % continue
      % end

      % if any(DD_SM == 0)
          % continue
      % end    

     DD_tt = timev{j,1} ; 
     DD_dSM_dt = diff(DD_SM) ; 
     %%%%%%%%%%%%%%%%  !!!!!!!!!!!!!!  %%%%%%%%%%%%%%%%%%%%%%%%%%

     % dt in model is 12 h, so multiply dSM/dt by 2 for loss rate per day
     % DD_dSM_dt = DD_dSM_dt .* 2 ; 

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % only use strictly decreasing dds
     if any(DD_dSM_dt > 0)
         continue
     end



     DD_SM_short = DD_SM(1:end-1) ;  

     if numel(DD_SM_short)~=numel(unique(DD_SM_short))
         continue
     end


     DD_dSM_dt_shallow_interpsm(count,:) = interp1(DD_SM_short,DD_dSM_dt,sminterp,'linear',NaN) ;   

      % DD_NDVI = NDVI_dummy(DD_tt(1:end-1)) ; 
      % DD_NDVI_interpsm(count,:) = interp1(DD_SM_short,DD_NDVI,sminterp,'linear',NaN) ;  

      % DD_mD = mD_dummy(DD_tt) ; 
      % DD_dmD_dt = diff(DD_mD) ; 
      % DD_dSM_dt_deep_interpsm(count,:) = interp1(DD_SM_short,DD_dmD_dt,sminterp,'linear',NaN) ;  

      % DD_Ts = Ts_dummy(DD_tt(1:end-1)) ; 
      % DD_Ts_interpsm(count,:) = interp1(DD_SM_short,DD_Ts,sminterp,'linear',NaN) ;  


       % DD_f_E_anomaly = f_E_anomaly_dummy(DD_tt(1:end-1)) ; 
       % DD_f_E_anomaly_interpsm(count,:) = interp1(DD_SM_short,DD_f_E_anomaly,sminterp,'linear',NaN) ;  

       % DD_tmp = tmp_dummy(DD_tt(1:end-1)) ; 
       % DD_tmp_CRUJRA_interpsm(count,:) = interp1(DD_SM_short,DD_tmp,sminterp,'linear',NaN) ;  

       % tmp_anomaly = tmp_anomaly_dummy(DD_tt(1:end-1)) ; 
       % DD_T2M_CRU_anomaly_interpsm(count,:) = interp1(DD_SM_short,tmp_anomaly,sminterp,'linear',NaN) ;  

        DD_NDVI_anomaly = NDVI_anomaly_dummy(DD_tt(1:end-1)) ; 
        DD_NDVI_anomaly_interpsm(count,:) = interp1(DD_SM_short,DD_NDVI_anomaly,sminterp,'linear',NaN) ;  

        tmp_anomaly = tmp_anomaly_dummy(DD_tt(1:end-1)) ; 
        DD_T2M_ERA_anomaly_interpsm(count,:) = interp1(DD_SM_short,tmp_anomaly,sminterp,'linear',NaN) ;  

      DD_cru_rows(count) = row_dummy ; 
      DD_cru_cols(count) = col_dummy ;
      DD_timev(count,1) = DD_tt(1) ;
      DD_timev(count,2) = DD_tt(end) ;      


     % if too short in SM space skip and not save
     if all(isnan(DD_dSM_dt_shallow_interpsm(count,:)))
        continue
     end



     count = count +1 ; 



end
i

end






% save 02 run

DD_dSM_dt_shallow_interpsm(isnan(DD_cru_rows),:)= [] ; 
% DD_T2M_CRU_anomaly_interpsm(isnan(DD_cru_cols),:) = [] ; 
% DD_f_E_anomaly_interpsm(isnan(DD_cru_cols),:) = [] ; 
% DD_tmp_CRUJRA_interpsm(isnan(DD_cru_cols),:) = [] ; 
DD_T2M_ERA_anomaly_interpsm(isnan(DD_cru_cols),:) = [] ; 
DD_NDVI_interpsm(isnan(DD_cru_cols),:) = [] ; 
% DD_NDVI_anomaly_interpsm(isnan(DD_cru_cols),:) = [] ; 
DD_cru_rows(isnan(DD_cru_rows)) = [] ; 
DD_cru_cols(isnan(DD_cru_cols)) = [] ; 




DD_dSM_dt_shallow_interpsm(isnan(DD_cru_rows),:)= [] ; 
DD_NDVI_anomaly_interpsm(isnan(DD_cru_rows),:)= [] ; 
DD_T2M_ERA_anomaly_interpsm(isnan(DD_cru_cols),:) = [] ; 
% DD_T2M_CRU_anomaly_interpsm(isnan(DD_cru_cols),:) = [] ; 
DD_timev(isnan(DD_cru_cols),:) = [] ; 
DD_cru_rows(isnan(DD_cru_rows)) = [] ; 
DD_cru_cols(isnan(DD_cru_cols)) = [] ; 



% cut to same length if needed
DD_dSM_dt_shallow_interpsm = DD_dSM_dt_shallow_interpsm(1:12000000,:) ; 
DD_NDVI_anomaly_interpsm = DD_NDVI_anomaly_interpsm(1:12000000,:) ; 
DD_T2M_ERA_anomaly_interpsm = DD_T2M_ERA_anomaly_interpsm(1:12000000,:) ; 
% DD_T2M_CRU_anomaly_interpsm = DD_T2M_CRU_anomaly_interpsmDD_timev(1:12000000,:) ;
DD_timev = DD_timev(1:12000000,:) ;
DD_cru_rows = DD_cru_rows(1:12000000,:) ;
DD_cru_cols = DD_cru_cols(1:12000000,:) ;




% clear model output after DD search?
clear mD ms TC Ts


save('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\DD_dSM_dt_shallow_interpsm','DD_dSM_dt_shallow_interpsm','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\DD_NDVI_anomaly_interpsm','DD_NDVI_anomaly_interpsm','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\DD_T2M_ERA_anomaly_interpsm','DD_T2M_ERA_anomaly_interpsm','-v7.3')

 % save('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\DD_T2M_CRU_anomaly_interpsm','DD_T2M_CRU_anomaly_interpsm','-v7.3')

save('F:\projects\SM_long_term_DDs\Zeppe_model_output_03\DD_tmp_CRUJRA_interpsm','DD_tmp_CRUJRA_interpsm','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_03\DD_NDVI_interpsm','DD_NDVI_interpsm','-v7.3')


save('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\DD_cru_rows','DD_cru_rows','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\DD_cru_cols','DD_cru_cols','-v7.3')

save('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\DD_timev','DD_timev','-v7.3')



DD_dSM_dt_shallow_interpsm = DD_dSM_dt_shallow_interpsm ./ 2 ; 





save('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\DD_dSM_dt_shallow_interpsm','DD_dSM_dt_shallow_interpsm','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\DD_NDVI_anomaly_interpsm','DD_NDVI_anomaly_interpsm','-v7.3')
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\DD_T2M_CRU_anomaly_interpsm','DD_T2M_CRU_anomaly_interpsm','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\DD_T2M_ERA_anomaly_interpsm','DD_T2M_ERA_anomaly_interpsm','-v7.3')

save('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\DD_cru_rows','DD_cru_rows','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\DD_cru_cols','DD_cru_cols','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\DD_timev','DD_timev','-v7.3')




% for tests runs

save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\DD_dSM_dt_shallow_interpsm','DD_dSM_dt_shallow_interpsm','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\DD_NDVI_anomaly_interpsm','DD_NDVI_anomaly_interpsm','-v7.3')
% save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\DD_T2M_CRU_anomaly_interpsm','DD_T2M_CRU_anomaly_interpsm','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\DD_T2M_ERA_anomaly_interpsm','DD_T2M_ERA_anomaly_interpsm','-v7.3')

save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\DD_cru_rows','DD_cru_rows','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\DD_cru_cols','DD_cru_cols','-v7.3')
save('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\DD_timev','DD_timev','-v7.3')


% %%
% clear
% 
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_03')
% 
% % load('DD_cru_cols.mat')
% % load('DD_cru_rows.mat')
% % load('DD_dSM_dt_shallow_interpsm.mat')
% % load('DD_f_E_anomaly_interpsm.mat')
% % load('DD_dSM_dt_deep_interpsm.mat')
% % load('DD_T2M_CRU_anomaly_interpsm.mat')
% % load('means_global.mat')
% 
% 
% 
% load('DD_cru_cols.mat')
% load('DD_cru_rows.mat')
% load('DD_dSM_dt_shallow_interpsm.mat')
% load('DD_NDVI_anomaly_interpsm.mat')
% load('DD_T2M_CRU_anomaly_interpsm.mat')
% load('means_ms_md.mat')
% 
% 
% cd('E:\Daten Baur\Matlab code\redblue')
% redblue_color = redblue(100) ; 
% 
% 
% 
% %% Sm and f_e anomalies
% 
% length_DD = sum(~isnan(DD_dSM_dt_shallow_interpsm),2,'omitnan') ; 
% 
% 
% 
%  imagesc(mean_ms)  ; % clim([0 400])
% [xs , ys] = getpts() ; xs = round(xs) ; ys= round(ys) ; 
% close
% %  xs = [1 720] ; ys = [1 360] ;
% 
% dummy_DD_dSM_dt = DD_dSM_dt_shallow_interpsm ; 
% dummy_DD_NDVI =   DD_NDVI_anomaly_interpsm   ;
% % dummy_DD_T2M_CRU =   DD_T2M_CRU_anomaly_interpsm   ;
% % DD_NDVI_anomaly_interpsm(DD_NDVI_anomaly_interpsm < 0.005 & DD_NDVI_anomaly_interpsm > -0.005) = NaN ; 
% 
% % kick out fast DD might be necessary
% dummy_DD_dSM_dt(dummy_DD_dSM_dt < -0.1) = NaN ; 
% % dummy_DD_dSM_dt(length_DD > 5,:) = NaN ; 
%  % dummy_DD_dSM_dt(dummy_DD_NDVI == 0) = NaN ; 
% 
% % % nother subset base don number of drydowns
% % NDry_extract = DD_NDry > 100 ; 
% % dummy_DD_dSM_dt(~NDry_extract,:) = NaN ; 
% % dummy_DD_dSM_dt_deep(~NDry_extract,:) = NaN ; 
% % dummy_DD_NDVI(~NDry_extract,:) = NaN ;
% 
% % spatial extract on subset
% spatial_extract = DD_cru_rows > ys(1) & DD_cru_rows < ys(2) & DD_cru_cols > xs(1) & DD_cru_cols < xs(2) ; 
% dummy_DD_dSM_dt(~spatial_extract,:) = NaN ; 
% % dummy_DD_dSM_dt_deep(~spatial_extract,:) = NaN ; 
% dummy_DD_NDVI(~spatial_extract,:) = NaN ; 
% % dummy_DD_T2M_CRU(~spatial_extract,:) = NaN ; 
% 
% 
% % % kick out mode SM = 1 pixels? 
% % filtermode1 = ms_mode_vec ==  1 ;
% % filtermode1 = filtermode1(:) ; 
% % col_cru_v_test = col_cru_v ; 
% % row_cru_v_test = row_cru_v ; 
% % col_cru_v_test(filtermode1) = [] ; 
% % row_cru_v_test(filtermode1) = [] ; 
% % 
% % spatial_mode1_extract = ismember([DD_cru_rows DD_cru_cols],[row_cru_v col_cru_v],'rows')  ; 
% % 
% % dummy_DD_dSM_dt(~spatial_mode1_extract,:) = NaN ; 
% % dummy_DD_NDVI(~spatial_mode1_extract,:) = NaN ; 
% % dummy_DD_T2M_CRU(~spatial_mode1_extract,:) = NaN ; 
% 
% 
% 
% mean_DD_dSM_dt = median(dummy_DD_dSM_dt,1,'omitnan') ; 
% 
% f_E_bins = linspace(-0.5,0.5,51) ; 
% sminterp_40 = linspace(0,1,51); 
% sminterp = linspace(0,1,50) ; 
% 
% DD_dSM_dt_shallow_2D = NaN(50,50,200000) ; 
% 
% 
% counter = NaN(50,50) ; 
% counter(:,:) = 1 ; 
% 
% for sm_bins = 1:49
% 
%     sm_bin_s = sminterp(sm_bins) ;
%     sm_bin_e = sminterp(sm_bins+1) ;
%     index_sm = sminterp > sm_bin_s & sminterp < sm_bin_e ; 
%     % this is if original sminterp is same size as binning here
%     index_sm = sm_bins ; 
% 
%     dSM_dt_shallow_extract = dummy_DD_dSM_dt(:,index_sm) ;  
%     f_e_extract = dummy_DD_NDVI(:,index_sm) ;  
%     % dSM_dt_deep_extract = dummy_DD_dSM_dt_deep(:,index_sm) ;     
% 
% 
%     for f_E_bin = 1:49
% 
%             f_E_bin_s = f_E_bins(f_E_bin) ;
%             f_E_bin_e = f_E_bins(f_E_bin+1) ;
% 
% 
% 
%             dSM_dt_shallow_extract_2 = dSM_dt_shallow_extract(f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e) ; 
%             dSM_dt_shallow_extract_2(isnan(dSM_dt_shallow_extract_2)) = [] ; 
% 
%             % dSM_dt_deep_extract_2 = dSM_dt_deep_extract(f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e) ; 
%             % dSM_dt_deep_extract_2(isnan(dSM_dt_deep_extract_2)) = [] ;            
% 
% 
%             DD_dSM_dt_shallow_2D(f_E_bin,sm_bins, 1:length(dSM_dt_shallow_extract_2)) = dSM_dt_shallow_extract_2 - mean_DD_dSM_dt(index_sm); 
%             % DD_dSM_dt_deep_2D(f_E_bin,sm_bins, 1:length(dSM_dt_deep_extract_2)) = dSM_dt_deep_extract_2 - mean_DD_dSM_dt(index_sm); 
% 
% 
% 
%     end
%     sm_bins
% end
% 
% 
% 
% 
% 
% 
% DD_dSM_dt_shallow_2D(DD_dSM_dt_shallow_2D == 0) = NaN ; 
% samples_3D_count = sum(~isnan(DD_dSM_dt_shallow_2D),3,'omitnan') ;  
% xmap = median(DD_dSM_dt_shallow_2D,3,'omitnan') ; 
% % xmap(samples_3D_count < 1000) = NaN ; 
% % xmap(samples_3D_count < 100) = NaN ; 
% % samples_3D_count(samples_3D_count < 100) = NaN ; 
% 
% % imagesc(samples_3D_count)
% DD_dSM_dt_shallow_2D_zeppe_global = xmap ; 
% samples_2D_count_zeppe_global = samples_3D_count ; 
% 
% save('F:\projects\SM_long_term_DDs\data_for_figures\DD_dSM_dt_shallow_2D_zeppe_global','DD_dSM_dt_shallow_2D_zeppe_global') ; 
% save('F:\projects\SM_long_term_DDs\data_for_figures\samples_2D_count_zeppe_global','samples_2D_count_zeppe_global') ; 
% 
% 
% cd('E:\Daten Baur\Matlab code\redblue')
% redblue_color = redblue(100) ; 
% 
% 
% Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
% h1 = pcolor(sminterp_40(2:end) - diff(sminterp_40(1:2))/2,   f_E_bins(2:end) - diff(f_E_bins(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
% clim([-4e-3 4e-3]) %
% % clim([-6e-4 0.0]) %
%  ylim([-0.4 0.4])
%  xlim([0 1])
%  colormap(redblue_color) %
% % colormap parula
% cbr = colorbar ;
% cbr.Label.String = "\DeltaSM/\Deltat shallow anomaly [m³/m³/day]";
% xlabel('SM [m³/m³]')
% ylabel('f_E anomaly [-]')
% title('Europe')
% set(gca,'FontSize',17)
% pbaspect([1 1 1])
% % fontsize(16,'points')
% 
% 
% saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\SM_loss_shallow_2D_SM_fE_anomalies_Europe_100samp','svg')
% close 
% 
% 
% 
% 
% test = matfile('DD_T2M_CRU_anomaly_interpsm.mat') ; 
% test = test.DD_T2M_CRU_anomaly_interpsm ;
% 
% 
% imagesc(DD_T2M_CRU_anomaly_interpsm(1:1000,:)) ; 
% figure
% imagesc(test(1:1000,:)) ; 
% 
% 
% plot(mean(DD_T2M_CRU_anomaly_interpsm,1,'omitnan'))
% hold on
% plot(mean(DD_T2M_CRU_anomaly_interpsm,1,'omitnan'))




%%  2D binning into SM and f_E values
% % do 60 by 60 reduce if needed
% 
% 
% imagesc(mean_Ts)  ; clim([0 400])
% [xs , ys] = getpts() ; xs = round(xs) ; ys= round(ys) ; 
% close
% %  xs = [1 720] ; ys([1 360]) ;
% 
% dummy_DD_dSM_dt = DD_dSM_dt_shallow_interpsm ; 
% dummy_DD_NDVI =   DD_NDVI_interpsm   ;
% spatial_extract = DD_cru_rows > ys(1) & DD_cru_rows < ys(2) & DD_cru_cols > xs(1) & DD_cru_cols < xs(2) ; 
% 
% dummy_DD_dSM_dt(~spatial_extract,:) = NaN ; 
% dummy_DD_NDVI(~spatial_extract,:) = NaN ; 
% 
% 
% 
% f_E_bins = linspace(0,1,51) ; 
% sminterp_40 = linspace(0,1,51); 
% sminterp = linspace(0,1,50) ; 
% 
% DD_dSM_dt_shallow_2D = NaN(50,50,100000) ; 
% 
% 
% 
% for sm_bins = 1:49
% 
%     sm_bin_s = sminterp(sm_bins) ;
%     sm_bin_e = sminterp(sm_bins+1) ;
%     index_sm = sminterp > sm_bin_s & sminterp < sm_bin_e ; 
%     % this is if original sminterp is same size as binning here
%     index_sm = sm_bins ; 
% 
%     dSM_dt_shallow_extract = dummy_DD_dSM_dt(:,index_sm) ;  
%     f_e_extract = dummy_DD_NDVI(:,index_sm) ;  
% 
% 
% 
%     for f_E_bin = 1:49
% 
%             f_E_bin_s = f_E_bins(f_E_bin) ;
%             f_E_bin_e = f_E_bins(f_E_bin+1) ;
% 
% 
% 
%             dSM_dt_shallow_extract_2 = dSM_dt_shallow_extract(f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e) ; 
%             dSM_dt_shallow_extract_2(isnan(dSM_dt_shallow_extract_2)) = [] ; 
% 
% 
%             DD_dSM_dt_shallow_2D(f_E_bin,sm_bins, 1:length(dSM_dt_shallow_extract_2)) = dSM_dt_shallow_extract_2 ; 
%             DD_dSM_dt_deep_2D(f_E_bin,sm_bins, 1:length(dSM_dt_deep_extract_2)) = dSM_dt_deep_extract_2 ; 
% 
% 
% 
%     end
%     sm_bins
% end
% 
% 
% 
% 
% DD_dSM_dt_shallow_2D(DD_dSM_dt_shallow_2D == 0) = NaN ; 
% samples_3D_count = sum(~isnan(DD_dSM_dt_shallow_2D),3,'omitnan') ;  
% xmap = mean(DD_dSM_dt_shallow_2D,3,'omitnan') ; 
% xmap(samples_3D_count < 20) = NaN ; 
% 
% 
% 
% Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
% h1 = pcolor(sminterp_40(2:end) - diff(sminterp_40(1:2))/2,   f_E_bins(2:end) - diff(f_E_bins(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
% clim([-0.015 0.0]) %
% % clim([-6e-4 0.0]) %
%  ylim([0 1])
%  xlim([0 1])
% % colormap(redblue_color) %
% colormap parula
% cbr = colorbar ;
% cbr.Label.String = "\DeltaSM/\Deltat shallow [m³/m³/day]";
% xlabel('SM [m³/m³]')
% ylabel('f_E [-]')
% title('India')
% set(gca,'FontSize',17)
% pbaspect([1 1 1])
% 
% 
% 
% saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model\SM_loss_shallow_2D_NDVI_forcing_India','svg')
% close 
% 
% 






%% Now f_E anomalies and temp
%  clear 
%  cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_03')
% 
% 
% load('DD_dSM_dt_shallow_interpsm.mat')
% load('DD_T2M_CRU_anomaly_interpsm.mat')
% % load('DD_f_E_anomaly_interpsm.mat')
% load('DD_NDVI_anomaly_interpsm.mat')
% load('DD_cru_rows')
% load('DD_cru_cols')
% load('means_ms_md.mat')
% 
% 
% imagesc(mode_ms)  ; %clim([0 400])
% [xs , ys] = getpts() ; xs = round(xs) ; ys= round(ys) ; 
% close
% %  xs = [1 720] ; ys = [1 360] ;
% 
% 
% % spatial extract on subset
% spatial_extract = DD_cru_rows > ys(1) & DD_cru_rows < ys(2) & DD_cru_cols > xs(1) & DD_cru_cols < xs(2) ; 
% 
% 
% f_E_bins = linspace(-0.2,0.2,51) ; 
% sminterp_40 = linspace(0,1,51); 
% % kelvin bins .. can be used for c and s
% temp_bins = linspace(-10,10,51) ; 
% 
% 
% dummy_dSM_dt = DD_dSM_dt_shallow_interpsm(:,10:20) ; 
% dummy_f_E = DD_NDVI_anomaly_interpsm(:,10:20) ; 
% dummy_T2M = DD_T2M_CRU_anomaly_interpsm(:,10:20) ; 
% 
% 
% % dummy_dSM_dt = DD_dSM_dt_shallow_interpsm; 
% % dummy_f_E = DD_f_E_anomaly_interpsm ; 
% % dummy_T2M = DD_T2M_CRU_anomaly_interpsm ; 
% % dummy_dSM_dt_mean = mean(dummy_dSM_dt(:),'omitnan') ;
% 
% 
% dummy_dSM_dt(~spatial_extract,:) = NaN ; 
% dummy_f_E(~spatial_extract,:) = NaN ; 
% dummy_T2M(~spatial_extract,:) = NaN ; 
% dummy_dSM_dt_mean = mean(dummy_dSM_dt(:),'omitnan') ;
% 
% 
% 
% 
% 
% DD_dSM_dt_shallow_2D_temp_f_E = NaN(50,50,100000) ; 
% 
% 
% for temp_bin = 1:50
%     % temp_bin = 23 ; 
% 
%     temp_bin_s = temp_bins(temp_bin) ;
%     temp_bin_e = temp_bins(temp_bin+1) ;
% 
% 
% 
%     for f_E_bin = 1:50
%    % f_E_bin = 25
% 
%             f_E_bin_s = f_E_bins(f_E_bin) ;
%             f_E_bin_e = f_E_bins(f_E_bin+1) ;
% 
%             % standard bining 
%             extract_index = dummy_f_E > f_E_bin_s & dummy_f_E < f_E_bin_e & ...
%                     dummy_T2M > temp_bin_s & dummy_T2M < temp_bin_e;
% 
% 
% 
%             dSM_dt_extract = dummy_dSM_dt(extract_index) ; 
% 
%             % get dSM/dt anomalies
%             dSM_dt_extract = dSM_dt_extract - dummy_dSM_dt_mean ; 
%             DD_dSM_dt_shallow_2D_temp_f_E(f_E_bin,temp_bin,1:length(dSM_dt_extract)) = dSM_dt_extract; 
% 
% 
% 
%     end
%     temp_bin
% end
% 
% 
% 
% 
% DD_dSM_dt_shallow_2D_temp_f_E(DD_dSM_dt_shallow_2D_temp_f_E == 0) = NaN ; 
% 
% 
% samples_3D_count = sum(~isnan(DD_dSM_dt_shallow_2D_temp_f_E),3,'omitnan') ;  
% xmap = mean(DD_dSM_dt_shallow_2D_temp_f_E,3,'omitnan') ; 
% xmap(samples_3D_count < 100) = NaN ; 
% % xmap(samples_3D_count < 1000) = NaN ; 
% 
% 
% cd('E:\Daten Baur\Matlab code\redblue')
% redblue_color = redblue(100) ; 
% 
% 
% Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
% h1 = pcolor(temp_bins(2:end) - diff(temp_bins(1:2))/2,   f_E_bins(2:end) - diff(f_E_bins(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
% clim([-5e-3 5e-3]) %
%  ylim([-0.2 0.2])
%  xlim([-20 20])
%   colormap(redblue_color) %
% %colormap parula
% cbr = colorbar ;
% cbr.Label.String = "\DeltaSM/\Deltat shallow anomaly [m³/m³/day]";
% xlabel('T forcing anomaly [°K]')
% ylabel('f_E anomaly [-]')
% title('SM = 0.2-0.3 South America')
% set(gca,'FontSize',17)
% pbaspect([1 1 1])
% yline(0)
% xline(0)
% 
% 
% pbaspect([0.8 0.8 0.8])
% axes = gca ; 
% axes.Units = 'centimeters' ; 
% axes_diff_Position = get(axes, 'Position');
% % calc positions
% %  3.6987    2.3080   22.0499   17.0999   are the pos 
%  arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',6,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
%  arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',6,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
% 
%  set(arrow1,'Position',[3.6987+22.0499+1 2.3080+(17.0999/2)+1   0  (17.0999/2)-1]) ; 
%  set(arrow2,'Position',[3.6987+22.0499+1 2.3080+(17.0999/2)-1   0 -(17.0999/2)+1]) ; 
% 
% textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
% 'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',17,'Units','centimeters');
% set(textbox1,'Position',[3.6987+22.0499+2, 2.3080+(17.0999/2)+1,    -(17.0999/2)+1, 0]) ; 
% set(gca,'Box','on');
% 
% textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
% 'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',17,'Units','centimeters');
% set(textbox2,'Position',[3.6987+22.0499+2, 2.3080+(17.0999/2)-1,    (17.0999/2)-1, 0]) ; 
% set(gca,'Box','on');
% set(gca,'FontSize',17)
% 
% 
% 
% saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\SM_loss_shallow_2D_fE_T2M_anomalies_SouthAmerica_100samp','svg')
% close 





%% SM and T2M anomalies
% clear 
% 
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests')
% 
% load('DD_cru_cols.mat')
% load('DD_cru_rows.mat')
% load('DD_dSM_dt_shallow_interpsm.mat')
% load('DD_NDVI_anomaly_interpsm.mat')
% % load('DD_T2M_CRU_anomaly_interpsm.mat')
% load('DD_T2M_ERA_anomaly_interpsm.mat')
% load('F:\projects\SM_long_term_DDs\Zeppe_model_output_02\means_global.mat')
% 
% cd('E:\Daten Baur\Matlab code\redblue')
% redblue_color = redblue(100) ; 
% 
% 
% imagesc(mean_Ts)  ; clim([0 400])
% [xs , ys] = getpts() ; xs = round(xs) ; ys= round(ys) ; 
% close
% %  xs = [1 720] ; ys = [1 360] ;
% 
% 
% dummy_DD_dSM_dt = DD_dSM_dt_shallow_interpsm ; 
% dummy_DD_NDVI =   DD_NDVI_anomaly_interpsm   ;
% dummy_DD_T2M_CRU =   DD_T2M_ERA_anomaly_interpsm   ;
% 
% 
% % % nother subset base don number of drydowns
% % NDry_extract = DD_NDry > 100 ; 
% % dummy_DD_dSM_dt(~NDry_extract,:) = NaN ; 
% % dummy_DD_dSM_dt_deep(~NDry_extract,:) = NaN ; 
% % dummy_DD_NDVI(~NDry_extract,:) = NaN ; 
% 
% % spatial extract on subset
% spatial_extract = DD_cru_rows > ys(1) & DD_cru_rows < ys(2) & DD_cru_cols > xs(1) & DD_cru_cols < xs(2) ; 
% dummy_DD_dSM_dt(~spatial_extract,:) = NaN ; 
% % dummy_DD_dSM_dt_deep(~spatial_extract,:) = NaN ; 
% dummy_DD_NDVI(~spatial_extract,:) = NaN ; 
% dummy_DD_T2M_CRU(~spatial_extract,:) = NaN ; 
% 
% 
% 
% mean_DD_dSM_dt = median(dummy_DD_dSM_dt,1,'omitnan') ; 
% 
% f_E_bins = linspace(-0.2,0.2,51) ; 
% T2M_bins = linspace(-10,10,51) ; 
% sminterp_40 = linspace(0,1,51); 
% sminterp = linspace(0,1,50) ; 
% 
% DD_dSM_dt_shallow_2D = NaN(50,50,100000) ; 
% % DD_dSM_dt_deep_2D = NaN(50,50,100000) ; 
% 
% 
% for sm_bins = 1:49
% 
%     sm_bin_s = sminterp(sm_bins) ;
%     sm_bin_e = sminterp(sm_bins+1) ;
%     index_sm = sminterp > sm_bin_s & sminterp < sm_bin_e ; 
%     % this is if original sminterp is same size as binning here
%     index_sm = sm_bins ; 
% 
%     dSM_dt_shallow_extract = dummy_DD_dSM_dt(:,index_sm) ;  
%     T2M_extract = dummy_DD_T2M_CRU(:,index_sm) ;  
%     % dSM_dt_deep_extract = dummy_DD_dSM_dt_deep(:,index_sm) ;     
% 
% 
%     for T2M_bin = 1:49
% 
%             T2M_bin_s = T2M_bins(T2M_bin) ;
%             T2M_bin_e = T2M_bins(T2M_bin+1) ;
% 
% 
% 
%             dSM_dt_shallow_extract_2 = dSM_dt_shallow_extract(T2M_extract > T2M_bin_s & T2M_extract < T2M_bin_e) ; 
%             dSM_dt_shallow_extract_2(isnan(dSM_dt_shallow_extract_2)) = [] ; 
% 
%             % dSM_dt_deep_extract_2 = dSM_dt_deep_extract(f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e) ; 
%             % dSM_dt_deep_extract_2(isnan(dSM_dt_deep_extract_2)) = [] ;            
% 
% 
%             DD_dSM_dt_shallow_2D(T2M_bin,sm_bins, 1:length(dSM_dt_shallow_extract_2)) = dSM_dt_shallow_extract_2 - mean_DD_dSM_dt(index_sm); 
%             % DD_dSM_dt_deep_2D(f_E_bin,sm_bins, 1:length(dSM_dt_deep_extract_2)) = dSM_dt_deep_extract_2 - mean_DD_dSM_dt(index_sm); 
% 
% 
% 
%     end
%     sm_bins
% end
% 
% 
% 
% 
% 
% 
% DD_dSM_dt_shallow_2D(DD_dSM_dt_shallow_2D == 0) = NaN ; 
% samples_3D_count = sum(~isnan(DD_dSM_dt_shallow_2D),3,'omitnan') ;  
% xmap = median(DD_dSM_dt_shallow_2D,3,'omitnan') ; 
% % xmap(samples_3D_count < 1000) = NaN ; 
% % xmap(samples_3D_count < 100) = NaN ; 
% 
% 
% % % imagesc(samples_3D_count)
% %  DD_dSM_dt_shallow_T2M_2D_zeppe_global = xmap ; 
% %  samples_3D_count_T2M_zeppe_global = samples_3D_count ; 
% % % 
% %  save('F:\projects\SM_long_term_DDs\data_for_figures\DD_dSM_dt_shallow_T2M_2D_zeppe_global','DD_dSM_dt_shallow_T2M_2D_zeppe_global') ; 
% %  save('F:\projects\SM_long_term_DDs\data_for_figures\samples_3D_count_T2M_zeppe_global','samples_3D_count_T2M_zeppe_global') ; 
% 
% 
% 
% Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
% h1 = pcolor(sminterp_40(2:end) - diff(sminterp_40(1:2))/2,   T2M_bins(2:end) - diff(T2M_bins(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
% clim([-1e-2 1e-2]) %
% % clim([-6e-4 0.0]) %
%  ylim([-15 15])
%  xlim([0 1])
%  colormap(redblue_color) %
% % colormap parula
% cbr = colorbar ;
% cbr.Label.String = "\DeltaSM/\Deltat shallow anomaly [m³/m³/day]";
% xlabel('SM [m³/m³]')
% ylabel('T2M ERA anomaly [-]')
% title('Europe')
% set(gca,'FontSize',17)
% pbaspect([1 1 1])
% % fontsize(16,'points')
% 
% 
% saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\SM_loss_shallow_2D_SM_fE_anomalies_Europe_100samp','svg')
% close 
% 







%% Now NDVI and temp fix SM 
% 
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_03')
% 
% load('DD_tmp_CRUJRA_interpsm.mat')
% load('DD_NDVI_interpsm.mat')
% load('DD_dSM_dt_shallow_interpsm.mat')
% 
% 
% f_E_bins = linspace(0,1,51) ; 
% sminterp_40 = linspace(0,1,51); 
% % kelvin bins .. can be used for c and s
% temp_bins = linspace(260,320,51) ; 
% 
% 
% dummy_dSM_dt = DD_dSM_dt_shallow_interpsm(:,30:40) ; 
% dummy_f_E = DD_NDVI_interpsm(:,30:40) ; 
% dummy_T2M = DD_tmp_CRUJRA_interpsm(:,30:40) ; 
% dummy_dSM_dt_mean = median(dummy_dSM_dt(:),'omitnan') ;
% 
% 
% DD_dSM_dt_shallow_2D_temp_f_E = NaN(40,40,100000) ; 
% 
% 
% for temp_bin = 1:50
% 
%     temp_bin_s = temp_bins(temp_bin) ;
%     temp_bin_e = temp_bins(temp_bin+1) ;
% 
% 
% 
%     for f_E_bin = 1:50
% 
%             f_E_bin_s = f_E_bins(f_E_bin) ;
%             f_E_bin_e = f_E_bins(f_E_bin+1) ;
% 
%             % standard bining 
%             extract_index = dummy_f_E > f_E_bin_s & dummy_f_E < f_E_bin_e & ...
%                     dummy_T2M > temp_bin_s & dummy_T2M < temp_bin_e;
% 
% 
% 
%             dSM_dt_extract = dummy_dSM_dt(extract_index) ; 
% 
%             % get dSM/dt anomalies
%             dSM_dt_extract = dSM_dt_extract - dummy_dSM_dt_mean ; 
%             DD_dSM_dt_shallow_2D_temp_f_E(f_E_bin,temp_bin,1:length(dSM_dt_extract)) = dSM_dt_extract; 
% 
% 
%     end
%     temp_bin
% end
% 
% 
% 
% 
% 
% DD_dSM_dt_shallow_2D_temp_f_E(DD_dSM_dt_shallow_2D_temp_f_E == 0) = NaN ; 
% samples_3D_count = sum(~isnan(DD_dSM_dt_shallow_2D_temp_f_E),3,'omitnan') ;  
% 
% xmap = median(DD_dSM_dt_shallow_2D_temp_f_E,3,'omitnan') ; 
% xmap(samples_3D_count < 100) = NaN ; 
% 
% 
% cd('E:\Daten Baur\Matlab code\redblue')
% redblue_color = redblue(100) ; 
% 
% 
% Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
% h1 = pcolor(temp_bins(2:end) - diff(temp_bins(1:2))/2,   f_E_bins(2:end) - diff(f_E_bins(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
% clim([-8e-3 8e-3]) %
%  ylim([0 1])
%  xlim([260 320])
%   colormap(redblue_color) %
% %colormap parula
% cbr = colorbar ;
% cbr.Label.String = "\DeltaSM/\Deltat shallow anomaly [m³/m³/day]";
% xlabel('T2M [°K]')
% ylabel('NDVI [-]')
% title('SM = 0.2-0.3')
% set(gca,'FontSize',17)
% pbaspect([1 1 1])
% 
% 
% 
% pbaspect([0.8 0.8 0.8])
% axes = gca ; 
% axes.Units = 'centimeters' ; 
% axes_diff_Position = get(axes, 'Position');
% % calc positions
% %  3.6987    2.3080   22.0499   17.0999   are the pos 
%  arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',6,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
%  arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',6,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
% 
%  set(arrow1,'Position',[3.6987+22.0499+1 2.3080+(17.0999/2)+1   0  (17.0999/2)-1]) ; 
%  set(arrow2,'Position',[3.6987+22.0499+1 2.3080+(17.0999/2)-1   0 -(17.0999/2)+1]) ; 
% 
% textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
% 'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',17,'Units','centimeters');
% set(textbox1,'Position',[3.6987+22.0499+2, 2.3080+(17.0999/2)+1,    -(17.0999/2)+1, 0]) ; 
% set(gca,'Box','on');
% 
% textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
% 'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',17,'Units','centimeters');
% set(textbox2,'Position',[3.6987+22.0499+2, 2.3080+(17.0999/2)-1,    (17.0999/2)-1, 0]) ; 
% set(gca,'Box','on');
% set(gca,'FontSize',17)
% 
% 
% 
% 
% saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\SM_loss_2D_NDVI_T2M_SM_02_03','svg')
% close 
% 
% 











%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze spatial and temporal distrubtiono of drydown frequency model vs
% observation
% clear 
% 
% cd('F:\ESA_CCI\2D_02')
% load('t_start_end_2D_array.mat')
% load('rowcol_2D_array.mat')
% % load('dSM_dt_interpsm_2D_array.mat')
% 
% ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
% 
% 
% mask = t_start_end_2D_array(:,1) < 4080 | t_start_end_2D_array(:,2) > 13210;
% rowcol_2D_array(mask,:) = [] ; 
% t_start_end_2D_array(mask,:) = [] ; 
% 
% 
% t_start_end_2D_array = t_start_end_2D_array - 4080 ; 
% % max(t_start_end_2D_array)  ; max(DD_timev)
% 
% 
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_03')
% load('DD_timev.mat')
% load('DD_cru_rows.mat')
% load('DD_cru_cols.mat')
% 
% 
% % do some kinda loop through unique cols and rows, find the matching SMAP
% % pixels and then compare the temporal distribution of drdowns
% 
% cd('F:\AVHRR_phenology')
% 
% 






%% MJB 23.07.2024 Calculate trends analysis for model results

% 
% 
% clear
% 
% cd('F:\ESA_CCI\analysis_datasets\spatial_02')
% % dd data
% dSM_dt_files = string(ls('*dSM_dt_interpsm*')) ;
% rowcol_files = string(ls('*rowcol_year_dummy*')) ;
% t_start_end_files = string(ls('*t_start*')) ;
% 
% % monthly LAI will be only until 2014
% cd('F:\ESA_CCI\analysis_datasets\NDVI_sminterp')
% NDVI_files_all = string(ls('*VIP_NDVI_daily*')) ;
% NDVI_files_121 = string(ls('*VIP_NDVI_daily_sminterp_121_*')) ;
% NDVI_files_5y = string(ls('*VIP_NDVI_daily_sminterp_anomaly*')) ;
% 
% NDVI_files_all = strtrim(NDVI_files_all) ;
% NDVI_files_121 = strtrim(NDVI_files_121) ;
% NDVI_files_5y = strtrim(NDVI_files_5y) ;
% 
% 
% [Lia Locb] = ismember(NDVI_files_121, NDVI_files_all,'rows') ; 
% NDVI_files = NDVI_files_all ; 
% NDVI_files(Locb) = [] ; 
% 
% 
% % find spatial locations of pos neg 
% load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_pos')
% load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_neg')
% load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_all')
% 
% load('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_mean_NDVI_FYslope.mat') 
% 
% % dSM_dt_mean_NDVI_pos_mask = dSM_dt_mean_NDVI_pos > 0 ; 
% 
% 
% 
% 
% ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
% tt = 1:15402 ; 
% 
% NDVI_binning = linspace(-0.20,0.20,41) ; 
% dSM_dt_binning = linspace(-0.08,0.08,41) ; 
% 
% sxity_to_thirty = 0:2:60 ; sxity_to_thirty(1) = 1 ; 
% 
% counter = NaN(40,40) ; 
% counter(:,:) = 1 ; 
% 
% 
% % ESA_CCI_datetime(4080)
% % ESA_CCI_datetime(end)
% years_vec = 1990:2014 ; 
% for i = 1:length(years_vec)
%     datetime_dummy_years(i) = datetime(strcat('01-Jan-',num2str(years_vec(i)))) ;
% end
% % datetime_dummy_years(end) = ESA_CCI_datetime(end) ; 
%  datetime_dummy_years(end) = datetime('31-Dec-2014') ; 
% 
% 
% 
% % spatial = 2450
% % spatial = 2550
% cd('F:\ESA_CCI\analysis_datasets')
% 
% 
% 
% 
% 
% for years_index = 10:length(datetime_dummy_years)
% 
% 
% counter = NaN(40,40) ; 
% counter(:,:) = 1 ; 
% % binning based on NDVI
% % NDVI_dSM_dt_binning = NaN(60,60,500000) ; 
% 
% 
% year_vector_start = datetime_dummy_years(years_index) ; 
% year_vector_end = datetime_dummy_years(years_index+1) ; 
% 
% 
% % binning based on NDVI
% dSM_dt_NDVI_binning = NaN(40,40,300000) ; 
% 
% 
% 
% 
% 
% % spatial = 2450
% % spatial = 2550
% 
% 
% 
% for spatial = 1:length(dSM_dt_files)
% 
% 
% % get drydown data
% 
% dummy_dSM_dt = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',dSM_dt_files(spatial))) ; 
% NDVI_dummy =  matfile(strcat('F:\ESA_CCI\analysis_datasets\NDVI_sminterp\',NDVI_files_5y(spatial))) ; 
% 
% % get 2.5 degree rowcol from name
% name_files = char(dSM_dt_files(spatial)) ; 
% % grep row and col
% pat = digitsPattern;
% row_col_name = extract(name_files,pat) ; 
% row_2_5 = str2double(row_col_name{1}) ; 
% col_2_5 = str2double(row_col_name{2}) ;     
% 
% 
% % logical condition to only analyze 2,5 degree boxes with extraction of
% % cover effect
%  % if (dSM_dt_mean_NDVI_FYslope(row_2_5,col_2_5) > 0)
%  %     continue
%  % end
% 
% 
% 
% 
% dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy   ; 
% % dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;
% 
% % less than 3k drydowns remove
% if size(dummy_dSM_dt,1) < 3000
%     continue
% end
% 
% 
% 
% % new condition removing very short DDs?
% % dummy_dSM_dt(sum(~isnan(dummy_dSM_dt),2) < 3 , :) = NaN ; 
% % kick out all drainage to save memory. reduced to 40 length from 60
% 
% 
% dummy_dSM_dt(:,41:60) = [] ;
% 
% 
%  % dummy_NDVI = NDVI_dummy.dummy_NDVI_monthly ; 
%  dummy_NDVI = NDVI_dummy.NDVI_spatial_sminterp ; 
%  dummy_NDVI(:,41:60) = [] ;
% 
% 
% % exclude dainage
% %dummy_dSM_dt(:,40:end) = NaN ; 
% 
% 
% dummy_rowcol = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',rowcol_files(spatial))) ; 
% dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
% 
% dummy_t_start_end = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',t_start_end_files(spatial))) ; 
% dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
% dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 
% 
% % pre filter for time if needed pre 1990 and post 2014
% mask = dummy_t_start_end(:,2) < 4080 |  dummy_t_start_end(:,1) > 13210; % 1-Jan-1990
% 
% 
% dummy_t_mean(mask) = NaN ; 
% dummy_dSM_dt(mask,:) = NaN ; 
% dummy_rowcol(mask,:) = NaN ; 
% dummy_NDVI(mask,:) = NaN ; 
% 
% % convert to mm/days loss
% dummy_dSM_dt = dummy_dSM_dt .* -100  ; 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % temporal condition to 1 year
% % time_start = find(ESA_CCI_datetime == datetime('01-Jan-2003')) ; 
% % time_end = find(ESA_CCI_datetime == datetime('01-Jan-2004')) ; 
% time_start = find(ESA_CCI_datetime == year_vector_start) ; 
% time_end = find(ESA_CCI_datetime == year_vector_end) ; 
% 
% 
% 
% mask = dummy_t_mean < time_start |  dummy_t_mean > time_end ; % 1-Jan-1990
% dummy_t_mean(mask) = NaN ; 
% dummy_dSM_dt(mask,:) = NaN ; 
% dummy_rowcol(mask,:) = NaN ; 
% dummy_NDVI(mask,:) = NaN ; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% % binning based on NDVI
% % dSM_dt_NDVI_binning = NaN(30,60) ; 
% %NDVI_mean = mean(dummy_NDVI,'omitnan') ; 
% NDVI_deviations = dummy_NDVI  ; 
% % get mean but for each SM bin
% dummy_dSM_dt_mean = median(dummy_dSM_dt,1,'omitnan') ; 
% 
% for bin = 1:40
% 
%     cur_bin_min = NDVI_binning(bin) ;
%     cur_bin_max = NDVI_binning(bin+1) ;   
%     extract_index = NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max ;
%     extract_index_row =  any(NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max,2) ; 
%     % get dSM/dt anomalies
%     extract_index = extract_index(extract_index_row,:) ; 
%     dSM_dt_extract = dummy_dSM_dt(extract_index_row,:) - dummy_dSM_dt_mean ; 
%     dSM_dt_extract(~extract_index) = NaN ; 
% 
%      % now do propper 2D binning with i and j to use memory efficiently.
%      % Binning array is filled with 3rd dim index specific for row and col
%     for j = 1:size(dSM_dt_extract,2)
% 
%         if all(isnan( dSM_dt_extract(:,j)))
%             continue
%         end
% 
%         for i =1:size(dSM_dt_extract,1)
% 
%         if (isnan( dSM_dt_extract(i,j)))
%             continue
%         end
% 
%     dSM_dt_NDVI_binning(bin,j,counter(bin,j)) = dSM_dt_extract(i,j); 
% 
%     counter(bin,j) = counter(bin,j) + 1 ;
% 
% 
%     if counter(bin,j) == 300000
%     return
%     end
% 
%         end % i
% 
%     end  % j 
% 
% 
% end
% 
%  spatial
% 
% 
% end
% 
% 
% 
% 
% % do binning and remove low sampling
% 
% dSM_dt_NDVI_binning(dSM_dt_NDVI_binning == 0) = NaN ; 
% samples_3D_count = sum(~isnan(dSM_dt_NDVI_binning),3,'omitnan') ;  
% 
% 
% xmap = median(dSM_dt_NDVI_binning,3,'omitnan') ; 
% % xmap(samples_3D_count < 100) = NaN ; 
%   % figure
%   % pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2,  (xmap)) ; colormap(redblue_color) ; clim([-0.5 0.5])  
% 
% save(strcat('F:\ESA_CCI\analysis_datasets\xmap_dSm_dt_03\dSM_dt_dev_',num2str(year(year_vector_start))),'xmap','samples_3D_count')
% 
% year_vector_start
% 
% 
% 
% end
% 
% 

















