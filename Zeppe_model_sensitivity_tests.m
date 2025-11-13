%% MJB 07-Aug-2025 Sensitivity tests on the model. Reply to reviewer comments. 
% set up model with data used in the paper but only run few pixels with a
% range of sensitvity options.




%% DATA
clear

% test = ncinfo("crujra.v2.4.5d.spfh.1990.365d.noc.nc") ; 
testinfo = ncinfo("F:\CRUJRA\crujra.v2.4.5d.tmp.2022.365d.noc.nc") ; 

cd('F:\CRUJRA\processed')

load('crujra_datetime_vector.mat')
load('datetime_1990_2014.mat')
load('row_cru_v.mat')
load('col_cru_v.mat')
load('lat_cru_v.mat')
load('lon_cru_v.mat')



SRB_file = matfile("Clara_SRB_full_array_linear.mat") ; 
tmp_file = matfile("tmp_full_array.mat") ;
pre_file = matfile("pre_full_array.mat") ;
spfh_file = matfile("spfh_full_array.mat") ;

% get VIP monthly NDVI for f_E
NDVI_file = matfile('F:\VIP_NDVI\VIP_NDVI_array_linear.mat') ; 


%tmp from ERA 5 monthly
cd('F:\ERA_5_monthly')
tmp_ERA_file = matfile("T2M_ERA_daily_interp.mat") ;


%% build 2D for pixel selection

load('F:\CRUJRA\Ancil_data\CRUJRA_land_mask.mat') ;
load('F:\CRUJRA\Ancil_data\CRUJRA_land_mask_vec.mat') ;

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d')
load('DD_dSM_dt_NDVI_gradient_zeppe_ERA.mat')


% figure
% imagesc(DD_dSM_dt_NDVI_gradient_zeppe_ERA)
% clim([-0.05 0.05]) ; 




figure
imagesc(CRUJRA_land_mask)
[xs, ys] = getpts() ; xs = round(xs) ; ys = round(ys) ; 
close




%% MODEL

Number_pixels = length(xs) ; 
Number_pixels_square = (xs(2)-xs(1)) .* (ys(2)-ys(1)) ; 


cd('F:\Matlab_Code\Zeppe_model')
% shorter with 1000 steps spinup
 Ts            = NaN(Number_pixels,9131+1000) ; 
 Tc            = NaN(Number_pixels,9131+1000) ;
 Tra_s         = NaN(Number_pixels,9131+1000) ;
 Tra_r         = NaN(Number_pixels,9131+1000) ;
 ms            = NaN(Number_pixels,9131+1000) ;
 mD            = NaN(Number_pixels,9131+1000) ;
 f_E           = NaN(Number_pixels,9131+1000) ;
 E_s           = NaN(Number_pixels,9131+1000) ;
 Cap_flux      = NaN(Number_pixels,9131+1000) ;

tmp_ERA_save   = NaN(Number_pixels,9131)      ;

steps_per_day = 1 ;

cd('F:\Matlab_Code\Zeppe_model')


% find i based on location
for i = 1:length(xs) 
[VIndex_rowcol(i)] = find(col_cru_v == xs(i) & row_cru_v == ys(i)) ; 
end


% find area of square
VIndex_rowcol = find((col_cru_v >= xs(1) & col_cru_v <= xs(2)) & (row_cru_v >= ys(1) & ...
    row_cru_v <= ys(2))) ; 
% VIndex_rowcol = sort(VIndex_rowcol) ; 


counter = 1 ; 



% i = VIndex_rowcol(1) ; 

for i = 1:length(VIndex_rowcol)

    cur_in = VIndex_rowcol(i) ; 

%  i = spatial_extract    
srb_test = SRB_file.Clara_SRB_full_array_linear(cur_in,:)        ;
tmp_test = tmp_file.tmp_full_array(cur_in,:)                     ;
pre_test = pre_file.pre_full_array(cur_in,:)                     ;
spfh_test = spfh_file.spfh_full_array(cur_in,:)                  ; 
NDVI_test = NDVI_file.VIP_NDVI_array(cur_in,:)                   ;
tmp_ERA_test = tmp_ERA_file.T2M_ERA_daily_interp(cur_in,:)       ;           
tmp_ERA_test = tmp_ERA_test + 273.15                        ; 
tmp_ERA_save(counter,:) = tmp_ERA_test                            ; 

 if mean(NDVI_test,'omitnan') < 0  || range(NDVI_test) < 0.1 || all(isnan(NDVI_test))
     continue
 end


srb_test =      cat(2,srb_test(1:1000),srb_test) ; 
tmp_test =     cat(2,tmp_test(1:1000),tmp_test) ; 
tmp_ERA_test =  cat(2,tmp_ERA_test(1:1000),tmp_ERA_test) ; 
spfh_test =     cat(2,spfh_test(1:1000),spfh_test) ;  
pre_test =      cat(2,pre_test(1:1000),pre_test) ; 
NDVI_test =     cat(2,NDVI_test(1:1000),NDVI_test) ; 



 [Ts(counter,:),Tc(counter,:),ms(counter,:),mD(counter,:),Tra_s(counter,:),Tra_r(counter,:),E_s(counter,:),Cap_flux(counter,:), f_E(counter,:)] = ...
         The_Model_f_E_vec(srb_test,tmp_ERA_test ,spfh_test,pre_test./86400,NDVI_test,steps_per_day) ; 


counter = counter + 1 ; 

counter
end


% function [Ts, TC, ms, mD, Tra_s_vec, Tra_r_vec,Es_vec, Cap_flux_vec, f_E] = The_Model_f_E_vec(F, T_R, q_R, P,NDVI,steps_per_day)
%                                % radiaiton Temp humidity precip              


% Ts = Ts(:,1001:end) ;
ms = ms(:,1001:end) ;
mD = mD(:,1001:end) ;
% rescale to same volumetric units this is critical and should rescale all
% SM values to roughly the same scale as SMAP. We can also it to 0-0.6?
% This would make it basically the same as SMAP might be useful
theta_max = 0.6 ; % pore space
ms = ms .* theta_max ;
mD = mD .* theta_max ;
tmp_ERA = tmp_ERA_test(:,1001:end) ; 




%% calculate anomalies for f_e, NDVI and temperature
load('F:\CRUJRA\processed\datetime_1990_2014.mat')



[y_NDVI, m_NDVI, d_NDVI] = ymd(datetime_1990_2014) ;

f_E_NDVI_daily_anomaly = NaN(size(f_E)) ; 
tmp_ERA_daily_anomaly =  NaN(size(tmp_ERA)) ; 

for i = 1:size(f_E,1)

    f_E_dummy = f_E(i,:) ; 
    tmp_ERA_dummy = tmp_ERA_save(i,:) ; 

        VIP_NDVI_5y_movmean = movmean(f_E_dummy,1825) ; 
        tmp_ERA_5y_movmean =  movmean(tmp_ERA_dummy,1825) ; 


        % % remove from mean year
        NDVI_extract_day_mean = NaN(365,1) ; 
        tmp_extract_day_mean = NaN(365,1) ; 

        for j = 1:365
            [y, m, d] = ymd(datetime_1990_2014(j)) ; 
            extract_index = find(d == d_NDVI & m == m_NDVI) ; 

            NDVI_extract_day_mean(j) = mean(f_E_dummy(extract_index),'omitnan') ; 
            tmp_extract_day_mean(j) = mean(tmp_ERA_dummy(extract_index),'omitnan') ;             

        end
        NDVI_extract_day_mean = repmat(NDVI_extract_day_mean,[34,1]) ; 
        NDVI_extract_day_mean = NDVI_extract_day_mean(1:length(f_E_dummy)) ; 

        tmp_extract_day_mean = repmat(tmp_extract_day_mean,[34,1]) ; 
        tmp_extract_day_mean = tmp_extract_day_mean(1:length(tmp_ERA_dummy)) ;         


        f_E_NDVI_daily_anomaly(i,:) = f_E_dummy - NDVI_extract_day_mean' ; 
        tmp_ERA_daily_anomaly(i,:) =  tmp_ERA_dummy - tmp_extract_day_mean' ;         

i
end





%% Use Model output to calculate strenght of extraction and cover effect
% do relative to varying variables of gs, root ditribution and NDVI_scaling
% basically sensitivity of the model to these variables
load('F:\CRUJRA\processed\datetime_1990_2014.mat')
sminterp = linspace(0,0.6,60) ; 


DD_dSM_dt_s_interpsm = NaN(400000,60) ; 
DD_tmp_ERA_anomaly_interpsm = NaN(400000,60) ; 
DD_NDVI_anomaly_interpsm = NaN(400000,60) ; 

DD_rows  =  NaN(400000,1) ; 
DD_cols  =  NaN(400000,1) ; 
DD_timev =  NaN(400000,2) ; 

cd('E:\Daten Baur\Matlab code\Project IGARSS multi frq tau\noodles_L_C_X')

 
DDLength = 4 ; 
NDry_count = 0 ; 
counter = 1 ; 

for i = 1:size(ms,1)

time = 1:9131 ; 
ms_dummy = ms(i,:) ; 

 row_dummy = row_cru_v(i,:) ; 
 col_dummy = col_cru_v(i,:) ; 
 mD_dummy = mD(i,:) ;
 mS_dummy = ms(i,:) ;
 NDVI_anomaly_dummy = f_E_NDVI_daily_anomaly(i,:) ; 
 tmp_anomaly_dummy = tmp_ERA_daily_anomaly(i,:) ; 


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

     DD_tt = timev{j,1} ; 
     DD_dSM_dt = diff(DD_SM) ; 

     % only use strictly decreasing dds
     if any(DD_dSM_dt > 0)
         continue
     end


     DD_SM_short = DD_SM(1:end-1) ;  

     if numel(DD_SM_short)~=numel(unique(DD_SM_short))
         continue
     end


     DD_dSM_dt_s_interpsm(counter,:) = interp1(DD_SM_short,DD_dSM_dt,sminterp,'linear',NaN) ;  

     DD_NDVI_anomaly = NDVI_anomaly_dummy(DD_tt(1:end-1)) ; 
     DD_NDVI_anomaly_interpsm(counter,:) = interp1(DD_SM_short,DD_NDVI_anomaly,sminterp,'linear',NaN) ;  

     tmp_anomaly = tmp_anomaly_dummy(DD_tt(1:end-1)) ; 
     DD_tmp_ERA_anomaly_interpsm(counter,:) = interp1(DD_SM_short,tmp_anomaly,sminterp,'linear',NaN) ;  

     DD_rows(counter) = row_dummy ; 
     DD_cols(counter) = col_dummy ;
     DD_timev(counter,1) = DD_tt(1) ;
     DD_timev(counter,2) = DD_tt(end) ;      


     % if too short in SM space skip and not save
     if all(isnan(DD_dSM_dt_s_interpsm(counter,:)))
        continue
     end

     counter = counter +1 ; 

end
i

end



DD_dSM_dt_s_interpsm(isnan(DD_rows),:) = [] ; 
DD_NDVI_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_tmp_ERA_anomaly_interpsm(isnan(DD_rows),:) = [] ; 


DD_cols(isnan(DD_rows),:) = [] ; 
DD_timev(isnan(DD_rows),:) = [] ; 
DD_rows(isnan(DD_rows),:) = [] ; 


% to mm/day
DD_dSM_dt_s_interpsm = DD_dSM_dt_s_interpsm.* -100 ; 



%% now bin drydown oriented variabiles relative to SM and anomalies to calculate extraction and cover effect
% do binning into 20 by 20 like for the 2.5 degree bins as we have slightly
% lower data depending on how many pixels we use.

DD_dSM_dt_s_interpsm(isnan(DD_rows),:) = [] ; 
DD_NDVI_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_tmp_ERA_anomaly_interpsm(isnan(DD_rows),:) = [] ; 


    f_E_bins = linspace(0,1,21) ; 
    SM_bins  = linspace(0,0.6,21) ; 
    NDVI_anomaly_bins = linspace(-0.2,0.2,21) ; 
    T2M_anomaly_bins = linspace(-10,10,21) ;    


DD_dSM_dt_s_interpsm_mean = mean(DD_dSM_dt_s_interpsm,1,'omitnan') ; 

dSM_dt_NDVI_binned_sampling = NaN(20,20) ; 
dSM_dt_NDVI_binned = NaN(20,20) ; 


for sm_bin = 1:20

    sm_bin_s = SM_bins(sm_bin) ;
    sm_bin_e = SM_bins(sm_bin+1) ;
    index_sm = sminterp > sm_bin_s & sminterp < sm_bin_e ; 

    % assign data
    dSM_dt_s_anomaly_extract = DD_dSM_dt_s_interpsm(:,index_sm)   - DD_dSM_dt_s_interpsm_mean(index_sm)  ;  
    NDVI_anomaly_extract  = DD_NDVI_anomaly_interpsm(:,index_sm) ;  
    tmp_ERA_anomaly_extract  = DD_tmp_ERA_anomaly_interpsm(:,index_sm) ;  


    for f_E_bin = 1:20

            NDVI_anomaly_bin_s = NDVI_anomaly_bins(f_E_bin) ;
            NDVI_anomaly_bin_e = NDVI_anomaly_bins(f_E_bin+1) ;

            extract_index = NDVI_anomaly_extract > NDVI_anomaly_bin_s & NDVI_anomaly_extract < NDVI_anomaly_bin_e ;         

            dSM_dt_extract = dSM_dt_s_anomaly_extract(extract_index) ; 

            dSM_dt_NDVI_binned_sampling(f_E_bin,sm_bin) = sum(~isnan(dSM_dt_extract)) ; 
            dSM_dt_NDVI_binned(f_E_bin,sm_bin) = median(dSM_dt_extract); 


    end

sm_bin
end



%% prelim figure for tests
SM_bins_02 = SM_bins(1:end-1) ; 
SM_bins_02 = SM_bins_02 + 0.0300 ; 

NDVI_anomaly_bins_02 = NDVI_anomaly_bins(1:end-1) ; 
NDVI_anomaly_bins_02 = NDVI_anomaly_bins_02 + 0.0200 ; 

cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 

dSM_dt_NDVI_binned(dSM_dt_NDVI_binned_sampling < 10) = NaN ;

xmap = dSM_dt_NDVI_binned ;
Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(SM_bins_02,   NDVI_anomaly_bins_02, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.5 0.5]) %
ylim([-0.25 0.25])
 xlim([0 0.5])
 colormap(redblue_color) %
% colormap parula
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat shallow anomaly [m³/m³/day]";
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
title('Europe')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% fontsize(16,'points')


set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 



[FX_NDVI, FY_NDVI] =  gradient(dSM_dt_NDVI_binned) ; 

dSM_dt_NDVI_binned_gradient(i) = median(FY_NDVI(:),'omitnan') ; 













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run model with varying input parameters, namely combinations of surface to vegetaiton conductance
clear






% test = ncinfo("crujra.v2.4.5d.spfh.1990.365d.noc.nc") ; 
testinfo = ncinfo("F:\CRUJRA\crujra.v2.4.5d.tmp.2022.365d.noc.nc") ; 

cd('F:\CRUJRA\processed')

load('crujra_datetime_vector.mat')
load('datetime_1990_2014.mat')
load('row_cru_v.mat')
load('col_cru_v.mat')
load('lat_cru_v.mat')
load('lon_cru_v.mat')

SRB_file = matfile("Clara_SRB_full_array_linear.mat") ; 
tmp_file = matfile("tmp_full_array.mat") ;
pre_file = matfile("pre_full_array.mat") ;
spfh_file = matfile("spfh_full_array.mat") ;

% get VIP monthly NDVI for f_E
NDVI_file = matfile('F:\VIP_NDVI\VIP_NDVI_array_linear.mat') ; 


%tmp from ERA 5 monthly
cd('F:\ERA_5_monthly')
tmp_ERA_file = matfile("T2M_ERA_daily_interp.mat") ;

% build 2D for pixel selection
load('F:\CRUJRA\Ancil_data\CRUJRA_land_mask.mat') ;
load('F:\CRUJRA\Ancil_data\CRUJRA_land_mask_vec.mat') ;

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d')
load('DD_dSM_dt_NDVI_gradient_zeppe_ERA.mat')


% figure
% imagesc(CRUJRA_land_mask)
% [xs, ys] = getpts() ; xs = round(xs) ; ys = round(ys) ; 
% close


% pixels used so far
% ys = [52 57] ; xs = [590 596] ; % boreal
% ys = [210 215] ; xs = [264 271] ; % Brazil
% ys = [139 145] ; xs = [512 520] ; % India
% ys = [156 160] ; xs = [388 395] ; % Sahel
% ys = [88 100] ; xs = [115 124] ; % US


% MODEL

Number_pixels = length(xs) ; 
Number_pixels_square = (xs(2)-xs(1)) .* (ys(2)-ys(1)) ; 


cd('F:\Matlab_Code\Zeppe_model')
% shorter with 1000 steps spinup
 Ts            = NaN(Number_pixels,9131+1000) ; 
 Tc            = NaN(Number_pixels,9131+1000) ;
 Tra_s         = NaN(Number_pixels,9131+1000) ;
 Tra_r         = NaN(Number_pixels,9131+1000) ;
 ms            = NaN(Number_pixels,9131+1000) ;
 mD            = NaN(Number_pixels,9131+1000) ;
 f_E           = NaN(Number_pixels,9131+1000) ;
 E_s           = NaN(Number_pixels,9131+1000) ;
 Cap_flux      = NaN(Number_pixels,9131+1000) ;

tmp_ERA_save   = NaN(Number_pixels,9131)      ;

steps_per_day = 1 ;

cd('F:\Matlab_Code\Zeppe_model')


% find i based on location
for i = 1:length(xs) 
[VIndex_rowcol(i)] = find(col_cru_v == xs(i) & row_cru_v == ys(i)) ; 
end


% find area of square
VIndex_rowcol = find((col_cru_v >= xs(1) & col_cru_v <= xs(2)) & (row_cru_v >= ys(1) & ...
    row_cru_v <= ys(2))) ; 
% VIndex_rowcol = sort(VIndex_rowcol) ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT array 

dSM_dt_gradient_sensitivity_gs_gc = NaN(30,30) ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameterization for sensitvity parameters
     % g_s = 1/2000;       % dry surface conductance [m/s] ORIGINAL
     % g_C = 1/800;

g_s_in_array = repmat(linspace(1/3000,1/100,30),[30,1]) ; 
g_C_in_array = repmat(linspace(1/3000,1/100,30)',[1,30]) ; 

g_s_in_array2D = repmat(linspace(1/3000,1/100,30),[30,1]) ; 
g_C_in_array2D = repmat(linspace(1/3000,1/100,30)',[1,30]) ; 


% test = (g_C_in_array./g_s_in_array) 

g_s_in_array = g_s_in_array(:) ; 
g_C_in_array = g_C_in_array(:) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_E_max = 0.8 ; 
root_frac_in = 0.7 ; 

% load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_gradient_sensitivity_gs_gc_US')
% dSM_dt_gradient_sensitivity_gs_gc = dSM_dt_gradient_sensitivity_gs_gc_US ;

% stopped US 751
for k = 751:900

cd('F:\Matlab_Code\Zeppe_model')
rootfrac = 0.7 ; 
g_s_in_dummy = g_s_in_array(k) ;
g_C_in_dummy = g_C_in_array(k) ;

counter = 1 ; 
% i = VIndex_rowcol(1) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Number_pixels = size(VIndex_rowcol,1) ;
Ts            = NaN(Number_pixels,9131+1000) ; 
 Tc            = NaN(Number_pixels,9131+1000) ;
 Tra_s         = NaN(Number_pixels,9131+1000) ;
 Tra_r         = NaN(Number_pixels,9131+1000) ;
 ms            = NaN(Number_pixels,9131+1000) ;
 mD            = NaN(Number_pixels,9131+1000) ;
 f_E           = NaN(Number_pixels,9131+1000) ;
 E_s           = NaN(Number_pixels,9131+1000) ;
 Cap_flux      = NaN(Number_pixels,9131+1000) ;
 tmp_ERA_save   = NaN(Number_pixels,9131)      ;

for i = 1:length(VIndex_rowcol)

    cur_in = VIndex_rowcol(i) ; 

%  i = spatial_extract    
srb_test = SRB_file.Clara_SRB_full_array_linear(cur_in,:)        ;
tmp_test = tmp_file.tmp_full_array(cur_in,:)                     ;
pre_test = pre_file.pre_full_array(cur_in,:)                     ;
spfh_test = spfh_file.spfh_full_array(cur_in,:)                  ; 
NDVI_test = NDVI_file.VIP_NDVI_array(cur_in,:)                   ;
tmp_ERA_test = tmp_ERA_file.T2M_ERA_daily_interp(cur_in,:)       ;           
tmp_ERA_test = tmp_ERA_test + 273.15                        ; 
tmp_ERA_save(counter,:) = tmp_ERA_test                            ; 

 if mean(NDVI_test,'omitnan') < 0  || range(NDVI_test) < 0.1 || all(isnan(NDVI_test))
     continue
 end


srb_test =      cat(2,srb_test(1:1000),srb_test) ; 
tmp_test =     cat(2,tmp_test(1:1000),tmp_test) ; 
tmp_ERA_test =  cat(2,tmp_ERA_test(1:1000),tmp_ERA_test) ; 
spfh_test =     cat(2,spfh_test(1:1000),spfh_test) ;  
pre_test =      cat(2,pre_test(1:1000),pre_test) ; 
NDVI_test =     cat(2,NDVI_test(1:1000),NDVI_test) ; 



 [Ts(counter,:),Tc(counter,:),ms(counter,:),mD(counter,:),Tra_s(counter,:),Tra_r(counter,:),E_s(counter,:),Cap_flux(counter,:), f_E(counter,:)] = ...
         The_Model_f_E_vec_sensitivity_study(srb_test,tmp_ERA_test ,spfh_test,pre_test./86400,NDVI_test,rootfrac,g_s_in_dummy,g_C_in_dummy,f_E_max,steps_per_day) ; 


counter = counter + 1 ; 

% counter
end


% function [Ts, TC, ms, mD, Tra_s_vec, Tra_r_vec,Es_vec, Cap_flux_vec, f_E] = The_Model_f_E_vec(F, T_R, q_R, P,NDVI,steps_per_day)
%                                % radiaiton Temp humidity precip              


% Ts = Ts(:,1001:end) ;
ms = ms(:,1001:end) ;
mD = mD(:,1001:end) ;
f_E = f_E(:,1001:end) ;


% rescale to same volumetric units this is critical and should rescale all
% SM values to roughly the same scale as SMAP. We can also it to 0-0.6?
% This would make it basically the same as SMAP might be useful
theta_max = 0.6 ; % pore space
ms = ms .* theta_max ;
mD = mD .* theta_max ;
tmp_ERA = tmp_ERA_test(:,1001:end) ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NDVI ANOMALY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('F:\CRUJRA\processed\datetime_1990_2014.mat')

[y_NDVI, m_NDVI, d_NDVI] = ymd(datetime_1990_2014) ;

f_E_NDVI_daily_anomaly = NaN(size(f_E)) ; 
tmp_ERA_daily_anomaly =  NaN(size(tmp_ERA)) ; 

for i = 1:size(f_E,1)

    f_E_dummy = f_E(i,:) ; 
    tmp_ERA_dummy = tmp_ERA_save(i,:) ; 

        VIP_NDVI_5y_movmean = movmean(f_E_dummy,1825) ; 
        tmp_ERA_5y_movmean =  movmean(tmp_ERA_dummy,1825) ; 


        % % remove from mean year
        NDVI_extract_day_mean = NaN(365,1) ; 
        tmp_extract_day_mean = NaN(365,1) ; 

        for j = 1:365
            [y, m, d] = ymd(datetime_1990_2014(j)) ; 
            extract_index = find(d == d_NDVI & m == m_NDVI) ; 

            NDVI_extract_day_mean(j) = mean(f_E_dummy(extract_index),'omitnan') ; 
            tmp_extract_day_mean(j) = mean(tmp_ERA_dummy(extract_index),'omitnan') ;             

        end
        NDVI_extract_day_mean = repmat(NDVI_extract_day_mean,[34,1]) ; 
        NDVI_extract_day_mean = NDVI_extract_day_mean(1:length(f_E_dummy)) ; 

        tmp_extract_day_mean = repmat(tmp_extract_day_mean,[34,1]) ; 
        tmp_extract_day_mean = tmp_extract_day_mean(1:length(tmp_ERA_dummy)) ;         


        f_E_NDVI_daily_anomaly(i,:) = f_E_dummy - NDVI_extract_day_mean' ; 
        tmp_ERA_daily_anomaly(i,:) =  tmp_ERA_dummy - tmp_extract_day_mean' ;         

% i
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRYDOWN DETECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('F:\CRUJRA\processed\datetime_1990_2014.mat')
sminterp = linspace(0,0.6,60) ; 

DD_dSM_dt_s_interpsm = NaN(400000,60) ; 
DD_tmp_ERA_anomaly_interpsm = NaN(400000,60) ; 
DD_NDVI_anomaly_interpsm = NaN(400000,60) ; 

DD_rows  =  NaN(400000,1) ; 
DD_cols  =  NaN(400000,1) ; 
DD_timev =  NaN(400000,2) ; 

cd('E:\Daten Baur\Matlab code\Project IGARSS multi frq tau\noodles_L_C_X')
DDLength = 4 ; 
NDry_count = 0 ; 
counter = 1 ;



for i = 1:size(ms,1)

time = 1:9131 ; 
ms_dummy = ms(i,:) ; 

 row_dummy = row_cru_v(i,:) ; 
 col_dummy = col_cru_v(i,:) ; 
 mD_dummy = mD(i,:) ;
 mS_dummy = ms(i,:) ;
 NDVI_anomaly_dummy = f_E_NDVI_daily_anomaly(i,:) ; 
 tmp_anomaly_dummy = tmp_ERA_daily_anomaly(i,:) ; 


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

     DD_tt = timev{j,1} ; 
     DD_dSM_dt = diff(DD_SM) ; 

     % only use strictly decreasing dds
     if any(DD_dSM_dt > 0)
         continue
     end

     DD_SM_short = DD_SM(1:end-1) ;  

     if numel(DD_SM_short)~=numel(unique(DD_SM_short))
         continue
     end

     DD_dSM_dt_s_interpsm(counter,:) = interp1(DD_SM_short,DD_dSM_dt,sminterp,'linear',NaN) ;  
     DD_NDVI_anomaly = NDVI_anomaly_dummy(DD_tt(1:end-1)) ; 
     DD_NDVI_anomaly_interpsm(counter,:) = interp1(DD_SM_short,DD_NDVI_anomaly,sminterp,'linear',NaN) ;  
     tmp_anomaly = tmp_anomaly_dummy(DD_tt(1:end-1)) ; 
     DD_tmp_ERA_anomaly_interpsm(counter,:) = interp1(DD_SM_short,tmp_anomaly,sminterp,'linear',NaN) ;  
     DD_rows(counter) = row_dummy ; 
     DD_cols(counter) = col_dummy ;
     DD_timev(counter,1) = DD_tt(1) ;
     DD_timev(counter,2) = DD_tt(end) ;      

     % if too short in SM space skip and not save
     if all(isnan(DD_dSM_dt_s_interpsm(counter,:)))
        continue
     end
     counter = counter +1 ; 
end
% i
end

DD_dSM_dt_s_interpsm(isnan(DD_rows),:) = [] ; 
DD_NDVI_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_tmp_ERA_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_cols(isnan(DD_rows),:) = [] ; 
DD_timev(isnan(DD_rows),:) = [] ; 
DD_rows(isnan(DD_rows),:) = [] ; 
% to mm/day
DD_dSM_dt_s_interpsm = DD_dSM_dt_s_interpsm.* -100 ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  fe binning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DD_dSM_dt_s_interpsm(isnan(DD_rows),:) = [] ; 
DD_NDVI_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_tmp_ERA_anomaly_interpsm(isnan(DD_rows),:) = [] ; 


    f_E_bins = linspace(0,1,21) ; 
    SM_bins  = linspace(0,0.6,21) ; 
    NDVI_anomaly_bins = linspace(-0.2,0.2,21) ; 
    T2M_anomaly_bins = linspace(-10,10,21) ;    


DD_dSM_dt_s_interpsm_mean = mean(DD_dSM_dt_s_interpsm,1,'omitnan') ; 

dSM_dt_NDVI_binned_sampling = NaN(20,20) ; 
dSM_dt_NDVI_binned = NaN(20,20) ; 


for sm_bin = 1:20

    sm_bin_s = SM_bins(sm_bin) ;
    sm_bin_e = SM_bins(sm_bin+1) ;
    index_sm = sminterp > sm_bin_s & sminterp < sm_bin_e ; 

    % assign data
    dSM_dt_s_anomaly_extract = DD_dSM_dt_s_interpsm(:,index_sm)   - DD_dSM_dt_s_interpsm_mean(index_sm)  ;  
    NDVI_anomaly_extract  = DD_NDVI_anomaly_interpsm(:,index_sm) ;  
    tmp_ERA_anomaly_extract  = DD_tmp_ERA_anomaly_interpsm(:,index_sm) ;  


    for f_E_bin = 1:20

            NDVI_anomaly_bin_s = NDVI_anomaly_bins(f_E_bin) ;
            NDVI_anomaly_bin_e = NDVI_anomaly_bins(f_E_bin+1) ;

            extract_index = NDVI_anomaly_extract > NDVI_anomaly_bin_s & NDVI_anomaly_extract < NDVI_anomaly_bin_e ;         

            dSM_dt_extract = dSM_dt_s_anomaly_extract(extract_index) ; 

            dSM_dt_NDVI_binned_sampling(f_E_bin,sm_bin) = sum(~isnan(dSM_dt_extract)) ; 
            dSM_dt_NDVI_binned(f_E_bin,sm_bin) = median(dSM_dt_extract); 


    end

% sm_bin
end


[FX_NDVI, FY_NDVI] =  gradient(wiener2(dSM_dt_NDVI_binned)) ; 

[rowd,cold] = ind2sub([30 30],k) ; 
dSM_dt_gradient_sensitivity_gs_gc(rowd,cold) = median(FY_NDVI(:),'omitnan') ; 



k


end



% g_s_in_array(1:10)
% g_C_in_array(1:10)
dSM_dt_gradient_sensitivity_gs_gc_US = dSM_dt_gradient_sensitivity_gs_gc ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_gradient_sensitivity_gs_gc_US','dSM_dt_gradient_sensitivity_gs_gc_US')


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(g_s_in_array2D,g_C_in_array2D,dSM_dt_gradient_sensitivity_gs_gc) ;
% h1 = imagesc(g_s_in_array2D(1,:),g_C_in_array2D(:,1),dSM_dt_gradient_sensitivity_gs_gc) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.1 0.1]) 
colormap(redblue_color) 
xlim([0 0.0101])
ylim([0 0.0101])
cbr = colorbar ;
cbr.Label.String = "\Delta SM loss / \Delta NDVI [mm/day]";
xlabel('soil conductance')
ylabel('canopy conductance')
pbaspect([1 1 1])
fontsize(16,'points')

axes_diff_Position = get(gca, 'Position');
% 3.627863592950994,2.307960416666667,21.62764834259246,17.09988854166667

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[28,    2+(17.1)/2+1,     0,     (17.1)/2-1  ]) ; 
 set(arrow2,'Position',[28,    2+(17.1)/2-1,     0,    -(17.1)/2+1 ]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss (Extraction)' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[28+1,    2+(17.1)/2+7.5,   10,     1  ]) ; 
set(gca,'Box','on');
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss (Cover)' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[28+1,     2+(17.1)/2-1,     10,  1 ]) ; 
set(gca,'Box','on');
% title(num2str([xs, ys]))
title('US')

saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures_revisions_01\dSM_dt_gradient_gsgc_2D_US','svg')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures_revisions_01\dSM_dt_gradient_gsgc_2D_US','png')
close






imagesc(dSM_dt_gradient_sensitivity_gs_gc)
colormap(redblue_color)
clim([-0.1 0.1])

cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 

imagesc(dSM_dt_NDVI_binned)
colormap(redblue_color)
clim([-1 1])

imagesc(wiener2(dSM_dt_NDVI_binned))
colormap(redblue_color)
clim([-1 1])







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup up to test sensitvity of f_e sclaing. Should be straightforward shift of cover or extraction effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% test = ncinfo("crujra.v2.4.5d.spfh.1990.365d.noc.nc") ; 
testinfo = ncinfo("F:\CRUJRA\crujra.v2.4.5d.tmp.2022.365d.noc.nc") ; 

cd('F:\CRUJRA\processed')

load('crujra_datetime_vector.mat')
load('datetime_1990_2014.mat')
load('row_cru_v.mat')
load('col_cru_v.mat')
load('lat_cru_v.mat')
load('lon_cru_v.mat')

SRB_file = matfile("Clara_SRB_full_array_linear.mat") ; 
tmp_file = matfile("tmp_full_array.mat") ;
pre_file = matfile("pre_full_array.mat") ;
spfh_file = matfile("spfh_full_array.mat") ;

% get VIP monthly NDVI for f_E
NDVI_file = matfile('F:\VIP_NDVI\VIP_NDVI_array_linear.mat') ; 


%tmp from ERA 5 monthly
cd('F:\ERA_5_monthly')
tmp_ERA_file = matfile("T2M_ERA_daily_interp.mat") ;

% build 2D for pixel selection
load('F:\CRUJRA\Ancil_data\CRUJRA_land_mask.mat') ;
load('F:\CRUJRA\Ancil_data\CRUJRA_land_mask_vec.mat') ;

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d')
load('DD_dSM_dt_NDVI_gradient_zeppe_ERA.mat')


figure
imagesc(CRUJRA_land_mask)
[xs, ys] = getpts() ; xs = round(xs) ; ys = round(ys) ; 
close

% pixels used so far
% ys = [52 57] ; xs = [590 596] ; % boreal
% ys = [210 215] ; xs = [264 271] ; % Brazil
% ys = [139 145] ; xs = [512 520] ; % India
% ys = [156 160] ; xs = [388 395] ; % Sahel
% ys = [88 100] ; xs = [115 124] ; % US


% MODEL

Number_pixels = length(xs) ; 
Number_pixels_square = (xs(2)-xs(1)) .* (ys(2)-ys(1)) ; 


cd('F:\Matlab_Code\Zeppe_model')
% shorter with 1000 steps spinup
 Ts            = NaN(Number_pixels,9131+1000) ; 
 Tc            = NaN(Number_pixels,9131+1000) ;
 Tra_s         = NaN(Number_pixels,9131+1000) ;
 Tra_r         = NaN(Number_pixels,9131+1000) ;
 ms            = NaN(Number_pixels,9131+1000) ;
 mD            = NaN(Number_pixels,9131+1000) ;
 f_E           = NaN(Number_pixels,9131+1000) ;
 E_s           = NaN(Number_pixels,9131+1000) ;
 Cap_flux      = NaN(Number_pixels,9131+1000) ;

tmp_ERA_save   = NaN(Number_pixels,9131)      ;

steps_per_day = 1 ;

cd('F:\Matlab_Code\Zeppe_model')


% find i based on location
for i = 1:length(xs) 
[VIndex_rowcol(i)] = find(col_cru_v == xs(i) & row_cru_v == ys(i)) ; 
end


% find area of square
VIndex_rowcol = find((col_cru_v >= xs(1) & col_cru_v <= xs(2)) & (row_cru_v >= ys(1) & ...
    row_cru_v <= ys(2))) ; 
% VIndex_rowcol = sort(VIndex_rowcol) ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT array 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        f_E = (NDVI - (0.0)) ./ (0.8 - (0.0)) ; 
f_E_scaling_v =linspace(0.1,1,20) ; % more than enough, visualize gradient effect as dots.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dSM_dt_gradient_sensitivity_fe_max = NaN(size(f_E_scaling_v)) ; 



for k = 1:length(f_E_scaling_v)

cd('F:\Matlab_Code\Zeppe_model')
rootfrac = 0.7 ; 
% keep default values
g_s_in_dummy = 1/2000;  
g_C_in_dummy = 1/800;
f_E_max = f_E_scaling_v(k) ; 
counter = 1 ; 


% i = VIndex_rowcol(1) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Number_pixels = size(VIndex_rowcol,1) ;
Ts            = NaN(Number_pixels,9131+1000) ; 
 Tc            = NaN(Number_pixels,9131+1000) ;
 Tra_s         = NaN(Number_pixels,9131+1000) ;
 Tra_r         = NaN(Number_pixels,9131+1000) ;
 ms            = NaN(Number_pixels,9131+1000) ;
 mD            = NaN(Number_pixels,9131+1000) ;
 f_E           = NaN(Number_pixels,9131+1000) ;
 E_s           = NaN(Number_pixels,9131+1000) ;
 Cap_flux      = NaN(Number_pixels,9131+1000) ;
 tmp_ERA_save   = NaN(Number_pixels,9131)      ;

for i = 1:length(VIndex_rowcol)

    cur_in = VIndex_rowcol(i) ; 

%  i = spatial_extract    
srb_test = SRB_file.Clara_SRB_full_array_linear(cur_in,:)        ;
tmp_test = tmp_file.tmp_full_array(cur_in,:)                     ;
pre_test = pre_file.pre_full_array(cur_in,:)                     ;
spfh_test = spfh_file.spfh_full_array(cur_in,:)                  ; 
NDVI_test = NDVI_file.VIP_NDVI_array(cur_in,:)                   ;
tmp_ERA_test = tmp_ERA_file.T2M_ERA_daily_interp(cur_in,:)       ;           
tmp_ERA_test = tmp_ERA_test + 273.15                        ; 
tmp_ERA_save(counter,:) = tmp_ERA_test                            ; 

 if mean(NDVI_test,'omitnan') < 0  || range(NDVI_test) < 0.1 || all(isnan(NDVI_test))
     continue
 end


srb_test =      cat(2,srb_test(1:1000),srb_test) ; 
tmp_test =     cat(2,tmp_test(1:1000),tmp_test) ; 
tmp_ERA_test =  cat(2,tmp_ERA_test(1:1000),tmp_ERA_test) ; 
spfh_test =     cat(2,spfh_test(1:1000),spfh_test) ;  
pre_test =      cat(2,pre_test(1:1000),pre_test) ; 
NDVI_test =     cat(2,NDVI_test(1:1000),NDVI_test) ; 



 [Ts(counter,:),Tc(counter,:),ms(counter,:),mD(counter,:),Tra_s(counter,:),Tra_r(counter,:),E_s(counter,:),Cap_flux(counter,:), f_E(counter,:)] = ...
         The_Model_f_E_vec_sensitivity_study(srb_test,tmp_ERA_test ,spfh_test,pre_test./86400,NDVI_test,rootfrac,g_s_in_dummy,g_C_in_dummy,f_E_max,steps_per_day) ; 


counter = counter + 1 ; 

% counter
end


% function [Ts, TC, ms, mD, Tra_s_vec, Tra_r_vec,Es_vec, Cap_flux_vec, f_E] = The_Model_f_E_vec(F, T_R, q_R, P,NDVI,steps_per_day)
%                                % radiaiton Temp humidity precip              


% Ts = Ts(:,1001:end) ;
ms = ms(:,1001:end) ;
mD = mD(:,1001:end) ;
f_E = f_E(:,1001:end) ;


% rescale to same volumetric units this is critical and should rescale all
% SM values to roughly the same scale as SMAP. We can also it to 0-0.6?
% This would make it basically the same as SMAP might be useful
theta_max = 0.6 ; % pore space
ms = ms .* theta_max ;
mD = mD .* theta_max ;
tmp_ERA = tmp_ERA_test(:,1001:end) ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NDVI ANOMALY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('F:\CRUJRA\processed\datetime_1990_2014.mat')

[y_NDVI, m_NDVI, d_NDVI] = ymd(datetime_1990_2014) ;

f_E_NDVI_daily_anomaly = NaN(size(f_E)) ; 
tmp_ERA_daily_anomaly =  NaN(size(tmp_ERA)) ; 

for i = 1:size(f_E,1)

    f_E_dummy = f_E(i,:) ; 
    tmp_ERA_dummy = tmp_ERA_save(i,:) ; 

        VIP_NDVI_5y_movmean = movmean(f_E_dummy,1825) ; 
        tmp_ERA_5y_movmean =  movmean(tmp_ERA_dummy,1825) ; 


        % % remove from mean year
        NDVI_extract_day_mean = NaN(365,1) ; 
        tmp_extract_day_mean = NaN(365,1) ; 

        for j = 1:365
            [y, m, d] = ymd(datetime_1990_2014(j)) ; 
            extract_index = find(d == d_NDVI & m == m_NDVI) ; 

            NDVI_extract_day_mean(j) = mean(f_E_dummy(extract_index),'omitnan') ; 
            tmp_extract_day_mean(j) = mean(tmp_ERA_dummy(extract_index),'omitnan') ;             

        end
        NDVI_extract_day_mean = repmat(NDVI_extract_day_mean,[34,1]) ; 
        NDVI_extract_day_mean = NDVI_extract_day_mean(1:length(f_E_dummy)) ; 

        tmp_extract_day_mean = repmat(tmp_extract_day_mean,[34,1]) ; 
        tmp_extract_day_mean = tmp_extract_day_mean(1:length(tmp_ERA_dummy)) ;         


        f_E_NDVI_daily_anomaly(i,:) = f_E_dummy - NDVI_extract_day_mean' ; 
        tmp_ERA_daily_anomaly(i,:) =  tmp_ERA_dummy - tmp_extract_day_mean' ;         

% i
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRYDOWN DETECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('F:\CRUJRA\processed\datetime_1990_2014.mat')
sminterp = linspace(0,0.6,60) ; 

DD_dSM_dt_s_interpsm = NaN(400000,60) ; 
DD_tmp_ERA_anomaly_interpsm = NaN(400000,60) ; 
DD_NDVI_anomaly_interpsm = NaN(400000,60) ; 

DD_rows  =  NaN(400000,1) ; 
DD_cols  =  NaN(400000,1) ; 
DD_timev =  NaN(400000,2) ; 

cd('E:\Daten Baur\Matlab code\Project IGARSS multi frq tau\noodles_L_C_X')
DDLength = 4 ; 
NDry_count = 0 ; 
counter = 1 ;



for i = 1:size(ms,1)

time = 1:9131 ; 
ms_dummy = ms(i,:) ; 

 row_dummy = row_cru_v(i,:) ; 
 col_dummy = col_cru_v(i,:) ; 
 mD_dummy = mD(i,:) ;
 mS_dummy = ms(i,:) ;
 NDVI_anomaly_dummy = f_E_NDVI_daily_anomaly(i,:) ; 
 tmp_anomaly_dummy = tmp_ERA_daily_anomaly(i,:) ; 


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

     DD_tt = timev{j,1} ; 
     DD_dSM_dt = diff(DD_SM) ; 

     % only use strictly decreasing dds
     if any(DD_dSM_dt > 0)
         continue
     end

     DD_SM_short = DD_SM(1:end-1) ;  

     if numel(DD_SM_short)~=numel(unique(DD_SM_short))
         continue
     end

     DD_dSM_dt_s_interpsm(counter,:) = interp1(DD_SM_short,DD_dSM_dt,sminterp,'linear',NaN) ;  
     DD_NDVI_anomaly = NDVI_anomaly_dummy(DD_tt(1:end-1)) ; 
     DD_NDVI_anomaly_interpsm(counter,:) = interp1(DD_SM_short,DD_NDVI_anomaly,sminterp,'linear',NaN) ;  
     tmp_anomaly = tmp_anomaly_dummy(DD_tt(1:end-1)) ; 
     DD_tmp_ERA_anomaly_interpsm(counter,:) = interp1(DD_SM_short,tmp_anomaly,sminterp,'linear',NaN) ;  
     DD_rows(counter) = row_dummy ; 
     DD_cols(counter) = col_dummy ;
     DD_timev(counter,1) = DD_tt(1) ;
     DD_timev(counter,2) = DD_tt(end) ;      

     % if too short in SM space skip and not save
     if all(isnan(DD_dSM_dt_s_interpsm(counter,:)))
        continue
     end
     counter = counter +1 ; 
end
% i
end

DD_dSM_dt_s_interpsm(isnan(DD_rows),:) = [] ; 
DD_NDVI_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_tmp_ERA_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_cols(isnan(DD_rows),:) = [] ; 
DD_timev(isnan(DD_rows),:) = [] ; 
DD_rows(isnan(DD_rows),:) = [] ; 
% to mm/day
DD_dSM_dt_s_interpsm = DD_dSM_dt_s_interpsm.* -100 ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  fe binning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DD_dSM_dt_s_interpsm(isnan(DD_rows),:) = [] ; 
DD_NDVI_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_tmp_ERA_anomaly_interpsm(isnan(DD_rows),:) = [] ; 


    f_E_bins = linspace(0,1,21) ; 
    SM_bins  = linspace(0,0.6,21) ; 
    NDVI_anomaly_bins = linspace(-0.2,0.2,21) ; 
    T2M_anomaly_bins = linspace(-10,10,21) ;    


DD_dSM_dt_s_interpsm_mean = mean(DD_dSM_dt_s_interpsm,1,'omitnan') ; 

dSM_dt_NDVI_binned_sampling = NaN(20,20) ; 
dSM_dt_NDVI_binned = NaN(20,20) ; 


for sm_bin = 1:20

    sm_bin_s = SM_bins(sm_bin) ;
    sm_bin_e = SM_bins(sm_bin+1) ;
    index_sm = sminterp > sm_bin_s & sminterp < sm_bin_e ; 

    % assign data
    dSM_dt_s_anomaly_extract = DD_dSM_dt_s_interpsm(:,index_sm)   - DD_dSM_dt_s_interpsm_mean(index_sm)  ;  
    NDVI_anomaly_extract  = DD_NDVI_anomaly_interpsm(:,index_sm) ;  
    tmp_ERA_anomaly_extract  = DD_tmp_ERA_anomaly_interpsm(:,index_sm) ;  


    for f_E_bin = 1:20

            NDVI_anomaly_bin_s = NDVI_anomaly_bins(f_E_bin) ;
            NDVI_anomaly_bin_e = NDVI_anomaly_bins(f_E_bin+1) ;

            extract_index = NDVI_anomaly_extract > NDVI_anomaly_bin_s & NDVI_anomaly_extract < NDVI_anomaly_bin_e ;         

            dSM_dt_extract = dSM_dt_s_anomaly_extract(extract_index) ; 

            dSM_dt_NDVI_binned_sampling(f_E_bin,sm_bin) = sum(~isnan(dSM_dt_extract)) ; 
            dSM_dt_NDVI_binned(f_E_bin,sm_bin) = median(dSM_dt_extract); 


    end

% sm_bin
end

% dSM_dt_NDVI_binned(dSM_dt_NDVI_binned_sampling < 10) ; 
[FX_NDVI, FY_NDVI] =  gradient(wiener2(dSM_dt_NDVI_binned)) ; 
% 
% figure
% imagesc(dSM_dt_NDVI_binned)
% colormap(redblue_color)
% clim([-0.3 0.3])
% title(num2str(k))

% imagesc(dSM_dt_NDVI_binned_sampling)


[rowd,cold] = ind2sub([30 30],k) ; 
dSM_dt_gradient_sensitivity_fe_max(k) = median(FY_NDVI(:),'omitnan') ; 



k

end




cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
plot(f_E_scaling_v,dSM_dt_gradient_sensitivity_fe_max,'o-','LineWidth',2,'Color',redblue_color(20,:))
yline(0,'--')
ylim([-0.08 0.08])
xlabel('NDVI max for scaling f_E')
ylabel('\Delta SM loss / \Delta NDVI [mm/day]')
pbaspect([1 1 1])

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[26,    2+(17.1)/2+1,     0,     (17.1)/2-1  ]) ; 
 set(arrow2,'Position',[26,    2+(17.1)/2-1,     0,    -(17.1)/2+1 ]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss (Extraction)' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[26+1,    2+(17.1)/2+7.5,   10,     1  ]) ; 
set(gca,'Box','on');
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss (Cover)' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[26+1,     2+(17.1)/2-1,     10,  1 ]) ; 
set(gca,'Box','on');
% title(num2str([xs, ys]))
 title('US')
fontsize(16,'points')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures_revisions_01\fe_scaling_dSM_dt_impact_US','svg')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures_revisions_01\fe_scaling_dSM_dt_impact_US','png')
close




%% Extra section make map with AOIS in rectangles
NDVI_3D_file = matfile('F:\VIP_NDVI\VIP_NDVI_monthly.mat') ; 
NDVI_mean = mean(NDVI_3D_file.VIP_NDVI_monthly,3,'omitnan') ; 

Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
imagesc(linspace(1,720,1440),linspace(1,360,720),NDVI_mean)
% add AOI rectangles
hold on
rectangle('Position',[590, 52, 596-590, 57-52],'EdgeColor','r','LineWidth',1.5) ; % boreal
rectangle('Position',[264, 210, 271-264, 215-210],'EdgeColor','r','LineWidth',1.5) ; % Brazil
rectangle('Position',[512, 139, 520-512, 145-139],'EdgeColor','r','LineWidth',1.5) ;  % India
rectangle('Position',[388, 156, 395-388, 160-156],'EdgeColor','r','LineWidth',1.5) ;  % Sahel
rectangle('Position',[115, 88, 124-115, 100-88],'EdgeColor','r','LineWidth',1.5) ; % US
cbr = colorbar ; 
cbr.Label.String = "Average NDVI [-]";
clim([0 1])
ylabel('rows relative to -90:90 degrees latitude')
xlabel('cols relative to -180:180 degrees longitude')
pbaspect([360 180 1])
fontsize(16,'points')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures_revisions_01\AOI_map_NDVI','svg')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures_revisions_01\AOI_map_NDVI','png')
close

cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Complete analysis with sensitivt of roothing fraction 0-1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear

% test = ncinfo("crujra.v2.4.5d.spfh.1990.365d.noc.nc") ; 
testinfo = ncinfo("F:\CRUJRA\crujra.v2.4.5d.tmp.2022.365d.noc.nc") ; 

cd('F:\CRUJRA\processed')

load('crujra_datetime_vector.mat')
load('datetime_1990_2014.mat')
load('row_cru_v.mat')
load('col_cru_v.mat')
load('lat_cru_v.mat')
load('lon_cru_v.mat')

SRB_file = matfile("Clara_SRB_full_array_linear.mat") ; 
tmp_file = matfile("tmp_full_array.mat") ;
pre_file = matfile("pre_full_array.mat") ;
spfh_file = matfile("spfh_full_array.mat") ;

% get VIP monthly NDVI for f_E
NDVI_file = matfile('F:\VIP_NDVI\VIP_NDVI_array_linear.mat') ; 


%tmp from ERA 5 monthly
cd('F:\ERA_5_monthly')
tmp_ERA_file = matfile("T2M_ERA_daily_interp.mat") ;

% build 2D for pixel selection
load('F:\CRUJRA\Ancil_data\CRUJRA_land_mask.mat') ;
load('F:\CRUJRA\Ancil_data\CRUJRA_land_mask_vec.mat') ;

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d')
load('DD_dSM_dt_NDVI_gradient_zeppe_ERA.mat')


figure
imagesc(CRUJRA_land_mask)
[xs, ys] = getpts() ; xs = round(xs) ; ys = round(ys) ; 
close

% pixels used so far
% ys = [52 57] ; xs = [590 596] ; % boreal
% ys = [210 215] ; xs = [264 271] ; % Brazil
% ys = [139 145] ; xs = [512 520] ; % India
% ys = [156 160] ; xs = [388 395] ; % Sahel
% ys = [88 100] ; xs = [115 124] ; % US


% MODEL

Number_pixels = length(xs) ; 
Number_pixels_square = (xs(2)-xs(1)) .* (ys(2)-ys(1)) ; 


cd('F:\Matlab_Code\Zeppe_model')
% shorter with 1000 steps spinup
 Ts            = NaN(Number_pixels,9131+1000) ; 
 Tc            = NaN(Number_pixels,9131+1000) ;
 Tra_s         = NaN(Number_pixels,9131+1000) ;
 Tra_r         = NaN(Number_pixels,9131+1000) ;
 ms            = NaN(Number_pixels,9131+1000) ;
 mD            = NaN(Number_pixels,9131+1000) ;
 f_E           = NaN(Number_pixels,9131+1000) ;
 E_s           = NaN(Number_pixels,9131+1000) ;
 Cap_flux      = NaN(Number_pixels,9131+1000) ;

tmp_ERA_save   = NaN(Number_pixels,9131)      ;

steps_per_day = 1 ;

cd('F:\Matlab_Code\Zeppe_model')


% find i based on location
for i = 1:length(xs) 
[VIndex_rowcol(i)] = find(col_cru_v == xs(i) & row_cru_v == ys(i)) ; 
end


% find area of square
VIndex_rowcol = find((col_cru_v >= xs(1) & col_cru_v <= xs(2)) & (row_cru_v >= ys(1) & ...
    row_cru_v <= ys(2))) ; 
% VIndex_rowcol = sort(VIndex_rowcol) ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT array 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        f_E = (NDVI - (0.0)) ./ (0.8 - (0.0)) ; 
root_fraction_v =linspace(0,1,20) ; % more than enough, visualize gradient effect as dots.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dSM_dt_gradient_sensitivity_rootfrac = NaN(size(root_fraction_v)) ; 



for k = 1:length(root_fraction_v)

cd('F:\Matlab_Code\Zeppe_model')
% rootfrac = 0.7 ; 
% keep default values
g_s_in_dummy = 1/2000;  
g_C_in_dummy = 1/800;
f_E_max = 0.8 ; 
root_frac_in = root_fraction_v(k) ; 
counter = 1 ; 


% i = VIndex_rowcol(1) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Number_pixels = size(VIndex_rowcol,1) ;
Ts            = NaN(Number_pixels,9131+1000) ; 
 Tc            = NaN(Number_pixels,9131+1000) ;
 Tra_s         = NaN(Number_pixels,9131+1000) ;
 Tra_r         = NaN(Number_pixels,9131+1000) ;
 ms            = NaN(Number_pixels,9131+1000) ;
 mD            = NaN(Number_pixels,9131+1000) ;
 f_E           = NaN(Number_pixels,9131+1000) ;
 E_s           = NaN(Number_pixels,9131+1000) ;
 Cap_flux      = NaN(Number_pixels,9131+1000) ;
 tmp_ERA_save   = NaN(Number_pixels,9131)      ;

for i = 1:length(VIndex_rowcol)

    cur_in = VIndex_rowcol(i) ; 

%  i = spatial_extract    
srb_test = SRB_file.Clara_SRB_full_array_linear(cur_in,:)        ;
tmp_test = tmp_file.tmp_full_array(cur_in,:)                     ;
pre_test = pre_file.pre_full_array(cur_in,:)                     ;
spfh_test = spfh_file.spfh_full_array(cur_in,:)                  ; 
NDVI_test = NDVI_file.VIP_NDVI_array(cur_in,:)                   ;
tmp_ERA_test = tmp_ERA_file.T2M_ERA_daily_interp(cur_in,:)       ;           
tmp_ERA_test = tmp_ERA_test + 273.15                        ; 
tmp_ERA_save(counter,:) = tmp_ERA_test                            ; 

 if mean(NDVI_test,'omitnan') < 0  || range(NDVI_test) < 0.1 || all(isnan(NDVI_test))
     continue
 end


srb_test =      cat(2,srb_test(1:1000),srb_test) ; 
tmp_test =     cat(2,tmp_test(1:1000),tmp_test) ; 
tmp_ERA_test =  cat(2,tmp_ERA_test(1:1000),tmp_ERA_test) ; 
spfh_test =     cat(2,spfh_test(1:1000),spfh_test) ;  
pre_test =      cat(2,pre_test(1:1000),pre_test) ; 
NDVI_test =     cat(2,NDVI_test(1:1000),NDVI_test) ; 



 [Ts(counter,:),Tc(counter,:),ms(counter,:),mD(counter,:),Tra_s(counter,:),Tra_r(counter,:),E_s(counter,:),Cap_flux(counter,:), f_E(counter,:)] = ...
         The_Model_f_E_vec_sensitivity_study(srb_test,tmp_ERA_test ,spfh_test,pre_test./86400,NDVI_test,root_frac_in,g_s_in_dummy,g_C_in_dummy,f_E_max,steps_per_day) ; 


counter = counter + 1 ; 

% counter
end


% function [Ts, TC, ms, mD, Tra_s_vec, Tra_r_vec,Es_vec, Cap_flux_vec, f_E] = The_Model_f_E_vec(F, T_R, q_R, P,NDVI,steps_per_day)
%                                % radiaiton Temp humidity precip              


% Ts = Ts(:,1001:end) ;
ms = ms(:,1001:end) ;
mD = mD(:,1001:end) ;
f_E = f_E(:,1001:end) ;


% rescale to same volumetric units this is critical and should rescale all
% SM values to roughly the same scale as SMAP. We can also it to 0-0.6?
% This would make it basically the same as SMAP might be useful
theta_max = 0.6 ; % pore space
ms = ms .* theta_max ;
mD = mD .* theta_max ;
tmp_ERA = tmp_ERA_test(:,1001:end) ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NDVI ANOMALY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('F:\CRUJRA\processed\datetime_1990_2014.mat')

[y_NDVI, m_NDVI, d_NDVI] = ymd(datetime_1990_2014) ;

f_E_NDVI_daily_anomaly = NaN(size(f_E)) ; 
tmp_ERA_daily_anomaly =  NaN(size(tmp_ERA)) ; 

for i = 1:size(f_E,1)

    f_E_dummy = f_E(i,:) ; 
    tmp_ERA_dummy = tmp_ERA_save(i,:) ; 

        VIP_NDVI_5y_movmean = movmean(f_E_dummy,1825) ; 
        tmp_ERA_5y_movmean =  movmean(tmp_ERA_dummy,1825) ; 


        % % remove from mean year
        NDVI_extract_day_mean = NaN(365,1) ; 
        tmp_extract_day_mean = NaN(365,1) ; 

        for j = 1:365
            [y, m, d] = ymd(datetime_1990_2014(j)) ; 
            extract_index = find(d == d_NDVI & m == m_NDVI) ; 

            NDVI_extract_day_mean(j) = mean(f_E_dummy(extract_index),'omitnan') ; 
            tmp_extract_day_mean(j) = mean(tmp_ERA_dummy(extract_index),'omitnan') ;             

        end
        NDVI_extract_day_mean = repmat(NDVI_extract_day_mean,[34,1]) ; 
        NDVI_extract_day_mean = NDVI_extract_day_mean(1:length(f_E_dummy)) ; 

        tmp_extract_day_mean = repmat(tmp_extract_day_mean,[34,1]) ; 
        tmp_extract_day_mean = tmp_extract_day_mean(1:length(tmp_ERA_dummy)) ;         


        f_E_NDVI_daily_anomaly(i,:) = f_E_dummy - NDVI_extract_day_mean' ; 
        tmp_ERA_daily_anomaly(i,:) =  tmp_ERA_dummy - tmp_extract_day_mean' ;         

% i
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRYDOWN DETECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('F:\CRUJRA\processed\datetime_1990_2014.mat')
sminterp = linspace(0,0.6,60) ; 

DD_dSM_dt_s_interpsm = NaN(400000,60) ; 
DD_tmp_ERA_anomaly_interpsm = NaN(400000,60) ; 
DD_NDVI_anomaly_interpsm = NaN(400000,60) ; 

DD_rows  =  NaN(400000,1) ; 
DD_cols  =  NaN(400000,1) ; 
DD_timev =  NaN(400000,2) ; 

cd('E:\Daten Baur\Matlab code\Project IGARSS multi frq tau\noodles_L_C_X')
DDLength = 4 ; 
NDry_count = 0 ; 
counter = 1 ;



for i = 1:size(ms,1)

time = 1:9131 ; 
ms_dummy = ms(i,:) ; 

 row_dummy = row_cru_v(i,:) ; 
 col_dummy = col_cru_v(i,:) ; 
 mD_dummy = mD(i,:) ;
 mS_dummy = ms(i,:) ;
 NDVI_anomaly_dummy = f_E_NDVI_daily_anomaly(i,:) ; 
 tmp_anomaly_dummy = tmp_ERA_daily_anomaly(i,:) ; 


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

     DD_tt = timev{j,1} ; 
     DD_dSM_dt = diff(DD_SM) ; 

     % only use strictly decreasing dds
     if any(DD_dSM_dt > 0)
         continue
     end

     DD_SM_short = DD_SM(1:end-1) ;  

     if numel(DD_SM_short)~=numel(unique(DD_SM_short))
         continue
     end

     DD_dSM_dt_s_interpsm(counter,:) = interp1(DD_SM_short,DD_dSM_dt,sminterp,'linear',NaN) ;  
     DD_NDVI_anomaly = NDVI_anomaly_dummy(DD_tt(1:end-1)) ; 
     DD_NDVI_anomaly_interpsm(counter,:) = interp1(DD_SM_short,DD_NDVI_anomaly,sminterp,'linear',NaN) ;  
     tmp_anomaly = tmp_anomaly_dummy(DD_tt(1:end-1)) ; 
     DD_tmp_ERA_anomaly_interpsm(counter,:) = interp1(DD_SM_short,tmp_anomaly,sminterp,'linear',NaN) ;  
     DD_rows(counter) = row_dummy ; 
     DD_cols(counter) = col_dummy ;
     DD_timev(counter,1) = DD_tt(1) ;
     DD_timev(counter,2) = DD_tt(end) ;      

     % if too short in SM space skip and not save
     if all(isnan(DD_dSM_dt_s_interpsm(counter,:)))
        continue
     end
     counter = counter +1 ; 
end
% i
end

DD_dSM_dt_s_interpsm(isnan(DD_rows),:) = [] ; 
DD_NDVI_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_tmp_ERA_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_cols(isnan(DD_rows),:) = [] ; 
DD_timev(isnan(DD_rows),:) = [] ; 
DD_rows(isnan(DD_rows),:) = [] ; 
% to mm/day
DD_dSM_dt_s_interpsm = DD_dSM_dt_s_interpsm.* -100 ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  fe binning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DD_dSM_dt_s_interpsm(isnan(DD_rows),:) = [] ; 
DD_NDVI_anomaly_interpsm(isnan(DD_rows),:) = [] ; 
DD_tmp_ERA_anomaly_interpsm(isnan(DD_rows),:) = [] ; 


    f_E_bins = linspace(0,1,21) ; 
    SM_bins  = linspace(0,0.6,21) ; 
    NDVI_anomaly_bins = linspace(-0.2,0.2,21) ; 
    T2M_anomaly_bins = linspace(-10,10,21) ;    


DD_dSM_dt_s_interpsm_mean = mean(DD_dSM_dt_s_interpsm,1,'omitnan') ; 

dSM_dt_NDVI_binned_sampling = NaN(20,20) ; 
dSM_dt_NDVI_binned = NaN(20,20) ; 


for sm_bin = 1:20

    sm_bin_s = SM_bins(sm_bin) ;
    sm_bin_e = SM_bins(sm_bin+1) ;
    index_sm = sminterp > sm_bin_s & sminterp < sm_bin_e ; 

    % assign data
    dSM_dt_s_anomaly_extract = DD_dSM_dt_s_interpsm(:,index_sm)   - DD_dSM_dt_s_interpsm_mean(index_sm)  ;  
    NDVI_anomaly_extract  = DD_NDVI_anomaly_interpsm(:,index_sm) ;  
    tmp_ERA_anomaly_extract  = DD_tmp_ERA_anomaly_interpsm(:,index_sm) ;  


    for f_E_bin = 1:20

            NDVI_anomaly_bin_s = NDVI_anomaly_bins(f_E_bin) ;
            NDVI_anomaly_bin_e = NDVI_anomaly_bins(f_E_bin+1) ;

            extract_index = NDVI_anomaly_extract > NDVI_anomaly_bin_s & NDVI_anomaly_extract < NDVI_anomaly_bin_e ;         

            dSM_dt_extract = dSM_dt_s_anomaly_extract(extract_index) ; 

            dSM_dt_NDVI_binned_sampling(f_E_bin,sm_bin) = sum(~isnan(dSM_dt_extract)) ; 
            dSM_dt_NDVI_binned(f_E_bin,sm_bin) = median(dSM_dt_extract); 


    end

% sm_bin
end

% dSM_dt_NDVI_binned(dSM_dt_NDVI_binned_sampling < 10) ; 
[FX_NDVI, FY_NDVI] =  gradient(wiener2(dSM_dt_NDVI_binned)) ; 
% 
% figure
% imagesc(dSM_dt_NDVI_binned)
% colormap(redblue_color)
% clim([-0.3 0.3])
% title(num2str(k))

% imagesc(dSM_dt_NDVI_binned_sampling)


[rowd,cold] = ind2sub([30 30],k) ; 
dSM_dt_gradient_sensitivity_rootfrac(k) = median(FY_NDVI(:),'omitnan') ; 



k

end




cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
plot(root_fraction_v,dSM_dt_gradient_sensitivity_rootfrac,'o-','LineWidth',2,'Color',redblue_color(20,:))
yline(0,'--')
ylim([-0.08 0.08])
xlabel('root fraction top layer')
ylabel('\Delta SM loss / \Delta NDVI [mm/day]')
pbaspect([1 1 1])

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[26,    2+(17.1)/2+1,     0,     (17.1)/2-1  ]) ; 
 set(arrow2,'Position',[26,    2+(17.1)/2-1,     0,    -(17.1)/2+1 ]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss (Extraction)' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[26+1,    2+(17.1)/2+7.5,   10,     1  ]) ; 
set(gca,'Box','on');
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss (Cover)' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[26+1,     2+(17.1)/2-1,     10,  1 ]) ; 
set(gca,'Box','on');
% title(num2str([xs, ys]))
title('Boreal')
fontsize(16,'points')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures_revisions_01\rootfrac_dSM_dt_impact_boreal','svg')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures_revisions_01\rootfrac_dSM_dt_impact_boreal','png')
close





