%% MJB 06-05-2024 Do Zeppe Model plots but for 2.5 degree bins 

clear

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05')
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_06_mm_d')
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_07_mm_d_low_SM')
cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08')
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests')



load('DD_dSM_dt_shallow_interpsm.mat')
load('DD_NDVI_anomaly_interpsm.mat')
% load('DD_T2M_CRU_anomaly_interpsm.mat')
load('DD_T2M_ERA_anomaly_interpsm.mat')
load('DD_cru_cols.mat')
load('DD_cru_rows.mat')
load('DD_timev.mat')



%%%%%%% small bit to test conversion to mm/d units %%%%%%%
% also make positive

DD_dSM_dt_shallow_interpsm = DD_dSM_dt_shallow_interpsm .* -100 ;  % 100 mm soil depth like in the model

%%%%%%%




cd('F:\CRUJRA\tmp')
tmp_info    = ncinfo('crujra.v2.4.5d.tmp.1990.365d.noc.nc')    ;
lon_cru =  ncread('crujra.v2.4.5d.tmp.1990.365d.noc.nc','lon')    ; 
lat_cru =  ncread('crujra.v2.4.5d.tmp.1990.365d.noc.nc','lat')    ; 
time_cru = ncread('crujra.v2.4.5d.tmp.1990.365d.noc.nc','time')    ;

lat_cru = flipud(lat_cru) ; 

lat_cru = repmat(lat_cru,[1 720]) ; 
lon_cru = repmat(lon_cru',[360 1]) ; 

lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 

lats_matching_vector = NaN(360,1) ; 
lons_matching_vector = NaN(720,1) ; 



for r = 1:360

    curlat = lat_cru(r,1) ; 
    [mins locs] =  min(abs(lats_2_5(:,1) - curlat) )  ; 
    
    lats_matching_vector(r) = locs ;  
        
  
end

  for c = 1:720
      
    curlon = lon_cru(1,c) ; 
    [mins locs] =  min(abs(lons_2_5(1,:) - curlon) )  ; 
    
    lons_matching_vector(c) = locs ;  

  end 
 
lats_Cru_2_5_match = repmat(lats_matching_vector ,[1, 1440]) ; 
lons_Cru_2_5_match = repmat(lons_matching_vector' ,[720, 1]) ; 




% profile on 
%  r = 17 ;  c = 75
for r = 1:size(lats_2_5,1)
    r
    for c = 1:size(lats_2_5,2)
        
%         if (isnan(mask_mean(r,c)))
%             continue
%         end
        
        % find matching EASE row and cols
        [EASE_row, ~ ]= find(lats_Cru_2_5_match(:,1) == r) ; 
        [~, EASE_col ]= find(lons_Cru_2_5_match(1,:) == c) ;  
            

        row_col_extract = DD_cru_cols > min(EASE_col) & DD_cru_cols < max(EASE_col) ...
            & DD_cru_rows > min(EASE_row) & DD_cru_rows < max(EASE_row) ; 
        
        if(all(~row_col_extract))
            continue
        end
       
        


        DD_dSM_dt_shallow_interpsm_spatial = DD_dSM_dt_shallow_interpsm(row_col_extract,:)   ;    
        % DD_dSM_dt_deep_interpsm_spatial = DD_dSM_dt_deep_interpsm(row_col_extract,:)   ;        
        DD_NDVI_anomaly_interpsm_spatial = DD_NDVI_anomaly_interpsm(row_col_extract,:)   ;  
        % DD_f_E_anomaly_interpsm_spatial = DD_f_E_anomaly_interpsm(row_col_extract,:)   ;         
        % DD_T2M_CRU_anomaly_interpsm_spatial = DD_T2M_CRU_anomaly_interpsm(row_col_extract,:)   ;  
         DD_T2M_ERA_anomaly_interpsm_spatial = DD_T2M_ERA_anomaly_interpsm(row_col_extract,:)   ;        
        DD_timev_spatial = DD_timev(row_col_extract,:)   ;          

 % save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d\','DD_dSM_dt_shallow_interpsm_spatial',num2str(r),'_',num2str(c))  ,'DD_dSM_dt_shallow_interpsm_spatial','-v7.3')
 % save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d\','DD_NDVI_anomaly_interpsm_spatial',num2str(r),'_',num2str(c))  ,'DD_NDVI_anomaly_interpsm_spatial','-v7.3')
 % % save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\','DD_T2M_CRU_anomaly_interpsm_spatial',num2str(r),'_',num2str(c))  ,'DD_T2M_CRU_anomaly_interpsm_spatial','-v7.3')
 % save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d\','DD_T2M_ERA_anomaly_interpsm_spatial',num2str(r),'_',num2str(c))  ,'DD_T2M_ERA_anomaly_interpsm_spatial','-v7.3')
 % save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d\','DD_timev_spatial',num2str(r),'_',num2str(c))  ,'DD_timev_spatial','-v7.3')


 % save to test here
 save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\','DD_dSM_dt_shallow_interpsm_spatial',num2str(r),'_',num2str(c))  ,'DD_dSM_dt_shallow_interpsm_spatial','-v7.3')
 save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\','DD_NDVI_anomaly_interpsm_spatial',num2str(r),'_',num2str(c))  ,'DD_NDVI_anomaly_interpsm_spatial','-v7.3')
 % save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\','DD_T2M_CRU_anomaly_interpsm_spatial',num2str(r),'_',num2str(c))  ,'DD_T2M_CRU_anomaly_interpsm_spatial','-v7.3')
 save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\','DD_T2M_ERA_anomaly_interpsm_spatial',num2str(r),'_',num2str(c))  ,'DD_T2M_ERA_anomaly_interpsm_spatial','-v7.3')
 save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\','DD_timev_spatial',num2str(r),'_',num2str(c))  ,'DD_timev_spatial','-v7.3')


                      
     % c  
    end
    r
end








clear

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\spatial_2_5')
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_06_mm_d\spatial_2_5')
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_07_mm_d_low_SM\spatial_2_5_mm_d')
cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5')



DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
DD_T2M_anomaly_files = string(ls('*T2M_ERA_anomaly*')) ; 
DD_f_E_anomaly_files = string(ls('*NDVI_anomaly*')) ; 
sminterp = linspace(0,1,50) ; 


DD_dSM_dt_NDVI_gradient_median_05 = NaN(72,144) ; 
DD_dSM_dt_T2M_gradient_median_05 = NaN(72,144) ; 
DD_drydown_sampling_zeppe_05  = NaN(72,144) ; 
DD_dSM_dt_loss_rate_2_5_05   = NaN(72,144,60) ; 


% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_03')
% load('DD_dSM_dt_shallow_interpsm.mat')
% mean_DD_dSM_dt_global = median(DD_dSM_dt_shallow_interpsm,1,'omitnan') ; 
% clear DD_dSM_dt_shallow_interpsm

lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 


% use;  spatial = 23 ;
for spatial = 1:length(DD_dSM_dt_shallow_files)

    dummy_dSM_dt_shallow = load(DD_dSM_dt_shallow_files(spatial)) ; 
    dummy_f_E_anomaly = load(DD_f_E_anomaly_files(spatial)) ;     
     dummy_T2M_anomaly = load(DD_T2M_anomaly_files(spatial)) ;    

    dummy_dSM_dt_shallow = dummy_dSM_dt_shallow.DD_dSM_dt_shallow_interpsm_spatial ;
    dummy_f_E_anomaly = dummy_f_E_anomaly.DD_NDVI_anomaly_interpsm_spatial ;
    dummy_T2M_anomaly = dummy_T2M_anomaly.DD_T2M_ERA_anomaly_interpsm_spatial ;       


  

        if size(dummy_dSM_dt_shallow,1) < 300
            continue
        end

    % remove high SM above 0.4 for slope analyses. Keep it for the average
    % loss rate
    dummy_dSM_dt_shallow(:,41:end) = [] ; 
    dummy_f_E_anomaly(:,41:end) = [] ;
    dummy_T2M_anomaly(:,41:end) = [] ;



    % binning based on f_E NDVI
    dSM_dt_NDVI_binning = NaN(20,20) ; 
    dSM_dt_NDVI_sampling = NaN(20,20) ;

     dSM_dt_T2M_binning = NaN(20,20) ; 
     dSM_dt_T2M_sampling = NaN(20,20) ;

    counter_NDVI = NaN(20,1) ; 
    counter_NDVI(:,1) = 1 ; 

    counter_T2M = NaN(20,1) ; 
    counter_T2M(:,1) = 1 ; 

    % get 2.5 degree rowcol from name
    name_files = char(DD_dSM_dt_shallow_files(spatial)) ; 
    % grep row and col
    pat = digitsPattern;
    row_col_name = extract(name_files,pat) ; 
    row_2_5 = str2double(row_col_name{1}) ; 
    col_2_5 = str2double(row_col_name{2}) ;  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new condition to check if too northern
    lat_2_5_dummy = lats_2_5(row_2_5,col_2_5) ; 
    lon_2_5_dummy = lons_2_5(row_2_5,col_2_5) ; 

    if lat_2_5_dummy > 70
     continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    f_E_bins = linspace(0,1,21) ; 
    SM_bins  = linspace(0,0.6,21) ; 
    f_E_anomaly_bins = linspace(-0.2,0.2,21) ; 
    T2M_anomaly_bins = linspace(-10,10,21) ;    

    dummy_dSM_dt_shallow_ref = median(dummy_dSM_dt_shallow,1,'omitnan') ; 



for sm_bin = 1:20

    sm_bin_s = SM_bins(sm_bin) ;
    sm_bin_e = SM_bins(sm_bin+1) ;
    index_sm = sminterp > sm_bin_s & sminterp < sm_bin_e ; 
    % this is if original sminterp is same size as binning here
    % index_sm = sm_bin ; 

    % now calc anomaly of loss rate
     dSM_dt_shallow_extract = dummy_dSM_dt_shallow(:,index_sm)   - dummy_dSM_dt_shallow_ref(index_sm)  ;  
    % dSM_dt_shallow_extract = dummy_dSM_dt_shallow(:,index_sm)   -  mean_DD_dSM_dt_global(index_sm) ; 
    f_E_anomaly_extract  = dummy_f_E_anomaly(:,index_sm) ;  
    % T2M_anomaly_extract  = dummy_T2M_anomaly(:,index_sm) ;  


    for f_E_bin = 1:20


            f_E_anomaly_bin_s = f_E_anomaly_bins(f_E_bin) ;
            f_E_anomaly_bin_e = f_E_anomaly_bins(f_E_bin+1) ;


            % standard bining 
            % extract_index = f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e ;
              extract_index = f_E_anomaly_extract > f_E_anomaly_bin_s & f_E_anomaly_extract < f_E_anomaly_bin_e ;         

            dSM_dt_extract = dSM_dt_shallow_extract(extract_index) ; 

            % sampling
            dSM_dt_NDVI_sampling(f_E_bin,sm_bin) = sum(~isnan(dSM_dt_extract)) ; 
            % get dSM/dt anomalies
            dSM_dt_NDVI_binning(f_E_bin,sm_bin) = median(dSM_dt_extract); 


        % pcolor(SM_bins(2:end) - diff(SM_bins(1:2))/2,   f_E_bins(2:end) - diff(f_E_bins(1:2))/2, (dSM_dt_NDVI_binning)) ;
        % pcolor(SM_bins(2:end) - diff(SM_bins(1:2))/2,f_E_anomaly_bins(2:end) - diff(f_E_anomaly_bins(1:2))/2, (dSM_dt_NDVI_binning)) ; colormap(redblue_color) ; clim([-4e-3 4e-3])

    end

end


%% add same binning but for T2M

for sm_bin = 1:20

    sm_bin_s = SM_bins(sm_bin) ;
    sm_bin_e = SM_bins(sm_bin+1) ;
    index_sm = sminterp > sm_bin_s & sminterp < sm_bin_e ; 
    % this is if original sminterp is same size as binning here
    % index_sm = sm_bin ; 

    % now calc anomaly of loss rate
    dSM_dt_shallow_extract = dummy_dSM_dt_shallow(:,index_sm)   - dummy_dSM_dt_shallow_ref(index_sm)  ;  
    % f_E_anomaly_extract  = dummy_f_E_anomaly(:,index_sm) ;  
    T2M_anomaly_extract  = dummy_T2M_anomaly(:,index_sm) ;  


    for T2M_bin = 1:20


            T2M_anomaly_bin_s = T2M_anomaly_bins(T2M_bin) ;
            T2M_anomaly_bin_e = T2M_anomaly_bins(T2M_bin+1) ;


            % standard bining 
            % extract_index = f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e ;
              extract_index = T2M_anomaly_extract > T2M_anomaly_bin_s & T2M_anomaly_extract < T2M_anomaly_bin_e ;         

            dSM_dt_extract = dSM_dt_shallow_extract(extract_index) ; 

            % sampling
            dSM_dt_T2M_sampling(T2M_bin,sm_bin) = sum(~isnan(dSM_dt_extract)) ; 
            % get dSM/dt anomalies
            dSM_dt_T2M_binning(T2M_bin,sm_bin) = median(dSM_dt_extract); 


        % pcolor(SM_bins(2:end) - diff(SM_bins(1:2))/2,   f_E_bins(2:end) - diff(f_E_bins(1:2))/2, (dSM_dt_NDVI_binning)) ;
       % pcolor(SM_bins(2:end) - diff(SM_bins(1:2))/2, f_E_anomaly_bins(2:end) - diff(f_E_anomaly_bins(1:2))/2,(dSM_dt_NDVI_binning)) ; colormap(redblue_color) ; 

    end

end





%%


 dSM_dt_NDVI_binning(dSM_dt_NDVI_sampling < 10) = NaN ; 
 dSM_dt_T2M_binning(dSM_dt_T2M_sampling < 10) = NaN ; 
% dSM_dt_NDVI_binning = wiener2(dSM_dt_NDVI_binning) ; 


[FX_NDVI, FY_NDVI] =  gradient(dSM_dt_NDVI_binning) ; 
 [FX_T2M, FY_T2M] =  gradient(dSM_dt_T2M_binning) ; 

DD_dSM_dt_NDVI_gradient_median_05(row_2_5,col_2_5) = median(FY_NDVI(:),'omitnan') ; 
 DD_dSM_dt_T2M_gradient_median_05(row_2_5,col_2_5) =  median(FY_T2M(:),'omitnan') ; 

DD_drydown_sampling_zeppe_05(row_2_5,col_2_5) =  size(dummy_dSM_dt_shallow,1) ;
DD_dSM_dt_loss_rate_2_5_05(row_2_5,col_2_5,1:40) = dummy_dSM_dt_shallow_ref ; 


spatial
end



DD_dSM_dt_NDVI_gradient_zeppe_ERA = DD_dSM_dt_NDVI_gradient_median_05 ; 
DD_dSM_dt_T2M_gradient_zeppe_ERA = DD_dSM_dt_T2M_gradient_median_05 ; 
DD_drydown_sampling_zeppe_ERA = DD_drydown_sampling_zeppe_05 ;
DD_dSM_dt_loss_rate_2_5_zeppe_ERA = DD_dSM_dt_loss_rate_2_5_05 ; 



save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\DD_dSM_dt_NDVI_gradient_zeppe_ERA','DD_dSM_dt_NDVI_gradient_zeppe_ERA') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\DD_dSM_dt_T2M_gradient_zeppe_ERA','DD_dSM_dt_T2M_gradient_zeppe_ERA') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\DD_drydown_sampling_zeppe_ERA','DD_drydown_sampling_zeppe_ERA') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\DD_dSM_dt_loss_rate_2_5_zeppe_ERA','DD_dSM_dt_loss_rate_2_5_zeppe_ERA') ; 

% load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\DD_dSM_dt_NDVI_gradient_zeppe_ERA')


Fig_Panel = figure('units','centimeters','position',[10 3 39 20])  ;
 h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, median(DD_dSM_dt_loss_rate_2_5_zeppe_ERA,3,'omitnan')); 
set(h,'LineStyle','none')
shading flat
colormap(flipud(bam_color(1:128,:))) 
clim([0 5]) ; 
cbr = colorbar ;  ylabel(cbr,'SM loss [mm/d]','FontSize',18)
hold on ; plot(CoastlineLon, CoastlineLat,'Color','k');
fontsize(Fig_Panel,15,'points')


cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')


lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 

 % h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = DD_dSM_dt_NDVI_gradient_median_05 ; 
%xmap = DD_dSM_dt_T2M_gradient_median_05 ; 

Fig_Panel = figure('units','centimeters','position',[10 3 39 20])  ;
 h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
colormap(redblue_color) 
% clim([-5e-4  5e-4])
clim([-0.1 0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss/\Deltaf_E [m³/m³/day/f_E]','FontSize',18)
% title('negative NDVI anomaly')
xlabel('longitude')
ylabel('latitude')
plot(CoastlineLon, CoastlineLat,'Color','k');
set(gca,'FontSize',15)


% % add arrow with description do it properly now based on position of axes
axes_diff_Position = get(gca, 'Position');
% calc positions
%  4.8682    2.2003   29.0217   16.3020

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow1,'Position',[4.8682+29.0217+3    2.2003+(16.3020)/2+1   0     (16.3020)/2-1  ]) ; 
 set(arrow2,'Position',[4.8682+29.0217+3    2.2003+(16.3020)/2-1     0  -(16.3020)/2+1 ]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[4.8682+29.0217+4.2    2.2003+(16.3020)/2+4.8   10     1  ]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[4.8682+29.0217+4.2   2.2003+(16.3020)/2-1     10  1 ]) ; 
set(gca,'Box','on');



% saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model\SM_loss_d_f_E_global_map','svg')
close 







%% use spatial to bin into 1 phase space anomaly relative to 2.5 loss rate
%% bin into NDVI anomalies and SM

clear

cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')


% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\spatial_2_5\')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_06_mm_d\spatial_2_5')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_07_mm_d_low_SM\spatial_2_5_mm_d')
cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')

DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
DD_T2M_anomaly_files = string(ls('*T2M_ERA_anomaly*')) ; 
DD_f_E_anomaly_files = string(ls('*NDVI_anomaly*')) ; 
DD_timev_files = string(ls('*timev*')) ; 


sminterp = linspace(0,0.6,60) ; 
f_E_bins = linspace(-0.2,0.2,41) ; 


lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 

dSM_dt_NDVI_binning = NaN(40,40,250000) ; 
counter =  NaN(40,40) ; 
counter(:,:) = 1 ; 


% load to check if pos or negative
load('F:\projects\SM_long_term_DDs\data_for_figures\DD_dSM_dt_NDVI_gradient_median_05.mat')
load('F:\projects\SM_long_term_DDs\data_for_figures\DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
DD_dSM_dt_NDVI_gradient_median_05(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% DD_dSM_dt_NDVI_gradient_median_05(DD_dSM_dt_NDVI_gradient_median_05 < 0) = NaN ; 

load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\DD_dSM_dt_NDVI_gradient_zeppe_ERA')
DD_dSM_dt_NDVI_gradient_zeppe_ERA(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 

% imagesc(DD_dSM_dt_NDVI_gradient_median_05) ; colormap(redblue_color) ; clim([-1e-3 1e-3])
% [xs, ys] = getpts() ;  xs  = round(xs) ;   ys  = round(ys) ; 
% close

xs = [1 144] ; ys = [1 72] ; 

%   spatial = 23 ;
for spatial = 1:length(DD_dSM_dt_shallow_files)

    % cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\spatial_2_5\')
    % cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\')
    % cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_07_mm_d_low_SM\spatial_2_5_mm_d')
    cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')

    dummy_dSM_dt_shallow = load(DD_dSM_dt_shallow_files(spatial)) ; 
    dummy_f_E_anomaly = load(DD_f_E_anomaly_files(spatial)) ;     
  

    dummy_dSM_dt_shallow = dummy_dSM_dt_shallow.DD_dSM_dt_shallow_interpsm_spatial ;
    dummy_f_E_anomaly = dummy_f_E_anomaly.DD_NDVI_anomaly_interpsm_spatial ;


    if size(dummy_dSM_dt_shallow,1) < 300
        continue
    end

    % get rid of SM higher than 0.4 so so it is comparable with SMAP
    dummy_dSM_dt_shallow(:,41:end) = [] ; 
    dummy_f_E_anomaly(:,41:end)  = [] ; 


    % get 2.5 degree rowcol from name
    name_files = char(DD_dSM_dt_shallow_files(spatial)) ; 
    % grep row and col
    pat = digitsPattern;
    row_col_name = extract(name_files,pat) ; 
    row_2_5 = str2double(row_col_name{1}) ; 
    col_2_5 = str2double(row_col_name{2}) ;   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new condition to check if too northern
    lat_2_5_dummy = lats_2_5(row_2_5,col_2_5) ; 
    lon_2_5_dummy = lons_2_5(row_2_5,col_2_5) ; 

    if lat_2_5_dummy > 70
     continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    % % % only analyze cover or extraction effect
    if DD_dSM_dt_NDVI_gradient_zeppe_ERA(row_2_5,col_2_5)  < 0 
        continue
    end

    % only analyze spatial subset
    if  ~(col_2_5 > xs(1) && col_2_5 < xs(2)  && row_2_5 > ys(1) && row_2_5 < ys(2))
        continue
    end




    median_spatial_loss = median(dummy_dSM_dt_shallow,1,'omitnan') ; 

    dSM_dt_shallow_anomaly = dummy_dSM_dt_shallow - median_spatial_loss ;  
    f_e_extract = dummy_f_E_anomaly;  
    % dSM_dt_deep_extract = dummy_DD_dSM_dt_deep(:,index_sm) ;     


    for f_E_bin = 1:40

            f_E_bin_s = f_E_bins(f_E_bin) ;
            f_E_bin_e = f_E_bins(f_E_bin+1) ;

            extract_row = any(f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e,2) ; 


            dSM_dt_shallow_extract_2 = dSM_dt_shallow_anomaly(extract_row,:) ; 
            f_e_extract_2 = f_e_extract(extract_row,:) ; 

            dSM_dt_shallow_extract_2(~(f_e_extract_2 > f_E_bin_s & f_e_extract_2 < f_E_bin_e)) = NaN ; 


    for j = 1:size(dSM_dt_shallow_extract_2,2)
        
        if all(isnan( dSM_dt_shallow_extract_2(:,j)))
            continue
        end

        for i =1:size(dSM_dt_shallow_extract_2,1)

        if (isnan( dSM_dt_shallow_extract_2(i,j)))
            continue
        end
    
    dSM_dt_NDVI_binning(f_E_bin,j,counter(f_E_bin,j)) = dSM_dt_shallow_extract_2(i,j); 

    counter(f_E_bin,j) = counter(f_E_bin,j) + 1 ;


    if counter(f_E_bin,j) == 250000
    return
    end

        end % i    
    end  % j 


    end
spatial
end



% visualize global phase space grow season Zeppe_model

sampling_2D_bin = sum(~isnan(dSM_dt_NDVI_binning),3) ; 
xmap = median(dSM_dt_NDVI_binning,3,'omitnan');
% xmap(sampling_2D_bin < 100) = NaN ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction_sampling = sum(~isnan(dSM_dt_NDVI_binning),3) ; 
dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction = median(dSM_dt_NDVI_binning,3,'omitnan');
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction_sampling','dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction_sampling') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction','dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction') ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp(1:40) - diff(sminterp(1:2))/2,   f_E_bins(2:end) - diff(f_E_bins(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.5 0.5]) %
% clim([-6e-4 0.0]) %
 ylim([-0.2 0.2])
 xlim([0 0.5])
 colormap(redblue_color) %
% colormap parula
cbr = colorbar ;
cbr.Label.String = "SM loss anomaly [mm/day]";
xlabel('SM [m³/m³]')
ylabel('f_E anomaly [-]')
title('Europe')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% fontsize(16,'points')


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\SM_loss_shallow_2D_SM_fE_anomalies_Europe_100samp','svg')
close 



save('F:\projects\SM_long_term_DDs\data_for_figures\DD_dSM_dt_T2M_gradient_zeppe_median','DD_dSM_dt_T2M_gradient_zeppe_median') ; 









%% do 2.5 or 5 degree but just mean DSM/dt map and mean loss rate curve. Comparison to observaitons
% multiply loss rate by 2, as dt was 12 hours. Ask Lucas about instability.
clear

cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_08\spatial_2_5\')


DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
DD_T2M_anomaly_files = string(ls('*T2M_CRU_anomaly*')) ; 
DD_f_E_anomaly_files = string(ls('*NDVI_anomaly*')) ; 
sminterp = linspace(0,1,50) ; 


DD_dSM_dt_zeppe_2_5_median = NaN(72,144) ; 
DD_dSM_dt_zeppe_2_5_mean = NaN(72,144) ; 
DD_drydown_sampling_zeppe  = NaN(72,144) ; 



% use;  spatial = 23 ;
for spatial = 1:length(DD_dSM_dt_shallow_files)

    dummy_dSM_dt_shallow = load(DD_dSM_dt_shallow_files(spatial)) ; 
    dummy_dSM_dt_shallow = dummy_dSM_dt_shallow.DD_dSM_dt_shallow_interpsm_spatial ;
    dummy_dSM_dt_shallow = dummy_dSM_dt_shallow .* 2 ; 

    % kick out high SM?
    % dummy_dSM_dt_shallow(:,40:end) = NaN ; 

       % if size(dummy_dSM_dt_shallow,1) < 300
           % continue
       % end


    name_files = char(DD_dSM_dt_shallow_files(spatial)) ; 
    % grep row and col
    pat = digitsPattern;
    row_col_name = extract(name_files,pat) ; 
    row_2_5 = str2double(row_col_name{1}) ; 
    col_2_5 = str2double(row_col_name{2}) ;   



DD_dSM_dt_zeppe_2_5_median(row_2_5,col_2_5) = median(median(dummy_dSM_dt_shallow,1,'omitnan'),'omitnan') ; 
DD_dSM_dt_zeppe_2_5_mean  (row_2_5,col_2_5) =  mean(mean(dummy_dSM_dt_shallow,1,'omitnan'),'omitnan') ; 
DD_drydown_sampling_zeppe(row_2_5,col_2_5) =  size(dummy_dSM_dt_shallow,1) ;


spatial
end




cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')


lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 



xmap = median(DD_dSM_dt_zeppe_2_5_median,3,'omitnan');
xmap(DD_drydown_sampling_zeppe < 300) = NaN ; 


Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
colormap(bam_color(1:128,:)) 
clim([-0.05 0])% dSM/dt diurnal
xticks(-160:80:160)
yticks(-80:40:80)
hcb2=colorbar;
set(hcb2, 'FontSize',16)
set(hcb2, 'FontSize',16,'YTick',-0.05:0.0125:0)
xlabel('longitude','FontSize',16)
ylabel('latitude','FontSize',16)
ylabel(hcb2,' median \DeltaSM/\Deltat model [m³/m³/day]','FontSize',16)
%ylabel(hcb2,'\DeltaSM/\Deltat change [%]','FontSize',16)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
% inset with pdf
histoaxes = axes('units','centimeters','Position',[ 4.5, 16.5-8-4.25+1, 3,  3]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(histoaxes, 'FontSize',12)
xlim([-0.05 0.0])




%% use 2_5 degree saves to bin into oone phase space but for growing season only 

% 
% 
% clear
% 
% 
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_03\spatial_2_5\')
% 
% 
% % calc global median curve
% load('F:\projects\SM_long_term_DDs\Zeppe_model_output_03\DD_dSM_dt_shallow_interpsm.mat')
% median_global_SM_loss = median(DD_dSM_dt_shallow_interpsm,1,'omitnan') ; 
% clear DD_dSM_dt_shallow_interpsm
% 
% DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
% DD_T2M_anomaly_files = string(ls('*T2M_CRU_anomaly*')) ; 
% DD_f_E_anomaly_files = string(ls('*NDVI_anomaly*')) ; 
% DD_timev_files = string(ls('*timev*')) ; 
% 
% sminterp = linspace(0,1,50) ; 
% f_E_bins = linspace(-0.2,0.2,50) ; 
% 
% cd('F:\AVHRR_phenology')
% load('Phenology_2D_EOS.mat')
% load('Phenology_2D_SOS.mat')
% ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
% Phenology_2D_EOS(Phenology_2D_EOS == 0) = NaN ; 
% Phenology_2D_SOS(Phenology_2D_SOS == 0) = NaN ; 
% 
% 
% % mask = t_start_end_2D_array(:,1) < 4080 | t_start_end_2D_array(:,2) > 13210;
% % rowcol_2D_array(mask,:) = [] ; 
% % t_start_end_2D_array(mask,:) = [] ; 
% 
% % cut to 1990-2014
% Phenology_2D_EOS = Phenology_2D_EOS(:,:,9:33) ; 
% Phenology_2D_SOS = Phenology_2D_SOS(:,:,9:33) ; 
% 
% Phenology_2D_EOS(Phenology_2D_EOS < 4080 | Phenology_2D_EOS > 13210) = NaN ; 
% Phenology_2D_SOS(Phenology_2D_SOS < 4080 | Phenology_2D_SOS > 13210) = NaN ; 
% 
% 
% DD_dSM_dt_NDVI_gradient_median = NaN(72,144) ; 
% DD_dSM_dt_T2M_gradient_median = NaN(72,144) ; 
% DD_drydown_sampling_zeppe  = NaN(72,144) ; 
% 
% 
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_03\spatial_2_5\')
% 
% 
% dSM_dt_NDVI_binning = NaN(50,50,50000) ; 
% counter =  NaN(50,50) ; 
% counter(:,:) = 1 ; 
% 
% 
% %   spatial = 23 ;
% for spatial = 1:length(DD_dSM_dt_shallow_files)
% 
%     cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_03\spatial_2_5\')
%     dummy_dSM_dt_shallow = load(DD_dSM_dt_shallow_files(spatial)) ; 
%     dummy_f_E_anomaly = load(DD_f_E_anomaly_files(spatial)) ;     
%     dummy_T2M_anomaly = load(DD_T2M_anomaly_files(spatial)) ;    
%     dummy_timev = load(DD_timev_files(spatial)) ; 
% 
%     dummy_dSM_dt_shallow = dummy_dSM_dt_shallow.DD_dSM_dt_shallow_interpsm_spatial ;
%     dummy_f_E_anomaly = dummy_f_E_anomaly.DD_NDVI_anomaly_interpsm_spatial ;
%     dummy_T2M_anomaly = dummy_T2M_anomaly.DD_T2M_CRU_anomaly_interpsm_spatial ;    
% 
%     dummy_timev = dummy_timev.DD_timev_spatial ;    
% 
% 
%     % get 2.5 degree rowcol from name
%     name_files = char(DD_dSM_dt_shallow_files(spatial)) ; 
%     % grep row and col
%     pat = digitsPattern;
%     row_col_name = extract(name_files,pat) ; 
%     row_2_5 = str2double(row_col_name{1}) ; 
%     col_2_5 = str2double(row_col_name{2}) ;   
% 
% 
%     % indexing in growing and non growing season
%     EOS_dummy = squeeze(Phenology_2D_EOS(row_2_5,col_2_5,:)) ; 
%     SOS_dummy = squeeze(Phenology_2D_SOS(row_2_5,col_2_5,:)) ; 
% 
%     % first criteria is that total sample number is over 10k
%     if (all(isnan(EOS_dummy)) || all(isnan(SOS_dummy)))
%       continue  
%     end
% 
% 
% 
% 
%  % build logical indices based on histcounts                
% cd('E:\Daten Baur\Matlab code')       
% dummy_t_mean = mean(dummy_timev,2,'omitnan') ; 
% [left, right] = IntervalUnion(SOS_dummy, EOS_dummy); % FEX
% edges = [left,right]';
% edges = edges(:);
% [~,~,loc] = histcounts(dummy_t_mean,edges);
% L = mod(loc,2)==1;          
% not_L = ~L ;                 
% 
% % get growing season
% dummy_dSM_dt_gseason =    dummy_dSM_dt_shallow(L,:) ;             
% dummy_dSM_dt_non_gseason =    dummy_dSM_dt_shallow(not_L,:) ;  
% 
% dummy_f_E_anomaly_gseason = dummy_f_E_anomaly(L,:) ;   
% dummy_f_E_anomaly_non_gseason = dummy_f_E_anomaly(not_L,:) ; 
% 
% 
%     dSM_dt_shallow_extract = dummy_dSM_dt_non_gseason - median_global_SM_loss ;  
%     f_e_extract = dummy_f_E_anomaly_non_gseason;  
%     % dSM_dt_deep_extract = dummy_DD_dSM_dt_deep(:,index_sm) ;     
% 
% 
%     for f_E_bin = 1:49
% 
%             f_E_bin_s = f_E_bins(f_E_bin) ;
%             f_E_bin_e = f_E_bins(f_E_bin+1) ;
% 
%             extract_row = any(f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e,2) ; 
% 
% 
%             dSM_dt_shallow_extract_2 = dSM_dt_shallow_extract(extract_row,:) ; 
%             f_e_extract_2 = f_e_extract(extract_row,:) ; 
% 
%             dSM_dt_shallow_extract_2(~(f_e_extract_2 > f_E_bin_s & f_e_extract_2 < f_E_bin_e)) = NaN ; 
% 
% 
%     for j = 1:size(dSM_dt_shallow_extract_2,2)
% 
%         if all(isnan( dSM_dt_shallow_extract_2(:,j)))
%             continue
%         end
% 
%         for i =1:size(dSM_dt_shallow_extract_2,1)
% 
%         if (isnan( dSM_dt_shallow_extract_2(i,j)))
%             continue
%         end
% 
%     dSM_dt_NDVI_binning(f_E_bin,j,counter(f_E_bin,j)) = dSM_dt_shallow_extract_2(i,j); 
% 
%     counter(f_E_bin,j) = counter(f_E_bin,j) + 1 ;
% 
% 
%     if counter(f_E_bin,j) == 50000
%     return
%     end
% 
%         end % i    
%     end  % j 
% 
% 
% 
% 
% 
% 
% 
%     end
% spatial
% end
% 
% 
% 
% % visualize global phase space grow season Zeppe_model
% 
% sampling_2D_bin = sum(~isnan(dSM_dt_NDVI_binning),3) ; 
% xmap = median(dSM_dt_NDVI_binning,3,'omitnan');
% xmap(sampling_2D_bin < 1000) = NaN ; 
% 
% 
% 
% 
% cd('E:\Daten Baur\Matlab code\redblue')
% redblue_color = redblue(100) ; 
% bam_color = crameri('bam') ;
% load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')
% 
% 
% 
% Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
% h1 = pcolor(sminterp(1:end) - diff(sminterp(1:2))/2,   f_E_bins(1:end) - diff(f_E_bins(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
% clim([-4e-3 4e-3]) %
% % clim([-6e-4 0.0]) %
%  ylim([-0.2 0.2])
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





%% spatial bin into 2D phase space for temperature


clear

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\spatial_2_5\')
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08_mm_d_low_SM\spatial_2_5_mm_d')
cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')

DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
DD_T2M_anomaly_files = string(ls('*T2M_ERA_anomaly*')) ; 
DD_f_E_anomaly_files = string(ls('*NDVI_anomaly*')) ; 
DD_timev_files = string(ls('*timev*')) ; 


sminterp = linspace(0,0.6,60) ; 
f_E_bins = linspace(-0.2,0.2,41) ; 
T2M_bins = linspace(-15,15,41) ; 



lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 

dSM_dt_T2M_binning = NaN(40,40,250000) ; 
counter =  NaN(40,40) ; 
counter(:,:) = 1 ; 


%   spatial = 23 ;
for spatial = 1:length(DD_dSM_dt_shallow_files)

    % cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\spatial_2_5\')
 cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')

    dummy_dSM_dt_shallow = load(DD_dSM_dt_shallow_files(spatial)) ; 
    dummy_f_E_anomaly = load(DD_f_E_anomaly_files(spatial)) ;     
    dummy_T2M_anomaly = load(DD_T2M_anomaly_files(spatial)) ;       
  

    dummy_dSM_dt_shallow = dummy_dSM_dt_shallow.DD_dSM_dt_shallow_interpsm_spatial ;
    dummy_f_E_anomaly = dummy_f_E_anomaly.DD_NDVI_anomaly_interpsm_spatial ;
    dummy_T2M_anomaly = dummy_T2M_anomaly.DD_T2M_ERA_anomaly_interpsm_spatial ;


    if size(dummy_dSM_dt_shallow,1) < 300
        continue
    end

    dummy_dSM_dt_shallow(:,41:end) = [] ; 
    dummy_f_E_anomaly(:,41:end) = [] ; 
    dummy_T2M_anomaly(:,41:end) = [] ; 


    % get 2.5 degree rowcol from name
    name_files = char(DD_dSM_dt_shallow_files(spatial)) ; 
    % grep row and col
    pat = digitsPattern;
    row_col_name = extract(name_files,pat) ; 
    row_2_5 = str2double(row_col_name{1}) ; 
    col_2_5 = str2double(row_col_name{2}) ;  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new condition to check if too northern
    lat_2_5_dummy = lats_2_5(row_2_5,col_2_5) ; 
    lon_2_5_dummy = lons_2_5(row_2_5,col_2_5) ; 

    if lat_2_5_dummy > 70
     continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    median_spatial_loss = median(dummy_dSM_dt_shallow,1,'omitnan') ; 

    dSM_dt_shallow_anomaly = dummy_dSM_dt_shallow - median_spatial_loss ;  
    f_e_extract = dummy_T2M_anomaly;  
    % dSM_dt_deep_extract = dummy_DD_dSM_dt_deep(:,index_sm) ;     


    for T2M_bin = 1:40

            T2M_bin_s = T2M_bins(T2M_bin) ;
            T2M_bin_e = T2M_bins(T2M_bin+1) ;

            extract_row = any(f_e_extract > T2M_bin_s & f_e_extract < T2M_bin_e,2) ; 


            dSM_dt_shallow_extract_2 = dSM_dt_shallow_anomaly(extract_row,:) ; 
            f_e_extract_2 = f_e_extract(extract_row,:) ; 

            dSM_dt_shallow_extract_2(~(f_e_extract_2 > T2M_bin_s & f_e_extract_2 < T2M_bin_e)) = NaN ; 


    for j = 1:size(dSM_dt_shallow_extract_2,2)
        
        if all(isnan( dSM_dt_shallow_extract_2(:,j)))
            continue
        end

        for i =1:size(dSM_dt_shallow_extract_2,1)

        if (isnan( dSM_dt_shallow_extract_2(i,j)))
            continue
        end
    
    dSM_dt_T2M_binning(T2M_bin,j,counter(T2M_bin,j)) = dSM_dt_shallow_extract_2(i,j); 

    counter(T2M_bin,j) = counter(T2M_bin,j) + 1 ;


    if counter(T2M_bin,j) == 250000
    return
    end

        end % i    
    end  % j 


    end
spatial
end



% visualize global phase space grow season Zeppe_model

sampling_2D_bin = sum(~isnan(dSM_dt_T2M_binning),3) ; 
xmap = median(dSM_dt_T2M_binning,3,'omitnan');
% xmap(sampling_2D_bin < 1000) = NaN ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dSM_dt_T2M_anomaly_2D_zeppe_ERA_sampling = sum(~isnan(dSM_dt_T2M_binning),3) ; 
dSM_dt_T2M_anomaly_2D_zeppe_ERA = median(dSM_dt_T2M_binning,3,'omitnan');
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_T2M_anomaly_2D_zeppe_ERA_sampling','dSM_dt_T2M_anomaly_2D_zeppe_ERA_sampling') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_T2M_anomaly_2D_zeppe_ERA','dSM_dt_T2M_anomaly_2D_zeppe_ERA') ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp(1:40) - diff(sminterp(1:2))/2,   T2M_bins(1:end) - diff(T2M_bins(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.5 0.5]) %
% clim([-6e-4 0.0]) %
 ylim([-10 10])
 xlim([0 0.5])
 colormap(redblue_color) %
% colormap parula
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat shallow anomaly [m³/m³/day]";
xlabel('SM [m³/m³]')
ylabel('f_E anomaly [-]')
title('Europe')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% fontsize(16,'points')


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\SM_loss_shallow_2D_SM_fE_anomalies_Europe_100samp','svg')
close 



save('F:\projects\SM_long_term_DDs\data_for_figures\DD_dSM_dt_T2M_gradient_zeppe_median','DD_dSM_dt_T2M_gradient_zeppe_median') ; 





%% spatial bin fixed SM into NDVI and T2M anomalies
clear

% cd('E:\Daten Baur\Matlab code\redblue')
% redblue_color = redblue(100) ; 
% bam_color = crameri('bam') ;
% load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')


cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')

DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
DD_T2M_anomaly_files = string(ls('*T2M_ERA_anomaly*')) ; 
DD_f_E_anomaly_files = string(ls('*NDVI_anomaly*')) ; 
DD_timev_files = string(ls('*timev*')) ; 


sminterp = linspace(0,0.6,60) ; 
f_E_bins = linspace(-0.2,0.2,41) ; 
T2M_bins = linspace(-15,15,41) ; 


lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 


dSM_dt_NDVI_binning = NaN(40,40,250000) ; 
counter =  NaN(40,40) ; 
counter(:,:) = 1 ; 


% load to check if pos or negative
load('F:\projects\SM_long_term_DDs\data_for_figures\DD_dSM_dt_NDVI_gradient_median_05.mat')
load('F:\projects\SM_long_term_DDs\data_for_figures\DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
DD_dSM_dt_NDVI_gradient_median_05(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% DD_dSM_dt_NDVI_gradient_median_05(DD_dSM_dt_NDVI_gradient_median_05 < 0) = NaN ; 

% imagesc(DD_dSM_dt_NDVI_gradient_median_05) ; colormap(redblue_color) ; clim([-1e-3 1e-3])
% [xs, ys] = getpts() ;  xs  = round(xs) ;   ys  = round(ys) ; 
% close
% global settings
xs = [1 144] ; ys = [1 72] ; 




%   spatial = 23 ;
for spatial = 1:length(DD_dSM_dt_shallow_files)

cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')
    dummy_dSM_dt_shallow = load(DD_dSM_dt_shallow_files(spatial)) ; 
    dummy_f_E_anomaly = load(DD_f_E_anomaly_files(spatial)) ;  
    dummy_T2M_anomaly = load(DD_T2M_anomaly_files(spatial)) ;    

  
    dummy_dSM_dt_shallow = dummy_dSM_dt_shallow.DD_dSM_dt_shallow_interpsm_spatial ;
    dummy_f_E_anomaly = dummy_f_E_anomaly.DD_NDVI_anomaly_interpsm_spatial ;
    dummy_T2M_anomaly = dummy_T2M_anomaly.DD_T2M_ERA_anomaly_interpsm_spatial  ;


    dummy_dSM_dt_shallow(:,41:end) = [] ; 
    dummy_f_E_anomaly(:,41:end) = [] ; 
    dummy_T2M_anomaly(:,41:end) = [] ; 
    sminterp(41:end) = [] ; 

    % if size(dummy_dSM_dt_shallow,1) < 300
    %     continue
    % end

    % get 2.5 degree rowcol from name
    name_files = char(DD_dSM_dt_shallow_files(spatial)) ; 
    % grep row and col
    pat = digitsPattern;
    row_col_name = extract(name_files,pat) ; 
    row_2_5 = str2double(row_col_name{1}) ; 
    col_2_5 = str2double(row_col_name{2}) ; 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new condition to check if too northern
    lat_2_5_dummy = lats_2_5(row_2_5,col_2_5) ; 
    lon_2_5_dummy = lons_2_5(row_2_5,col_2_5) ; 

    if lat_2_5_dummy > 70
     continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % % only analyze cover or extraction effect
    % if DD_dSM_dt_NDVI_gradient_median_05(row_2_5,col_2_5)  < 0 
    %     continue
    % end


    % % only analyze spatial subset
    % if  ~(col_2_5 > xs(1) && col_2_5 < xs(2)  && row_2_5 > ys(1) && row_2_5 < ys(2))
    %     continue
    % end


    % only analyze a fixed range of SM 0.2-0.3    
    dummy_dSM_dt_shallow =  dummy_dSM_dt_shallow(:,sminterp > 0.3 & sminterp < 0.4) ;
    dummy_f_E_anomaly    =  dummy_f_E_anomaly(:,sminterp > 0.3 & sminterp < 0.4) ;
    dummy_T2M_anomaly    =  dummy_T2M_anomaly(:,sminterp > 0.3 & sminterp < 0.4) ;

                

    % now just 1 valie
    median_spatial_loss = median(dummy_dSM_dt_shallow,1,'omitnan') ; 
    dSM_dt_shallow_anomaly = dummy_dSM_dt_shallow - median_spatial_loss ;  




    for f_E_bin = 1:40

            f_E_bin_s = f_E_bins(f_E_bin) ;
            f_E_bin_e = f_E_bins(f_E_bin+1) ;

            for T2M_bin = 1:40

               T2M_bin_s = T2M_bins(T2M_bin) ;
               T2M_bin_e = T2M_bins(T2M_bin+1) ;

               extract_dummy = (dummy_f_E_anomaly > f_E_bin_s & dummy_f_E_anomaly < f_E_bin_e & ...
                          dummy_T2M_anomaly > T2M_bin_s & dummy_T2M_anomaly < T2M_bin_e  ) ;

               dSM_dt_extract = dSM_dt_shallow_anomaly(extract_dummy)  ;


               dSM_dt_NDVI_binning(f_E_bin,T2M_bin,counter(f_E_bin,T2M_bin):counter(f_E_bin,T2M_bin) + length(dSM_dt_extract) -1) = dSM_dt_extract ; 

               counter(f_E_bin,T2M_bin) =  counter(f_E_bin,T2M_bin) + length(dSM_dt_extract) ; 

               if counter(f_E_bin,T2M_bin) > 250000
                   return
               end


            end
    end



spatial


end


dSM_dt_NDVI_binning(dSM_dt_NDVI_binning == 0) = NaN ; 
sampling_2D_bin = sum(~isnan(dSM_dt_NDVI_binning),3) ; 
xmap = median(dSM_dt_NDVI_binning,3,'omitnan');
% xmap(sampling_2D_bin < 100) = NaN ; 
% xmap(sampling_2D_bin < 1000) = NaN ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA_sampling = sum(~isnan(dSM_dt_NDVI_binning),3) ; 
dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA = median(dSM_dt_NDVI_binning,3,'omitnan');
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA_sampling','dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA_sampling') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA','dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA') ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(T2M_bins(1:end) - diff(T2M_bins(1:2))/2,   f_E_bins(1:end) - diff(f_E_bins(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-1 1]) %
% clim([-6e-4 0.0]) %
 ylim([-0.2 0.2])
 xlim([-15 15])
 colormap(redblue_color) %
% colormap parula
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat shallow anomaly [m³/m³/day]";
xlabel('T2M anomaly [°K]')
ylabel('NDVI anomaly [-]')
% title('Europe')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% fontsize(16,'points')

saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_04\_dSM_dt_anomaly_T2M_NDVI_2D','svg')
close 





%% do spatial anaysis for single phase space. Try linear correction for temperature effects
% temp is a problem especially for the model. Seems to be less problem for
% SMAP. Maybe a simple linear correction of Temp effects on dSM/dt per 2.5
% degree box is enough to help.
% clear
% 
% 
% cd('E:\Daten Baur\Matlab code\redblue')
% redblue_color = redblue(100) ; 
% bam_color = crameri('bam') ;
% load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')
% 
% 
% 
% 
%  cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\spatial_2_5\')
% % cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\')
% 
% DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
% DD_T2M_anomaly_files = string(ls('*T2M_CRU_anomaly*')) ; 
% DD_f_E_anomaly_files = string(ls('*NDVI_anomaly*')) ; 
% DD_timev_files = string(ls('*timev*')) ; 
% 
% 
% sminterp = linspace(0,1,50) ; 
% sminterp_volumetric = linspace(0,0.5,50) ; 
% f_E_bins = linspace(-0.2,0.2,50) ; 
% 
% 
% 
% 
% dSM_dt_NDVI_binning = NaN(50,50,150000) ; 
% counter =  NaN(50,50) ; 
% counter(:,:) = 1 ; 
% f_E_T2M_corrcoef_vec = NaN(72,144) ;
% 
% 
% % load to check if pos or negative
% load('F:\projects\SM_long_term_DDs\data_for_figures\DD_dSM_dt_NDVI_gradient_median_05.mat')
% load('F:\projects\SM_long_term_DDs\data_for_figures\DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
% DD_dSM_dt_NDVI_gradient_median_05(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% % DD_dSM_dt_NDVI_gradient_median_05(DD_dSM_dt_NDVI_gradient_median_05 < 0) = NaN ; 
% 
% 
% % imagesc(DD_dSM_dt_NDVI_gradient_median_05) ; colormap(redblue_color) ; clim([-1e-3 1e-3])
% % [xs, ys] = getpts() ;  xs  = round(xs) ;   ys  = round(ys) ; 
% % close
% 
% xs = [1 144] ; ys = [1 72] ;
% 
% 
% 
% %   spatial = 23 ;
% for spatial = 1:length(DD_dSM_dt_shallow_files)
% 
%     % cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\spatial_2_5\')
%     cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\')
% 
% 
%     dummy_dSM_dt_shallow = load(DD_dSM_dt_shallow_files(spatial)) ; 
%     dummy_f_E_anomaly = load(DD_f_E_anomaly_files(spatial)) ;  
%     dummy_T2M_anomaly = load(DD_T2M_anomaly_files(spatial)) ;  
% 
% 
%     dummy_dSM_dt_shallow = dummy_dSM_dt_shallow.DD_dSM_dt_shallow_interpsm_spatial ;
%     dummy_f_E_anomaly = dummy_f_E_anomaly.DD_NDVI_anomaly_interpsm_spatial ;
%     dummy_T2M_anomaly = dummy_T2M_anomaly.DD_T2M_CRU_anomaly_interpsm_spatial ;
% 
%     % correlate dSM_dt and dummy T2M anomaly
%     % sminterp_array = repmat(sminterp,[size(dummy_dSM_dt_shallow,1),1]) ;
%     % scat = scatter(sminterp_array(:), dummy_dSM_dt_shallow(:),20,dummy_T2M_anomaly(:),'filled','MarkerFaceAlpha',0.5) ; 
%     %  scat = scatter(dummy_f_E_anomaly(:),dummy_T2M_anomaly(:),20,'filled','MarkerFaceAlpha',0.5) ; 
%     dummy_f_E_anomaly_vec = dummy_f_E_anomaly(:) ; 
%     dummy_T2M_anomaly_vec = dummy_T2M_anomaly(:) ; 
%     dummy_f_E_anomaly_vec(isnan(dummy_f_E_anomaly_vec)) = [] ; 
%     dummy_T2M_anomaly_vec(isnan(dummy_T2M_anomaly_vec)) = [] ; 
% 
%     corr_f_E_NDVI_anomalies = corrcoef(dummy_f_E_anomaly_vec,dummy_T2M_anomaly_vec) ;
% 
% for i = 1:50
% scatter(dummy_T2M_anomaly(:,i),dSM_dt_shallow_anomaly(:,i))
% hold on
% end
% 
% 
% 
% 
% % dummy_f_E_anomaly_standard = (dummy_f_E_anomaly - mean(dummy_f_E_anomaly(:),'omitnan')) ./ std(dummy_f_E_anomaly(:),'omitnan') ; 
% % dummy_T2M_anomaly_standard = (dummy_T2M_anomaly - mean(dummy_T2M_anomaly(:),'omitnan')) ./ std(dummy_T2M_anomaly(:),'omitnan') ; 
% % dummy_dSM_dt_shallow_standard = (dummy_dSM_dt_shallow - mean(dummy_dSM_dt_shallow(:),'omitnan')) ./ std(dummy_dSM_dt_shallow(:),'omitnan') ; 
% 
% 
%  datatable_test = table(dummy_f_E_anomaly(:), dummy_T2M_anomaly(:), dSM_dt_shallow_anomaly(:),...
%      'VariableNames',{'NDVI_anomaly', 'T2M_anomaly', 'dSM_dt'}) ; 
% 
% % include SM state as well and standardize? 
% dSM_dt_model = fitlm(datatable_test,'dSM_dt ~ 1 + T2M_anomaly + NDVI_anomaly') ; 
% 
% 
% 
%     if size(dummy_dSM_dt_shallow,1) < 300
%         continue
%     end
% 
% 
%     % get 2.5 degree rowcol from name
%     name_files = char(DD_dSM_dt_shallow_files(spatial)) ; 
%     % grep row and col
%     pat = digitsPattern;
%     row_col_name = extract(name_files,pat) ; 
%     row_2_5 = str2double(row_col_name{1}) ; 
%     col_2_5 = str2double(row_col_name{2}) ;   
% 
%     f_E_T2M_corrcoef_vec(row_2_5,col_2_5) = corr_f_E_NDVI_anomalies(1,2) ; 
% 
% 
% 
%     % % only analyze cover or extraction effect
%     % if DD_dSM_dt_NDVI_gradient_median_05(row_2_5,col_2_5)  < 0 
%     %     continue
%     % end
% 
%     % only analyze spatial subset
%     if  ~(col_2_5 > xs(1) && col_2_5 < xs(2)  && row_2_5 > ys(1) && row_2_5 < ys(2))
%         continue
%     end
% 
% 
% 
% 
%     median_spatial_loss = median(dummy_dSM_dt_shallow,1,'omitnan') ; 
%     dSM_dt_shallow_anomaly = dummy_dSM_dt_shallow - median_spatial_loss ;  
%     f_e_extract = dummy_f_E_anomaly;  
%     % dSM_dt_deep_extract = dummy_DD_dSM_dt_deep(:,index_sm) ;     
% 
% 
%     for f_E_bin = 1:49
% 
%             f_E_bin_s = f_E_bins(f_E_bin) ;
%             f_E_bin_e = f_E_bins(f_E_bin+1) ;
% 
%             extract_row = any(f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e,2) ; 
% 
% 
%             dSM_dt_shallow_extract_2 = dSM_dt_shallow_anomaly(extract_row,:) ; 
%             f_e_extract_2 = f_e_extract(extract_row,:) ; 
% 
%             dSM_dt_shallow_extract_2(~(f_e_extract_2 > f_E_bin_s & f_e_extract_2 < f_E_bin_e)) = NaN ; 
% 
% 
%     for j = 1:size(dSM_dt_shallow_extract_2,2)
% 
%         if all(isnan( dSM_dt_shallow_extract_2(:,j)))
%             continue
%         end
% 
%         for i =1:size(dSM_dt_shallow_extract_2,1)
% 
%         if (isnan( dSM_dt_shallow_extract_2(i,j)))
%             continue
%         end
% 
%     dSM_dt_NDVI_binning(f_E_bin,j,counter(f_E_bin,j)) = dSM_dt_shallow_extract_2(i,j); 
% 
%     counter(f_E_bin,j) = counter(f_E_bin,j) + 1 ;
% 
% 
%     if counter(f_E_bin,j) == 150000
%     return
%     end
% 
%         end % i    
%     end  % j 
% 
% 
%     end
% spatial
% end
% 
% 
% 
% % visualize global phase space grow season Zeppe_model
% 
% sampling_2D_bin = sum(~isnan(dSM_dt_NDVI_binning),3) ; 
% xmap = median(dSM_dt_NDVI_binning,3,'omitnan');
% % xmap(sampling_2D_bin < 1000) = NaN ; 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dSM_dt_NDVI_anomaly_Zeppe_05_2D_sampling = sum(~isnan(dSM_dt_NDVI_binning),3) ; 
% dSM_dt_NDVI_anomaly_Zeppe_05_2D = median(dSM_dt_NDVI_binning,3,'omitnan');
% save('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_NDVI_anomaly_Zeppe_05_2D_sampling','dSM_dt_NDVI_anomaly_Zeppe_05_2D_sampling') ; 
% save('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_NDVI_anomaly_Zeppe_05_2D','dSM_dt_NDVI_anomaly_Zeppe_05_2D') ; 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% cd('E:\Daten Baur\Matlab code\redblue')
% redblue_color = redblue(100) ; 
% bam_color = crameri('bam') ;
% load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')
% 
% 
% Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
% h1 = pcolor(sminterp(1:end) - diff(sminterp(1:2))/2,   f_E_bins(1:end) - diff(f_E_bins(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
% clim([-5e-3 5e-3]) %
% % clim([-6e-4 0.0]) %
%  ylim([-0.2 0.2])
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
% imagesc(f_E_T2M_corrcoef_vec) ; colormap(redblue_color) ; clim([-1 1]) ;
% 
% 
% 
% 






%% dostandard 2D phase space into SM and NDVI anomalies but for every year. Then calc trend through the surface






clear

cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')




% cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\spatial_2_5\')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_06_mm_d\spatial_2_5')

cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5')

DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
DD_T2M_anomaly_files = string(ls('*T2M_ERA_anomaly*')) ; 
DD_f_E_anomaly_files = string(ls('*NDVI_anomaly*')) ; 
DD_timev_files = string(ls('*timev*')) ; 


sminterp = linspace(0,0.6,60) ; 
f_E_bins = linspace(-0.2,0.2,41) ; 







% load to check if pos or negative
load('F:\projects\SM_long_term_DDs\data_for_figures\DD_dSM_dt_NDVI_gradient_median_05.mat')
load('F:\projects\SM_long_term_DDs\data_for_figures\DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
DD_dSM_dt_NDVI_gradient_median_05(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% DD_dSM_dt_NDVI_gradient_median_05(DD_dSM_dt_NDVI_gradient_median_05 < 0) = NaN ; 

load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\DD_dSM_dt_NDVI_gradient_zeppe_ERA')
DD_dSM_dt_NDVI_gradient_zeppe_ERA(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 

% imagesc(DD_dSM_dt_NDVI_gradient_median_05) ; colormap(redblue_color) ; clim([-1e-3 1e-3])
% [xs, ys] = getpts() ;  xs  = round(xs) ;   ys  = round(ys) ; 
% close


xs = [1 144] ; ys = [1 72] ; 



years_vec = 1990:2021 ; 
for i = 1:length(years_vec)
    datetime_dummy_years(i) = datetime(strcat('01-Jan-',num2str(years_vec(i)))) ;
end



Zeppe_datetime = datetime('01-Jan-1990'):days(1):datetime('31-Dec-2014') ;




for years_index = 1:24



year_vector_start = datetime_dummy_years(years_index) ; 
year_vector_end = datetime_dummy_years(years_index+1) ; 


dSM_dt_NDVI_binning = NaN(40,40,40000) ; 
counter =  NaN(40,40) ; 
counter(:,:) = 1 ; 




%   spatial = 23 ;
for spatial = 1:length(DD_dSM_dt_shallow_files)

    % cd('F:\projects\SM_long_term_DDs\Zeppe_model_outputs_05\spatial_2_5\')
    % cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5\')
    cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')


    dummy_dSM_dt_shallow = load(DD_dSM_dt_shallow_files(spatial)) ; 
    dummy_f_E_anomaly = load(DD_f_E_anomaly_files(spatial)) ;     
    dummy_timev = load(DD_timev_files(spatial)) ;  

    dummy_dSM_dt_shallow = dummy_dSM_dt_shallow.DD_dSM_dt_shallow_interpsm_spatial ;
    dummy_f_E_anomaly = dummy_f_E_anomaly.DD_NDVI_anomaly_interpsm_spatial ;
    dummy_timev = dummy_timev.DD_timev_spatial ; 




    if size(dummy_dSM_dt_shallow,1) < 300
        continue
    end

    % get rid of SM higher than 0.4 so so it is comparable with SMAP
    dummy_dSM_dt_shallow(:,41:end) = [] ; 
    dummy_f_E_anomaly(:,41:end)  = [] ; 


    % get 2.5 degree rowcol from name
    name_files = char(DD_dSM_dt_shallow_files(spatial)) ; 
    % grep row and col
    pat = digitsPattern;
    row_col_name = extract(name_files,pat) ; 
    row_2_5 = str2double(row_col_name{1}) ; 
    col_2_5 = str2double(row_col_name{2}) ;   



    % % only analyze cover or extraction effect
    if DD_dSM_dt_NDVI_gradient_zeppe_ERA(row_2_5,col_2_5)  < 0 
        continue
    end

    % only analyze spatial subset
    if  ~(col_2_5 > xs(1) && col_2_5 < xs(2)  && row_2_5 > ys(1) && row_2_5 < ys(2))
        continue
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % temporal condition to 1 year
    % time_start = find(ESA_CCI_datetime == datetime('01-Jan-2003')) ; 
    % time_end = find(ESA_CCI_datetime == datetime('01-Jan-2004')) ; 
    time_start = find(Zeppe_datetime == year_vector_start) ; 
    time_end = find(Zeppe_datetime == year_vector_end) ; 



    mask = dummy_timev(:,2) < time_start |  dummy_timev(:,1) > time_end ; % 1-Jan-1990
    dummy_timev(mask) = NaN ; 
    dummy_dSM_dt_shallow(mask,:) = NaN ; 
    dummy_f_E_anomaly(mask,:) = NaN ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    median_spatial_loss = median(dummy_dSM_dt_shallow,1,'omitnan') ; 
    dSM_dt_shallow_anomaly = dummy_dSM_dt_shallow - median_spatial_loss ;  
    f_e_extract = dummy_f_E_anomaly;  
    % dSM_dt_deep_extract = dummy_DD_dSM_dt_deep(:,index_sm) ;     


    for f_E_bin = 1:40

            f_E_bin_s = f_E_bins(f_E_bin) ;
            f_E_bin_e = f_E_bins(f_E_bin+1) ;

            extract_row = any(f_e_extract > f_E_bin_s & f_e_extract < f_E_bin_e,2) ; 


            dSM_dt_shallow_extract_2 = dSM_dt_shallow_anomaly(extract_row,:) ; 
            f_e_extract_2 = f_e_extract(extract_row,:) ; 

            dSM_dt_shallow_extract_2(~(f_e_extract_2 > f_E_bin_s & f_e_extract_2 < f_E_bin_e)) = NaN ; 


    for j = 1:size(dSM_dt_shallow_extract_2,2)
        
        if all(isnan( dSM_dt_shallow_extract_2(:,j)))
            continue
        end

        for i =1:size(dSM_dt_shallow_extract_2,1)

        if (isnan( dSM_dt_shallow_extract_2(i,j)))
            continue
        end
    
    dSM_dt_NDVI_binning(f_E_bin,j,counter(f_E_bin,j)) = dSM_dt_shallow_extract_2(i,j); 

    counter(f_E_bin,j) = counter(f_E_bin,j) + 1 ;


    if counter(f_E_bin,j) == 250000
    return
    end

        end % i    
    end  % j 


    end
spatial
end


% do binning and remove low sampling

dSM_dt_NDVI_binning(dSM_dt_NDVI_binning == 0) = NaN ; 
samples_3D_count = sum(~isnan(dSM_dt_NDVI_binning),3,'omitnan') ;  
xmap = median(dSM_dt_NDVI_binning,3,'omitnan') ; 


save(strcat('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\xmap_annual\xmap_annual_',num2str(year(year_vector_start))),'xmap','samples_3D_count')




years_index
end






cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')

cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_07_mm_d_low_SM\xmap_annual')

xmap(samples_3D_count < 10) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp(1:40) - diff(sminterp(1:2))/2,   f_E_bins(1:end) - diff(f_E_bins(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.5 0.5]) %
% clim([-6e-4 0.0]) %
 ylim([-0.2 0.2])
 xlim([0 0.5])
 colormap(redblue_color) %
% colormap parula
cbr = colorbar ;
cbr.Label.String = "SM loss anomaly [mm/day]";
xlabel('SM [m³/m³]')
ylabel('f_E anomaly [-]')
title('Europe')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% fontsize(16,'points')


% saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\SM_loss_shallow_2D_SM_fE_anomalies_Europe_100samp','svg')
% close 


figure
h1 = pcolor(sminterp(1:40) - diff(sminterp(1:2))/2,   f_E_bins(1:end) - diff(f_E_bins(1:2))/2, (samples_3D_count)) ;
clim([0 400])





clear

cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\xmap_annual')

filenames = string(ls('*xmap*')) ; 


xmap_array = NaN(40,40,24) ; 

for i = 1:24

    dummy = matfile(filenames(i)) ; 
    dummy_samples = dummy.samples_3D_count ;     
    dummy = dummy.xmap ;
    dummy(dummy_samples < 10) = NaN ; 
    xmap_array(:,:,i) = dummy ; 

i
end


Xmap_space_slope = NaN(40,40) ; 
Xmap_space_intercept= NaN(40,40) ; 
Xmap_space_slope_p = NaN(40,40) ; 

for r = 1:40
    for c = 1:40

        dummy_ts = squeeze(xmap_array(r,c,:)) ; 
        % plot(dummy_ts)
        if (all(~isnan(dummy_ts)))

              
              cd('E:\Daten Baur\Matlab code')
              [m_cur b_cur] = TheilSen([(1:length(dummy_ts))' , dummy_ts]) ; 
              Xmap_space_slope(r,c) = m_cur ;    
              Xmap_space_intercept(r,c) = b_cur ;               
              [H,p_value] = Mann_Kendall(dummy_ts,0.05)  ;   
              Xmap_space_slope_p(r,c) = p_value ; 


        else
            continue
        end


    end
    r
end




Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp(1:40),f_E_bins(1:end) - diff(f_E_bins(1:2))/2,Xmap_space_slope) ;
hold on
set(h1,'LineStyle','none')
shading flat
clim([-0.025 0.025]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.6])
colormap(bam_color) %
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat anomaly trend [m³/m³/days/year]";
xticks([0:0.2:0.6])
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
title('1990-2014 Theil Sen slope')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% plot significance
 [rowfind colfind] = find(Xmap_space_slope_p < 0.05) ; 
 points1 = plot(sminterp(colfind)+ 0.0100/2,f_E_bins(rowfind+1) ,'.k','MarkerSize',5) ; 

% print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSM_dt_anomaly_NDVI_tile_trend','-dbmp','-r150')
% close


SM_loss_mm_d_Zeppe_trend = Xmap_space_slope ;
SM_loss_mm_d_Zeppe_trend_p =  Xmap_space_slope_p ; 

save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\SM_loss_mm_d_Zeppe_trend','SM_loss_mm_d_Zeppe_trend')
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\SM_loss_mm_d_Zeppe_trend_p','SM_loss_mm_d_Zeppe_trend_p')










%% 2.5° degree map of SM loss trends. Only do binning in SM axis

clear

cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')


cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5')

DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
DD_T2M_anomaly_files = string(ls('*T2M_ERA_anomaly*')) ; 
DD_f_E_anomaly_files = string(ls('*NDVI_anomaly*')) ; 
DD_timev_files = string(ls('*timev*')) ; 


sminterp = linspace(0,0.6,60) ; 
f_E_bins = linspace(-0.2,0.2,40) ; 




% load to check if pos or negative
load('F:\projects\SM_long_term_DDs\data_for_figures\DD_dSM_dt_NDVI_gradient_median_05.mat')
load('F:\projects\SM_long_term_DDs\data_for_figures\DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
DD_dSM_dt_NDVI_gradient_median_05(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% DD_dSM_dt_NDVI_gradient_median_05(DD_dSM_dt_NDVI_gradient_median_05 < 0) = NaN ; 

load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\DD_dSM_dt_NDVI_gradient_zeppe_ERA')
DD_dSM_dt_NDVI_gradient_zeppe_ERA(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 



lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 

xs = [1 144] ; ys = [1 72] ; 


years_vec = 1990:2021 ; 
for i = 1:length(years_vec)
    datetime_dummy_years(i) = datetime(strcat('01-Jan-',num2str(years_vec(i)))) ;
end

Zeppe_datetime = datetime('01-Jan-1990'):days(1):datetime('31-Dec-2014') ;


dSM_dt_TS_MK_trend_p = NaN(72,144) ; 
dSM_dt_mean_all= NaN(72,144,24) ; 
dSM_dt_TS_slope_all = NaN(72,144) ; 

% spatial = 2050
for spatial = 1:size(DD_T2M_anomaly_files,1)

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')
cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5')

dummy_dSM_dt = matfile(DD_dSM_dt_shallow_files(spatial)) ; 
dummy_dSM_dt = dummy_dSM_dt.DD_dSM_dt_shallow_interpsm_spatial     ; 
% dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% exclude dainage
dummy_dSM_dt(:,40:end) = NaN ; 
% cut off clumns?
dummy_dSM_dt(:,41:end) = [] ; 

    
dummy_t_start_end = matfile(DD_timev_files(spatial)) ; 
dummy_t_start_end = dummy_t_start_end.DD_timev_spatial ;    




% get 2.5 degree rowcol from name
name = char(DD_dSM_dt_shallow_files(spatial)) ; 
% grep row and col
pat = digitsPattern;
row_col_name = extract(name,pat) ; 
row_2_5 = str2double(row_col_name{1}) ; 
col_2_5 = str2double(row_col_name{2}) ;   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new condition to check if too northern
    lat_2_5_dummy = lats_2_5(row_2_5,col_2_5) ; 
    lon_2_5_dummy = lons_2_5(row_2_5,col_2_5) ; 

    if lat_2_5_dummy > 70
     continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   



 % get average annual and then trend   
 dSM_dt_year_dummy = NaN(24,40) ; 

% year_vector = 1990:2:2019 ; 
 year_vector = 1990:2014 ; 
 
for i =  1:length(year_vector)-1% 1990:2:2019
    
    year_cur = year_vector(i) ;
    years_overlap = find(year(Zeppe_datetime) == year_cur) ; 
   
     select_year = dummy_t_start_end(:,1) > min(years_overlap) & dummy_t_start_end(:,2) < max(years_overlap) ; 

     dummy_select = dummy_dSM_dt(select_year,:) ; 
     dummy_select(:,sum(~isnan(dummy_select)) < 10) = NaN ; 

    %imagesc(dummy_select)
    dSM_dt_year_dummy(i,:) = median( dummy_select,1,'omitnan') ;   

    
end



% imagesc(dSM_dt_year_dummy)
dSM_dt_year_dummy(dSM_dt_year_dummy == 0) = NaN ; 
length_SM_loss_function = sum(~isnan(dSM_dt_year_dummy),2) ; 
all_valid_SM_loss_function = all(~isnan(dSM_dt_year_dummy),1) ; 


% check new condition, remove SM loss functions shortert than 10. Arbitrary
% but should get rid of mini loss functions that are annomalous.
% dSM_dt_year_dummy(length_SM_loss_function < 10,:) = NaN ; 


dSM_dt_year_dummy = median(dSM_dt_year_dummy,2,'omitnan') ; 


cd('E:\Daten Baur\Matlab code')
% correlation between loss rates and NDVI could do correlation across SM
% conditions and then mean or distirbution?? 
  [m_cur, b_cur] = TheilSen([(1:24)' , dSM_dt_year_dummy]) ; 



mask_2 = isnan(dSM_dt_year_dummy) ; 

dSM_dt_year_dummy_short = dSM_dt_year_dummy ; 
dSM_dt_year_dummy_short(mask_2) = [] ; 


 if all(isnan(dSM_dt_year_dummy))
     continue
 end



   [H,p_value] = Mann_Kendall(dSM_dt_year_dummy,0.05)  ;   
   dSM_dt_TS_MK_trend_p(row_2_5,col_2_5) = p_value ; 

    dSM_dt_mean_all(row_2_5,col_2_5,:) = dSM_dt_year_dummy ; 
    dSM_dt_TS_slope_all(row_2_5,col_2_5,:) = m_cur ; 



   spatial
   
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dSM_dt_slope_1990_2014_zeppe = dSM_dt_TS_slope_all ; 
dSM_dt_MK_trend_p_1990_2014_zeppe = dSM_dt_TS_MK_trend_p ; 

save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_slope_1990_2014_zeppe','dSM_dt_slope_1990_2014_zeppe') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_MK_trend_p_1990_2014_zeppe','dSM_dt_MK_trend_p_1990_2014_zeppe') ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
load('E:\Daten Baur\Matlab files\means_über_zeitreihe\Coastlines.mat')


lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 



xmap = dSM_dt_TS_slope_all ;
figure('units','centimeters','position',[10 3 39 20])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
%colormap(redblue_color) 
colormap(bam_color)
clim([-0.05  0.05])
  % clim([-1  1]) 
hcb=colorbar;
ylabel(hcb,'dSM/dt trend [m³/m³/day/yr]','FontSize',17)
% title('non growing season','FontSize',17)
xlabel('longitude')
ylabel('latitude')
%title('ESA CCI SM 1990-2020 trend')
plot(CoastlineLon, CoastlineLat,'Color','k');
hold on
[rowfind colfind] = find(dSM_dt_TS_MK_trend_p < 0.05) ; 
 crosses1 = plot(lons_2_5(1,colfind),lats_2_5(rowfind,1),'.k','MarkerSize',7) ; 
set(gca,'FontSize',17)





%% do a trend of SM loss function just over SM dimension


clear


cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')
% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5')

DD_dSM_dt_shallow_files = string(ls('*shallow*')) ; 
DD_timev_files = string(ls('*timev*')) ; 


sminterp = linspace(0,0.6,60) ; 
f_E_bins = linspace(-0.2,0.2,40) ; 


cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d')
load('SM_loss_function_array_ESACCI.mat')
load('SM_loss_function_array_samples_ESACCI.mat')
load('SM_loss_function_array_Zeppe.mat')
load('SM_loss_function_array_Zeppe_samples.mat')


figure
plot(SM_loss_function_array_Zeppe(:,24))
hold on
plot(SM_loss_function_array_ESACCI(:,24))


figure
plot(SM_loss_function_array_Zeppe(12,:))
hold on
plot(SM_loss_function_array_ESACCI(12,:))


figure
plot(mean(SM_loss_function_array_Zeppe,2,'omitnan'))
hold on
plot(mean(SM_loss_function_array_ESACCI,2,'omitnan'))




figure
plot(SM_loss_function_array_Zeppe(:,22))
hold on
plot(SM_loss_function_array_ESACCI(:,22))








load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
DD_dSM_dt_NDVI_gradient_median_05(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% DD_dSM_dt_NDVI_gradient_median_05(DD_dSM_dt_NDVI_gradient_median_05 < 0) = NaN ; 
load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\DD_dSM_dt_NDVI_gradient_zeppe_ERA')
DD_dSM_dt_NDVI_gradient_zeppe_ERA(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 

lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 

years_vec = 1990:2021 ; 
for i = 1:length(years_vec)
    datetime_dummy_years(i) = datetime(strcat('01-Jan-',num2str(years_vec(i)))) ;
end

Zeppe_datetime = datetime('01-Jan-1990'):days(1):datetime('31-Dec-2014') ;
load('F:\CRUJRA\processed\datetime_1990_2014.mat')


SM_loss_function_trend = NaN (1,40) ;
SM_loss_function_trend_p = NaN (1,40) ;

% spatial = 2050

counter = 1 ; 
SM_loss_function_array = NaN(24,40) ; 
SM_loss_function_array_samples = NaN(24,40) ; 

for years_index = 1:24


year_vector_start = datetime_dummy_years(years_index) ; 
year_vector_end = datetime_dummy_years(years_index+1) ; 

year_vector_start_index = find(Zeppe_datetime == year_vector_start) ; 
year_vector_end_index = find(Zeppe_datetime == year_vector_end) ; 

dSM_dt_dummy_year_array = NaN(70000,40) ; 
counter = 1 ; 
% spatial = 345
for spatial = 1:size(DD_dSM_dt_shallow_files,1)

% cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\spatial_2_5_mm_d')
 cd('F:\projects\SM_long_term_DDs\Zeppe_model_output_tests\spatial_2_5')

dummy_dSM_dt = matfile(DD_dSM_dt_shallow_files(spatial)) ; 
dummy_dSM_dt = dummy_dSM_dt.DD_dSM_dt_shallow_interpsm_spatial     ; 
% dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% exclude dainage
dummy_dSM_dt(:,41:end) = [] ; 

% get 2.5 degree rowcol from name
name = char(DD_dSM_dt_shallow_files(spatial)) ; 
% grep row and col
pat = digitsPattern;
row_col_name = extract(name,pat) ; 
row_2_5 = str2double(row_col_name{1}) ; 
col_2_5 = str2double(row_col_name{2}) ;   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new condition to check if too northern
    lat_2_5_dummy = lats_2_5(row_2_5,col_2_5) ; 
    lon_2_5_dummy = lons_2_5(row_2_5,col_2_5) ; 

    if lat_2_5_dummy > 70
     continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
dummy_t_start_end = matfile(DD_timev_files(spatial)) ; 
dummy_t_start_end = dummy_t_start_end.DD_timev_spatial ;    

% size(dummy_t_start_end)
% size(dummy_dSM_dt)


select_year = dummy_t_start_end(:,1) > year_vector_start_index & dummy_t_start_end(:,2) < year_vector_end_index ; 

dummy_select = dummy_dSM_dt(select_year,:) ; 


dSM_dt_dummy_year_array(counter:counter+size(dummy_select,1)-1,:) = dummy_select ; 
       counter  = counter + length(dummy_select) ; 

 spatial

end


dSM_dt_dummy_year_array_median = median(dSM_dt_dummy_year_array,1,'omitnan') ;
dSM_dt_dummy_year_array_median_samples = sum(~isnan(dSM_dt_dummy_year_array),1) ; 


SM_loss_function_array(years_index,:) = dSM_dt_dummy_year_array_median ;
SM_loss_function_array_samples(years_index,:) = dSM_dt_dummy_year_array_median_samples ;

years_index
end



SM_loss_function_array_Zeppe = SM_loss_function_array ;
SM_loss_function_array_Zeppe_samples = SM_loss_function_array_samples ;




save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\SM_loss_function_array_Zeppe_nu20test','SM_loss_function_array_Zeppe') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\SM_loss_function_array_Zeppe_nu20test_samples','SM_loss_function_array_Zeppe_samples') ; 








SM_loss_trend = NaN(1,60) ; 
SM_loss_trend_p = NaN(1,60) ; 
cd('E:\Daten Baur\Matlab code')

for i = 1:40

              [m_cur b_cur] = TheilSen([(1:24)' , SM_loss_function_array(:,i)]) ; 
              SM_loss_trend(:,i) = m_cur ;    
              % Xmap_space_intercept(r,c) = b_cur ;               
              [H,p_value] = Mann_Kendall(SM_loss_function_array(:,i),0.05)  ;   
              SM_loss_trend_p(:,i) = p_value ; 



end

Global_SM_loss_trend    = SM_loss_trend   ; 
Global_SM_loss_trend_p  = SM_loss_trend_p ;

save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\Global_SM_loss_trend','Global_SM_loss_trend') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\Global_SM_loss_trend_p','Global_SM_loss_trend_p') ; 







