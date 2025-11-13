%% MJB 12-AUg-2023 .. bin into LAI percentiles month and make a tile plot. check if higher LAI slows down
% drydowns below threshold and speeds them up above. Explore roles on
% veggetaiton
clear

%% for plotting


cd('E:\Daten Baur\Matlab code')
col_L = hex2rgb('808080') ; 
col_L_light = hex2rgb('CDCDCD') ; 

col_C = hex2rgb('3399FF') ;
col_C_light = hex2rgb('80E6FF') ; 

col_X = hex2rgb('FF3333') ;
col_X_light = hex2rgb('FF8080') ; 


lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 

cd('E:\Daten Baur\Matlab files\means_über_zeitreihe')
load('Coastlines.mat')
cd('E:\Daten Baur\Matlab code\redblue')
redblue_color = redblue(100) ; 
bam_color = crameri('bam') ;
tokyo_color = crameri('tokyo') ;
imola_color = crameri('imola') ;
cork_color = crameri('cork') ;
batlow_color = crameri('batlow') ;
sminterp = (0.01:0.01:0.6)' ; 
sminterp_30 = (0.01:0.02:0.6)' ; 
sminterp_40 = linspace(0.01,0.6,40)' ; 

% NDVI_binning = linspace(-0.15,0.15,41) ; 
% dSM_dt_binning = linspace(-0.08,0.08,41) ; 




%% here just use all data
clear




cd('F:\ESA_CCI\analysis_datasets\spatial_02')
% dd data
dSM_dt_files = string(ls('*dSM_dt_interpsm*')) ;
rowcol_files = string(ls('*rowcol_year_dummy*')) ;
t_start_end_files = string(ls('*t_start*')) ;

% monthly LAI will be only until 2014
cd('F:\ESA_CCI\analysis_datasets\NDVI_sminterp')
NDVI_files_all = string(ls('*VIP_NDVI_daily*')) ;
NDVI_files_121 = string(ls('*VIP_NDVI_daily_sminterp_121_*')) ;
NDVI_files_5y = string(ls('*VIP_NDVI_daily_sminterp_anomaly*')) ;

NDVI_files_all = strtrim(NDVI_files_all) ;
NDVI_files_121 = strtrim(NDVI_files_121) ;
NDVI_files_5y = strtrim(NDVI_files_5y) ;



% T2M for correction
cd('F:\ESA_CCI\analysis_datasets\T2M_sminterp')
T2M_files = string(ls('*T2M*')) ;
T2M_121_files = string(ls('*T2M_daily_sminterp_121_*')) ;
[Lia Locb] = ismember(T2M_121_files,T2M_files) ; 
T2M_sample_files = T2M_files ; 
T2M_sample_files(Locb) = [] ; 



[Lia Locb] = ismember(NDVI_files_121, NDVI_files_all,'rows') ; 
NDVI_files = NDVI_files_all ; 
NDVI_files(Locb) = [] ; 



ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 

NDVI_binning = linspace(-0.20,0.20,41) ; 
dSM_dt_binning = linspace(-0.08,0.08,41) ; 

sxity_to_thirty = 0:2:60 ; sxity_to_thirty(1) = 1 ; 

counter = NaN(40,40) ; 
counter(:,:) = 1 ; 


% binning based on NDVI
%dSM_dt_NDVI_binning = NaN(40,40,1000000) ; 
dSM_dt_NDVI_binning = NaN(40,40,2000000) ; 




% spatial = 2450
% spatial = 2550
cd('F:\ESA_CCI\analysis_datasets')

for spatial = 1:length(dSM_dt_files)


% get drydown data

dummy_dSM_dt = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',dSM_dt_files(spatial))) ; 
NDVI_dummy =  matfile(strcat('F:\ESA_CCI\analysis_datasets\NDVI_sminterp\',NDVI_files_5y(spatial))) ; 
dummy_T2M =  matfile(strcat('F:\ESA_CCI\analysis_datasets\T2M_sminterp\',T2M_sample_files(spatial))) ; 


% get 2.5 degree rowcol from name
name_files = char(dSM_dt_files(spatial)) ; 
% grep row and col
pat = digitsPattern;
row_col_name = extract(name_files,pat) ; 
row_2_5 = str2double(row_col_name{1}) ; 
col_2_5 = str2double(row_col_name{2}) ;     


dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy   ; 
dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% new condition removing very short DDs?
dummy_dSM_dt(sum(~isnan(dummy_dSM_dt),2) < 3 , :) = NaN ; 
% kick out all drainage to save memory. reduced to 40 length from 60
dummy_dSM_dt(:,41:60) = [] ;


 % dummy_NDVI = NDVI_dummy.dummy_NDVI_monthly ; 
 dummy_NDVI = NDVI_dummy.NDVI_spatial_sminterp ; 
 dummy_NDVI(:,41:60) = [] ;

 dummy_T2M = dummy_T2M.T2M_spatial_sminterp ; 
 dummy_T2M(:,41:60) = [] ;




% exclude dainage
%dummy_dSM_dt(:,40:end) = NaN ; 


dummy_rowcol = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',rowcol_files(spatial))) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',t_start_end_files(spatial))) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 

% pre filter for time if needed
mask = dummy_t_mean < 4080 ; % 1-Jan-1990
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_NDVI(mask,:) = NaN ; 
dummy_T2M(mask,:) = NaN ; 


% MJB could make faster by not indexing for pos/neg T2M but setting dSM/dt
% to NAN when T2M pos or neg her... maybe use for better performance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporal condition to 1 year
% time_start = find(ESA_CCI_datetime == datetime('01-Jan-2003')) ; 
% time_end = find(ESA_CCI_datetime == datetime('01-Jan-2004')) ; 
% % time_start = find(ESA_CCI_datetime == year_vector_start) ; 
% % time_end = find(ESA_CCI_datetime == year_vector_end) ; 



% mask = dummy_t_mean < time_start |  dummy_t_mean > time_end ; % 1-Jan-1990
% dummy_t_mean(mask) = NaN ; 
% dummy_dSM_dt(mask,:) = NaN ; 
% dummy_rowcol(mask,:) = NaN ; 
% dummy_NDVI(mask,:) = NaN ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(any(~isnan(dummy_dSM_dt)))
% return
% end



% binning based on NDVI
% dSM_dt_NDVI_binning = NaN(30,60) ; 
%NDVI_mean = mean(dummy_NDVI,'omitnan') ; 
NDVI_deviations = dummy_NDVI  ; 
T2M_deviations = dummy_T2M ; 
% get mean but for each SM bin
dummy_dSM_dt_mean = mean(dummy_dSM_dt,1,'omitnan') ; 

for bin = 1:40

    cur_bin_min = NDVI_binning(bin) ;
    cur_bin_max = NDVI_binning(bin+1) ;   
    % standard bining 
    % extract_index = NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max ;
    % extract_index_row =  any(NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max,2) ; 

    % add condition of positive/negative T2M anomaly
    extract_index = NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max & T2M_deviations > 0;
    extract_index_row =  any(NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max & T2M_deviations > 0,2) ; 



    % get dSM/dt anomalies
    extract_index = extract_index(extract_index_row,:) ; 
    dSM_dt_extract = dummy_dSM_dt(extract_index_row,:) - dummy_dSM_dt_mean ; 
    dSM_dt_extract(~extract_index) = NaN ; 

     % now do propper 2D binning with i and j to use memory efficiently.
     % Binning array is filled with 3rd dim index specific for row and col
    for j = 1:size(dSM_dt_extract,2)

        if all(isnan( dSM_dt_extract(:,j)))
            continue
        end

        for i =1:size(dSM_dt_extract,1)

        if (isnan( dSM_dt_extract(i,j)))
            continue
        end
    
    dSM_dt_NDVI_binning(bin,j,counter(bin,j)) = dSM_dt_extract(i,j); 

    counter(bin,j) = counter(bin,j) + 1 ;


    if counter(bin,j) == 2000000
    return
    end

        end % i
    
    end  % j 


end

 spatial


end




dSM_dt_NDVI_binning(dSM_dt_NDVI_binning == 0) = NaN ; 
samples_3D_count = sum(~isnan(dSM_dt_NDVI_binning),3,'omitnan') ;  



xmap = mean(dSM_dt_NDVI_binning,3,'omitnan') ; 
 xmap(samples_3D_count < 1000) = NaN ; 
% xmap(samples_3D_count < 100) = NaN ; 

% imagesc(samples_3D_count)


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.005 0.005]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(redblue_color) %
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat anomaly [m³/m³/day]";
cbr.Label.FontSize = 17 ; 
xticks([0:0.2:0.6])
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
title('T2M anomalies are > 0')
pbaspect([0.8 0.8 0.8])
axes = gca ; 
axes.Units = 'centimeters' ; 
axes_diff_Position = get(axes, 'Position');
% calc positions
%  3.6987    2.3080   22.0499   17.0999   are the pos 
 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',6,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',6,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow1,'Position',[3.6987+22.0499+1 2.3080+(17.0999/2)+1   0  (17.0999/2)-1]) ; 
 set(arrow2,'Position',[3.6987+22.0499+1 2.3080+(17.0999/2)-1   0 -(17.0999/2)+1]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',17,'Units','centimeters');
set(textbox1,'Position',[3.6987+22.0499+2, 2.3080+(17.0999/2)+1,    -(17.0999/2)+1, 0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',17,'Units','centimeters');
set(textbox2,'Position',[3.6987+22.0499+2, 2.3080+(17.0999/2)-1,    (17.0999/2)-1, 0]) ; 
set(gca,'Box','on');
set(gca,'FontSize',17)



print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSM_dt_anomaly_NDVI_T2M_pos_only_all','-dbmp','-r150')
close




%% do annual plot 





clear

cd('F:\ESA_CCI\analysis_datasets\spatial_02')
% dd data
dSM_dt_files = string(ls('*dSM_dt_interpsm*')) ;
rowcol_files = string(ls('*rowcol_year_dummy*')) ;
t_start_end_files = string(ls('*t_start*')) ;

% monthly LAI will be only until 2014
cd('F:\ESA_CCI\analysis_datasets\NDVI_sminterp')
NDVI_files_all = string(ls('*VIP_NDVI_daily*')) ;
NDVI_files_121 = string(ls('*VIP_NDVI_daily_sminterp_121_*')) ;
NDVI_files_5y = string(ls('*VIP_NDVI_daily_sminterp_anomaly*')) ;

NDVI_files_all = strtrim(NDVI_files_all) ;
NDVI_files_121 = strtrim(NDVI_files_121) ;
NDVI_files_5y = strtrim(NDVI_files_5y) ;

% find spatial locations of pos neg 
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_pos')
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_neg')
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_all')
dSM_dt_mean_NDVI_pos_mask = dSM_dt_mean_NDVI_pos > 0 ; 


ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 

NDVI_binning = linspace(-0.20,0.20,41) ; 
dSM_dt_binning = linspace(-0.08,0.08,41) ; 

sxity_to_thirty = 0:2:60 ; sxity_to_thirty(1) = 1 ; 

counter = NaN(60,1) ; 
counter(:,1) = 1 ; 

% binning based on NDVI
dSM_dt_NDVI_binning = NaN(40,40,100000) ; 







% ESA_CCI_datetime(4080)
% ESA_CCI_datetime(end)
years_vec = 1990:2021 ; 
for i = 1:length(years_vec)
    datetime_dummy_years(i) = datetime(strcat('01-Jan-',num2str(years_vec(i)))) ;
end
datetime_dummy_years(end) = ESA_CCI_datetime(end) ; 




% stopped 10
cd('F:\ESA_CCI\analysis_datasets')
for years_index = 1:length(datetime_dummy_years)-1


counter = NaN(40,40) ; 
counter(:,:) = 1 ; 
% binning based on NDVI
% NDVI_dSM_dt_binning = NaN(60,60,500000) ; 


year_vector_start = datetime_dummy_years(years_index) ; 
year_vector_end = datetime_dummy_years(years_index+1) ; 


% binning based on NDVI
dSM_dt_NDVI_binning = NaN(40,40,400000) ; 





% spatial = 2450
% spatial = 2550

for spatial = 1: length(dSM_dt_files)


% get drydown data
dummy_dSM_dt = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',dSM_dt_files(spatial))) ; 
NDVI_dummy =  matfile(strcat('F:\ESA_CCI\analysis_datasets\NDVI_sminterp\',NDVI_files_5y(spatial))) ; 


% get 2.5 degree rowcol from name
name_files = char(dSM_dt_files(spatial)) ; 
% grep row and col
pat = digitsPattern;
row_col_name = extract(name_files,pat) ; 
row_2_5 = str2double(row_col_name{1}) ; 
col_2_5 = str2double(row_col_name{2}) ;     


% this checks whether we are on cover or extraction effect dominated
if (dSM_dt_mean_NDVI_pos_mask(row_2_5,col_2_5) < 0)
    continue
end



dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy   ; 
dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% new condition removing very short DDs?
% dummy_dSM_dt(sum(~isnan(dummy_dSM_dt),2) < 3 , :) = NaN ; 
% kick out all drainage to save memory. reduced to 40 length from 60
dummy_dSM_dt(:,41:60) = [] ;


 % dummy_NDVI = NDVI_dummy.dummy_NDVI_monthly ; 
 dummy_NDVI = NDVI_dummy.NDVI_spatial_sminterp ; 
 dummy_NDVI(:,41:60) = [] ;


% exclude dainage
%dummy_dSM_dt(:,40:end) = NaN ; 

dummy_rowcol = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',rowcol_files(spatial))) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',t_start_end_files(spatial))) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 

% pre filter for time if needed
mask = dummy_t_mean < 4080 ; % 1-Jan-1990
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_NDVI(mask,:) = NaN ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporal condition to 1 year
% time_start = find(ESA_CCI_datetime == datetime('01-Jan-2003')) ; 
% time_end = find(ESA_CCI_datetime == datetime('01-Jan-2004')) ; 
time_start = find(ESA_CCI_datetime == year_vector_start) ; 
time_end = find(ESA_CCI_datetime == year_vector_end) ; 



mask = dummy_t_mean < time_start |  dummy_t_mean > time_end ; % 1-Jan-1990
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_NDVI(mask,:) = NaN ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(any(~isnan(dummy_dSM_dt)))
% return
% end



% binning based on NDVI
% dSM_dt_NDVI_binning = NaN(30,60) ; 
%NDVI_mean = mean(dummy_NDVI,'omitnan') ; 
NDVI_deviations = dummy_NDVI  ; 
% get mean but for each SM bin
dummy_dSM_dt_mean = mean(dummy_dSM_dt,1,'omitnan') ; 

for bin = 1:40

    cur_bin_min = NDVI_binning(bin) ;
    cur_bin_max = NDVI_binning(bin+1) ;   
    extract_index = NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max ;
    extract_index_row =  any(NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max,2) ; 
    % get dSM/dt anomalies
    extract_index = extract_index(extract_index_row,:) ; 
    dSM_dt_extract = dummy_dSM_dt(extract_index_row,:) - dummy_dSM_dt_mean ; 
    dSM_dt_extract(~extract_index) = NaN ; 

     % now do propper 2D binning with i and j to use memory efficiently.
     % Binning array is filled with 3rd dim index specific for row and col
    for j = 1:size(dSM_dt_extract,2)

        if all(isnan( dSM_dt_extract(:,j)))
            continue
        end

        for i =1:size(dSM_dt_extract,1)

        if (isnan( dSM_dt_extract(i,j)))
            continue
        end
    
    dSM_dt_NDVI_binning(bin,j,counter(bin,j)) = dSM_dt_extract(i,j); 

    counter(bin,j) = counter(bin,j) + 1 ;


    if counter(bin,j) == 1000000
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


xmap = mean(dSM_dt_NDVI_binning,3,'omitnan') ; 
% xmap(samples_3D_count < 100) = NaN ; 
  % figure
  % pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2,  (xmap)) ; colormap(redblue_color) ; clim([-0.003 0.003])  

save(strcat('F:\ESA_CCI\analysis_datasets\xmap_dSm_dt_02\dSM_dt_dev_from_seas_anomaly5y_cover_only_',num2str(year(year_vector_start))),'xmap','samples_3D_count')


year_vector_start






end





%% use xmaps to do trends of dSM/dt NDVI sens ...


%% get annual tileplots for trend analysis
clear

cd('F:\ESA_CCI\analysis_datasets\xmap_dSm_dt_03')

filenames = string(ls('*dSM_dt_dev_from_seas_anomaly5y*')) ; 
filenames = string(ls('*dSM_dt*')) ; 


xmap_array = NaN(40,40,10) ; 

for i = 1:24

    dummy = matfile(filenames(i)) ; 
    dummy_samples = dummy.samples_3D_count ;     
    dummy = dummy.xmap ;
    dummy(dummy_samples < 100) = NaN ; 
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
h1 = pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2,Xmap_space_slope) ;
hold on
set(h1,'LineStyle','none')
shading flat
clim([-0.01 0.01]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.6])
% colormap(redblue_color) %
colormap(bam_color)
cbr = colorbar ;
cbr.Label.String = "SM loss trend [mm/day/year]";
cbr.Ticks = -0.01:0.005:0.01 ;
xticks([0:0.2:0.6]) ;
yticks([-0.25:0.125:0.25])
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
title('1990-2014 Theil Sen slope')
title('')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% plot significance
 [rowfind colfind] = find(Xmap_space_slope_p < 0.05) ; 
 points1 = plot(sminterp(colfind)+ 0.0100/2,NDVI_binning(rowfind+1) ,'.k','MarkerSize',5) ; 
fontsize(18,'points')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures\Panel_XX_SM_loss_trends_phase_space','svg')
close 





Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp(1:40),dSM_dt_binning(2:end) - diff(dSM_dt_binning(1:2))/2,xmap_array(:,:,12)) ;
hold on
set(h1,'LineStyle','none')
shading flat
clim([-0.01 0.01]) %
 ylim([-0.08 0.08])
 xlim([-0.00 0.6])
colormap(redblue_color) %





%% Do binning into SM and abs(NDVI) for just dSM/dt no anomaly
clear



cd('F:\ESA_CCI\analysis_datasets\spatial_02')
% dd data
dSM_dt_files = string(ls('*dSM_dt_interpsm*')) ;
rowcol_files = string(ls('*rowcol_year_dummy*')) ;
t_start_end_files = string(ls('*t_start*')) ;


% monthly LAI will be only until 2014 exchange this for new SE NDVI
cd('F:\ESA_CCI\analysis_datasets\NDVI_sminterp_noanomaly')
NDVI_files_all = string(ls('*VIP_NDVI_daily*')) ;
NDVI_files_all = strtrim(NDVI_files_all) ;



ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 

NDVI_binning = linspace(0,1,41) ; 
counter = NaN(40,40) ; 
counter(:,:) = 1 ; 


% binning based on NDVI
% dSM_dt_NDVI_binning = NaN(60,60,1000000) ; 
dSM_dt_NDVI_binning = NaN(40,40,2000000) ; 




% spatial = 2450
% spatial = 2550
cd('F:\ESA_CCI\analysis_datasets')

for spatial = 1:length(dSM_dt_files)


% get drydown data

dummy_dSM_dt = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',dSM_dt_files(spatial))) ; 
NDVI_dummy =  matfile(strcat('F:\ESA_CCI\analysis_datasets\NDVI_sminterp_noanomaly\',NDVI_files_all(spatial))) ; 

% get 2.5 degree rowcol from name
name_files = char(dSM_dt_files(spatial)) ; 
% grep row and col
pat = digitsPattern;
row_col_name = extract(name_files,pat) ; 
row_2_5 = str2double(row_col_name{1}) ; 
col_2_5 = str2double(row_col_name{2}) ;     


dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy    ; 
dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% new condition removing very short DDs?
dummy_dSM_dt(sum(~isnan(dummy_dSM_dt),2) < 3 , :) = NaN ; 
% kick out all drainage to save memory. reduced to 40 length from 60
 dummy_dSM_dt(:,41:60) = [] ;


 % dummy_NDVI = NDVI_dummy.dummy_NDVI_monthly ; 
 dummy_NDVI = NDVI_dummy.NDVI_spatial_sminterp ; 
 % add this once we use NDVI interpsm
 dummy_NDVI(:,41:60) = [] ;


% exclude dainage
%dummy_dSM_dt(:,40:end) = NaN ; 


dummy_rowcol = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',rowcol_files(spatial))) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',t_start_end_files(spatial))) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 

% pre filter for time if needed
mask = dummy_t_mean < 4080 ; % 1-Jan-1990
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_NDVI(mask,:) = NaN ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporal condition to 1 year
% time_start = find(ESA_CCI_datetime == datetime('01-Jan-2003')) ; 
% time_end = find(ESA_CCI_datetime == datetime('01-Jan-2004')) ; 
% % time_start = find(ESA_CCI_datetime == year_vector_start) ; 
% % time_end = find(ESA_CCI_datetime == year_vector_end) ; 



% mask = dummy_t_mean < time_start |  dummy_t_mean > time_end ; % 1-Jan-1990
% dummy_t_mean(mask) = NaN ; 
% dummy_dSM_dt(mask,:) = NaN ; 
% dummy_rowcol(mask,:) = NaN ; 
% dummy_NDVI(mask,:) = NaN ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(any(~isnan(dummy_dSM_dt)))
% return
% end



% binning based on NDVI
% dSM_dt_NDVI_binning = NaN(30,60) ; 
%NDVI_mean = mean(dummy_NDVI,'omitnan') ; 
NDVI_deviations = dummy_NDVI  ; 
% get mean but for each SM bin
% dummy_dSM_dt_mean = mean(dummy_dSM_dt,1,'omitnan') ; 

for bin = 1:40

    cur_bin_min = NDVI_binning(bin) ;
    cur_bin_max = NDVI_binning(bin+1) ;   
    extract_index = NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max ;
    extract_index_row =  any(NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max,2) ; 
    % get dSM/dt anomalies
    extract_index = extract_index(extract_index_row,:) ; 
    dSM_dt_extract = dummy_dSM_dt(extract_index_row,:) ; 
    dSM_dt_extract(~extract_index) = NaN ; 

     % now do propper 2D binning with i and j to use memory efficiently.
     % Binning array is filled with 3rd dim index specific for row and col
    for j = 1:size(dSM_dt_extract,2)

        if all(isnan( dSM_dt_extract(:,j)))
            continue
        end

        for i =1:size(dSM_dt_extract,1)

        if (isnan( dSM_dt_extract(i,j)))
            continue
        end
    
    dSM_dt_NDVI_binning(bin,j,counter(bin,j)) = dSM_dt_extract(i,j); 

    counter(bin,j) = counter(bin,j) + 1 ;


    if counter(bin,j) == 2000000
    return
    end

        end % i
    
    end  % j 


end

 spatial


end






dSM_dt_NDVI_binning(dSM_dt_NDVI_binning == 0) = NaN ; 
samples_3D_count = sum(~isnan(dSM_dt_NDVI_binning),3,'omitnan') ;  



xmap = mean(dSM_dt_NDVI_binning,3,'omitnan') ; 
 xmap(samples_3D_count < 1000) = NaN ; 
% xmap(samples_3D_count < 100) = NaN ; 

% imagesc(samples_3D_count)


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
caxis([-0.06 0]) %
 ylim([0 1])
 xlim([-0.00 0.5])
colormap(bam_color(1:128,:)) %
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat [m³/m³/day]";
cbr.Label.FontSize = 17 ; 
xticks([0:0.2:0.6])
xlabel('SM [m³/m³]')
ylabel('NDVI [-]')
%title('2000-2001')
pbaspect([0.8 0.8 0.8])
set(gca,'FontSize',17)

print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSM_dt_binned_NDVI_sminterp_01','-dbmp','-r150')
close




%% boxplots of main loss rate



sub5 = subplot(2,3,6) ;
set(sub5, 'DefaultTextFontSize', 16);
Post_cur = boxplot(boxplot_data1NBRSM  ,'Labels',severity_chars_NBR_final, ...
    'Whisker', 1.5, 'Jitter', 0.0001, 'Positions',1:7,'Width',0.65) ;
boxes1 = findobj(gcf,'tag','Box','-and','DisplayName','') ; 
set(Post_cur(7,:),'Visible','off')
set(findobj(gcf,'tag','Box','-and','DisplayName',''), 'Color', [0.25 0.25 0.25],'DisplayName','xx','LineWidth',1.3);
set(findobj(gcf,'tag','Upper Whisker','-and','DisplayName',''), 'Color', [0 0 0],'DisplayName','xx');
set(findobj(gcf,'tag','Lower Whisker','-and','DisplayName',''), 'Color', [0 0 0],'DisplayName','xx');
h_post = findobj(gca,'Tag','Box','-and','DisplayName','');
set(h_post,'DisplayName','MB')
hold on
 for j=1:length(boxes1)
   patches1 =  patch(get(boxes1(j),'XData'),get(boxes1(j),'YData'),col_L,'FaceAlpha',.5);
 end
 hold on

Post_cur = boxplot(boxplot_data1NBRSM  ,'Labels',severity_chars_NBR_final, ...
     'Whisker', 1.5, 'Jitter', 0.0001, 'Positions',1:7,'Width',0.65) ;
set(Post_cur(7,:),'Visible','off')
set(findobj(gcf,'tag','Box','-and','DisplayName',''), 'Color', [0.25 0.25 0.25],'DisplayName','xx','LineWidth',1.3);
set(findobj(gcf,'tag','Upper Whisker','-and','DisplayName',''), 'Color', [0 0 0],'DisplayName','xx');
set(findobj(gcf,'tag','Lower Whisker','-and','DisplayName',''), 'Color', [0 0 0],'DisplayName','xx');
h_post = findobj(gca,'Tag','Box','-and','DisplayName','');
set(h_post,'DisplayName','MB')
xlim([0 8])
 ylim([-0.02 0.02])
%ylim([-50 50])
xtickangle(45)
ylabel('\DeltaSM/\Deltat change [m³/m³/day]','FontSize',16)
%ylabel('\DeltaSM/\Deltat change [%]','FontSize',16)  % % change
hold on
set(gca,'FontSize',16)
% set position
set(sub5,'units','centimeters','position',  [ 2.5+18+11, 15-8-4.25+1, 8,  8])
xlabel('fire severity \DeltaNBR [-]','FontSize',16)








%% use all tiles to get an average SM loss every year. then calculate trends




clear


cd('F:\ESA_CCI\analysis_datasets\spatial_02')
% dd data
dSM_dt_files = string(ls('*dSM_dt_interpsm*')) ;
rowcol_files = string(ls('*rowcol_year_dummy*')) ;
t_start_end_files = string(ls('*t_start*')) ;



sminterp = linspace(0,0.6,60) ; 
f_E_bins = linspace(-0.2,0.2,40) ; 




load('F:\projects\SM_long_term_DDs\data_for_figures\DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
DD_dSM_dt_NDVI_gradient_median_05(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% DD_dSM_dt_NDVI_gradient_median_05(DD_dSM_dt_NDVI_gradient_median_05 < 0) = NaN ; 
load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\DD_dSM_dt_NDVI_gradient_zeppe_ERA')
DD_dSM_dt_NDVI_gradient_zeppe_ERA(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 




years_vec = 1990:2021 ; 
for i = 1:length(years_vec)
    datetime_dummy_years(i) = datetime(strcat('01-Jan-',num2str(years_vec(i)))) ;
end

ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ;

SM_loss_function_trend = NaN (1,40) ;
SM_loss_function_trend_p = NaN (1,40) ;

% spatial = 2050

counter = 1 ; 
SM_loss_function_array = NaN(24,40) ; 
SM_loss_function_array_samples = NaN(24,40) ; 

for years_index = 1:24




    
year_vector_start = datetime_dummy_years(years_index) ; 
year_vector_end = datetime_dummy_years(years_index+1) ; 

year_vector_start_index = find(ESA_CCI_datetime == year_vector_start) ; 
year_vector_end_index = find(ESA_CCI_datetime == year_vector_end) ; 

dSM_dt_dummy_year_array = NaN(700000,40) ; 
counter = 1 ; 

% spatial = 345
for spatial = 1:size(dSM_dt_files,1)

cd('F:\ESA_CCI\analysis_datasets\spatial_02')
  
dummy_dSM_dt = matfile(dSM_dt_files(spatial)) ; 
dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy    ; 
% dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% exclude dainage
dummy_dSM_dt(:,41:end) = [] ; 
dummy_dSM_dt = dummy_dSM_dt .* -100 ; 


% get 2.5 degree rowcol from name
name_files = char(dSM_dt_files(spatial)) ; 
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


   
dummy_t_start_end = matfile(t_start_end_files(spatial)) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    

select_year = dummy_t_start_end(:,1) > year_vector_start_index & dummy_t_start_end(:,2) < year_vector_end_index ; 

dummy_select = dummy_dSM_dt(select_year,:) ; 


dSM_dt_dummy_year_array(counter:counter+size(dummy_select,1)-1,:) = dummy_select ; 
       counter  = counter + length(dummy_select) ; 

       if counter ==  700000
           return
       end

spatial

end


dSM_dt_dummy_year_array_median = median(dSM_dt_dummy_year_array,1,'omitnan') ;
dSM_dt_dummy_year_array_median_samples = sum(~isnan(dSM_dt_dummy_year_array),1) ; 


SM_loss_function_array(years_index,:) = dSM_dt_dummy_year_array_median ;
SM_loss_function_array_samples(years_index,:) = dSM_dt_dummy_year_array_median_samples ;

years_index
end


SM_loss_function_array_ESACCI = SM_loss_function_array ; 
SM_loss_function_array_samples_ESACCI = SM_loss_function_array_samples ; 





save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\SM_loss_function_array_ESACCI','SM_loss_function_array_ESACCI') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\SM_loss_function_array_samples_ESACCI','SM_loss_function_array_samples_ESACCI') ; 





load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\SM_loss_function_array_ESACCI','SM_loss_function_array_ESACCI') ; 
load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\SM_loss_function_array_samples_ESACCI','SM_loss_function_array_samples_ESACCI') ; 
load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\SM_loss_function_array_Zeppe','SM_loss_function_array_Zeppe') ; 



SM_loss_trend = NaN(1,40) ; 
SM_loss_trend_p = NaN(1,40) ; 
cd('E:\Daten Baur\Matlab code')

for i = 1:40

              [m_cur b_cur] = TheilSen([(1:24)' , SM_loss_function_array_ESACCI(:,i)]) ; 
              SM_loss_trend(:,i) = m_cur ;    
              SM_loss_trend_b(:,i) = b_cur ; 
              % Xmap_space_intercept(r,c) = b_cur ;               
              [H,p_value] = Mann_Kendall(SM_loss_function_array_ESACCI(:,i),0.05)  ;   
              SM_loss_trend_p(:,i) = p_value ; 
              


              [m_cur b_cur] = TheilSen([(1:24)' , SM_loss_function_array_Zeppe(:,i)]) ; 
              SM_loss_trend_zeppe(:,i) = m_cur ;    
              % Xmap_space_intercept(r,c) = b_cur ;               
              [H,p_value] = Mann_Kendall(SM_loss_function_array_Zeppe(:,i),0.05)  ;   
              SM_loss_trend_p_zeppe(:,i) = p_value ; 
              SM_loss_trend_b_zeppe(:,i) = b_cur ; 


end

Global_SM_loss_trend_ESACCI    = SM_loss_trend   ; 
Global_SM_loss_trend_ESACCI_b    = SM_loss_trend_b   ; 
Global_SM_loss_trend_ESCACCI_p  = SM_loss_trend_p ;


Global_SM_loss_trend_Zeppe   = SM_loss_trend_zeppe   ; 
Global_SM_loss_trend_Zeppe_b    = SM_loss_trend_b_zeppe   ; 
Global_SM_loss_trend_Zeppe_p  = SM_loss_trend_p_zeppe ;








i = 30

figure
plot(SM_loss_function_array_Zeppe(:,i))
hold on
plot(SM_loss_function_array_ESACCI(:,i))
plot(Global_SM_loss_trend_ESACCI(i) .* (1:24) + Global_SM_loss_trend_ESACCI_b(i))
plot(Global_SM_loss_trend_Zeppe(i) .* (1:24) + Global_SM_loss_trend_Zeppe_b(i))


sminterp(i)


Global_SM_loss_trend_ESACCI(i) 
Global_SM_loss_trend_Zeppe(i)


% check years of ESACCI loss function relative to TWS studies suggesting
% big drying in 2003
years = 1990:2014;

for i = 1:24 

plot(SM_loss_function_array_ESACCI(i,:))
hold on
title(strcat("start year is  ",num2str(years(i))))
pause(2)

end


test = mean(SM_loss_function_array_ESACCI,2,'omitnan') ; 

plot(1990:2013,test) ; 


save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\Global_SM_loss_trend_ESACCI','Global_SM_loss_trend_ESACCI') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\Global_SM_loss_trend_ESCACCI_p','Global_SM_loss_trend_ESCACCI_p') ; 



Fig_Panel = figure('units','centimeters','position',[10 2 30 20])  ;

for i = 1:24
plot(sminterp_40, SM_loss_function_array_Zeppe(i,:))
%  'Color',bam_color(3*i,:)
hold on

end

ylim([0 3])
xlabel('SM [m³/m³]')
ylabel('SM loss [mm/day]')
fontsize(18,'points')

saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Analysis_trends\Zeppe_SM_loss_24year','svg')
close 





Fig_Panel = figure('units','centimeters','position',[10 2 30 20])  ;

for i = 1:24
plot(sminterp_40, SM_loss_function_array_ESACCI(i,:))
%  'Color',bam_color(3*i,:)
hold on

end

ylim([0 3])
xlabel('SM [m³/m³]')
ylabel('SM loss [mm/day]')
fontsize(18,'points')

saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Analysis_trends\ESACCI_SM_loss_24year','svg')
close 






