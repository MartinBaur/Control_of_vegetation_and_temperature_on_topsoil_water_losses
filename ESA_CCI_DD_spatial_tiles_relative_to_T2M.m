%% MJB 26.07.2023 ESA CCI dSM/dt relative to T2M anomalies if effect try to remove from NDVI effect

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
sminterp40 = (0.01:0.01:0.4)' ; 


T2M_binning = linspace(-8,8,61) ; 
dSM_dt_binning = linspace(-0.08,0.08,61) ; 



%%
clear

cd('F:\ESA_CCI\analysis_datasets\spatial_02')
% dd data


dSM_dt_files = string(ls('*dSM_dt_interpsm*')) ;
rowcol_files = string(ls('*rowcol_year_dummy*')) ;
t_start_end_files = string(ls('*t_start*')) ;

% monthly LAI will be only until 2014
cd('F:\ESA_CCI\analysis_datasets\T2M_sminterp')

T2M_files_all = string(ls('*T2M*')) ;
T2M_files_121 = string(ls('*T2M*sminterp_121_*')) ;
[Lia Locb] = ismember(T2M_files_121,T2M_files_all);

T2M_files = T2M_files_all ; 
T2M_files(Locb) = [] ; 


ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 


counter = NaN(60,1) ; 
counter(:,1) = 1 ; 



ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 


T2M_binning = linspace(-15,15,41) ; 
dSM_dt_binning = linspace(-0.08,0.08,41) ; 


% ESA_CCI_datetime(4080)
% ESA_CCI_datetime(end)
years_vec = 1990:2021 ; 
for i = 1:length(years_vec)
    datetime_dummy_years(i) = datetime(strcat('01-Jan-',num2str(years_vec(i)))) ;
end
datetime_dummy_years(end) = ESA_CCI_datetime(end) ; 


counter = NaN(40,40) ; 
counter(:,:) = 1 ; 
% binning based on NDVI
dSM_dt_T2M_binning = NaN(40,40,600000) ; 



% spatial = 2450
for spatial = 1:length(dSM_dt_files)




% get drydown data
cd('F:\ESA_CCI\analysis_datasets\spatial_02')   
dummy_dSM_dt = matfile(dSM_dt_files(spatial)) ; 
cd('F:\ESA_CCI\analysis_datasets\T2M_sminterp')
T2M_dummy =  matfile(T2M_files(spatial)) ; 

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


dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy   ; 
% dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% new condition removing very short DDs?
% dummy_dSM_dt(sum(~isnan(dummy_dSM_dt),2) < 3 , :) = NaN ; 
% kick out all drainage to save memory. reduced to 40 length from 60
dummy_dSM_dt(:,41:60) = [] ;

 % dummy_NDVI = NDVI_dummy.dummy_NDVI_monthly ; 
 dummy_T2M = T2M_dummy.T2M_spatial_sminterp ; 
 dummy_T2M(:,41:60) = [] ;

% exclude dainage
%dummy_dSM_dt(:,40:end) = NaN ; 

cd('F:\ESA_CCI\analysis_datasets\spatial_02')   
dummy_rowcol = matfile(rowcol_files(spatial)) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(t_start_end_files(spatial)) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 

% pre filter for time if needed double check dates again .. 
mask = dummy_t_start_end(:,2) < 4080 | dummy_t_start_end(:,1) > 13210;
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_T2M(mask,:) = NaN ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert to mm/day by multiplying with 100 mmm soil depth
dummy_dSM_dt = dummy_dSM_dt .* 100 ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % temporal condition to 1 year
% time_start = find(ESA_CCI_datetime == datetime('01-Jan-2013')) ; 
% time_end = find(ESA_CCI_datetime == datetime('01-Jan-2014')) ; 
% % time_start = find(ESA_CCI_datetime == year_vector_start) ; 
% % time_end = find(ESA_CCI_datetime == year_vector_end) ; 
% 
% 
% 
% mask = dummy_t_mean < time_start |  dummy_t_mean > time_end ; % 1-Jan-1990
% dummy_t_mean(mask) = NaN ; 
% dummy_dSM_dt(mask,:) = NaN ; 
% dummy_rowcol(mask,:) = NaN ; 
% dummy_T2M(mask,:) = NaN ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(any(~isnan(dummy_dSM_dt)))
% return
% end



% binning based on NDVI
% binning based on NDVI
% dSM_dt_NDVI_binning = NaN(30,60) ; 
%NDVI_mean = mean(dummy_NDVI,'omitnan') ; 
NDVI_deviations = dummy_T2M  ; 
% get mean but for each SM bin
dummy_dSM_dt_median = median(dummy_dSM_dt,1,'omitnan') ; 


for bin = 1:40

    cur_bin_min = T2M_binning(bin) ;
    cur_bin_max = T2M_binning(bin+1) ;   
    extract_index = NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max ;
    extract_index_row =  any(NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max,2) ; 
    % get dSM/dt anomalies
    extract_index = extract_index(extract_index_row,:) ; 
    dSM_dt_extract = dummy_dSM_dt(extract_index_row,:) - dummy_dSM_dt_median ; 
    dSM_dt_extract(~extract_index) = NaN ; 



    for j = 1:size(dSM_dt_extract,2)

        if all(isnan( dSM_dt_extract(:,j)))
            continue
        end

        for i =1:size(dSM_dt_extract,1)

        if (isnan( dSM_dt_extract(i,j)))
            continue
        end
    
    dSM_dt_T2M_binning(bin,j,counter(bin,j)) = dSM_dt_extract(i,j); 

    counter(bin,j) = counter(bin,j) + 1 ;


    if counter(bin,j) == 600000
    return
    end

        end % i
    
    end  % j 



end

% dSM_dt_T2M_binning(dSM_dt_T2M_binning == 0) = NaN ; 
% samples_3D_count = sum(~isnan(dSM_dt_T2M_binning),3,'omitnan') ;  
% 
% xmap = mean(dSM_dt_T2M_binning,3,'omitnan') ; 
% %xmap(samples_3D_count < 100) = NaN ; 
% save(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_spatial\dSM_dt_xmap_spatial_', num2str(row_2_5),'_',num2str(col_2_5))  ,'xmap')
% save(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_spatial\dSM_dt_xmap_samples_', num2str(row_2_5),'_',num2str(col_2_5))  ,'samples_3D_count')



 spatial

end





dSM_dt_T2M_binning(dSM_dt_T2M_binning == 0) = NaN ; 
samples_3D_count = sum(~isnan(dSM_dt_T2M_binning),3,'omitnan') ;  
% imagesc(samples_3D_count)


dSM_dt_T2M_binning_median = median(dSM_dt_T2M_binning,3,'omitnan') ; 
dSM_dt_T2M_binning_median_samples = samples_3D_count ; 


save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_T2M_binning_median','dSM_dt_T2M_binning_median')
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_T2M_binning_median_samples','dSM_dt_T2M_binning_median_samples')



xmap = median(dSM_dt_T2M_binning,3,'omitnan') ; 
xmap(samples_3D_count < 1000) = NaN ; 



Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp40,T2M_binning(2:end) - diff(T2M_binning(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.2 0.2]) %
 ylim([-12 12])
 xlim([-0.00 0.5])
colormap(redblue_color) %
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat deviation from mean loss rate [m³/m³/day]";
xticks([0:0.2:0.6])
xlabel('SM [m³/m³]')
ylabel('T2M deviation from pixel seasonal mean [K]')
% title('tiles 500-1500')
set(gca,'FontSize',17)
pbaspect([1 1 1])


print('-image','F:\projects\SM_long_term_DDs\figures\02\dSM_dt_dev_from_mean_T2M_dev_2013_2014','-dbmp','-r150')


close








%% xmap for T2M anomalies






clear

cd('F:\ESA_CCI\analysis_datasets\spatial_02')
% dd data
dSM_dt_files = string(ls('*dSM_dt_interpsm*')) ;
rowcol_files = string(ls('*rowcol_year_dummy*')) ;
t_start_end_files = string(ls('*t_start*')) ;


% monthly LAI will be only until 2014
cd('F:\ESA_CCI\analysis_datasets\T2M_sminterp')
T2M_files_all = string(ls('*T2M*')) ;
T2M_files_121 = string(ls('*T2M*sminterp_121_*')) ;
[Lia Locb] = ismember(T2M_files_121,T2M_files_all);

T2M_files = T2M_files_all ; 
T2M_files(Locb) = [] ; 



ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 

% NDVI_binning = linspace(-0.20,0.20,41) ; 
dSM_dt_binning = linspace(-0.08,0.08,41) ; 
T2M_binning = linspace(-15,15,21) ; 


sxity_to_thirty = 0:2:60 ; sxity_to_thirty(1) = 1 ; 



% ESA_CCI_datetime(4080)
% ESA_CCI_datetime(end)
years_vec = 1990:2021 ; 
for i = 1:length(years_vec)
    datetime_dummy_years(i) = datetime(strcat('01-Jan-',num2str(years_vec(i)))) ;
end
datetime_dummy_years(end) = ESA_CCI_datetime(end) ; 





% spatial = 2450
for spatial = 1:length(dSM_dt_files)


% binning based on NDVI
dSM_dt_T2M_binning = NaN(20,40,200000) ; 

counter = NaN(40,1) ; 
counter(:,1) = 1 ; 


dummy_dSM_dt = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',dSM_dt_files(spatial))) ; 
T2M_dummy =  matfile(strcat('F:\ESA_CCI\analysis_datasets\T2M_sminterp\',T2M_files(spatial))) ; 



% get 2.5 degree rowcol from name
name_files = char(dSM_dt_files(spatial)) ; 
% grep row and col
pat = digitsPattern;
row_col_name = extract(name_files,pat) ; 
row_2_5 = str2double(row_col_name{1}) ; 
col_2_5 = str2double(row_col_name{2}) ;     


dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy   ; 
% dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% new condition removing very short DDs?
% dummy_dSM_dt(sum(~isnan(dummy_dSM_dt),2) < 3 , :) = NaN ; 
% kick out all drainage to save memory. reduced to 40 length from 60
dummy_dSM_dt(:,41:60) = [] ;


 % dummy_NDVI = NDVI_dummy.dummy_NDVI_monthly ; 
 dummy_T2M = T2M_dummy.T2M_spatial_sminterp ; 
 dummy_T2M(:,41:60) = [] ;

% exclude dainage
%dummy_dSM_dt(:,40:end) = NaN ; 


dummy_rowcol = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',rowcol_files(spatial))) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',t_start_end_files(spatial))) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 

% pre filter for time if needed
mask = dummy_t_start_end(:,2) < 4080 | dummy_t_start_end(:,2) > 13210 ;
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_T2M(mask,:) = NaN ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporal condition to 1 year
% time_start = find(ESA_CCI_datetime == datetime('01-Jan-2013')) ; 
% time_end = find(ESA_CCI_datetime == datetime('01-Jan-2014')) ; 
% time_start = find(ESA_CCI_datetime == year_vector_start) ; 
% time_end = find(ESA_CCI_datetime == year_vector_end) ; 



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
NDVI_deviations = dummy_T2M  ; 
% get mean but for each SM bin
dummy_dSM_dt_median = median(dummy_dSM_dt,1,'omitnan') ; 

for bin = 1:20

    cur_bin_min = T2M_binning(bin) ;
    cur_bin_max = T2M_binning(bin+1) ;   
    extract_index = NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max ;
    extract_index_row =  any(NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max,2) ; 
    % get dSM/dt anomalies
    extract_index = extract_index(extract_index_row,:) ; 
    dSM_dt_extract = dummy_dSM_dt(extract_index_row,:) - dummy_dSM_dt_median ; 
    dSM_dt_extract(~extract_index) = NaN ; 


    for i =1:size(dSM_dt_extract,1)

    % dSM_dt_extract_40 = interp1(sminterp,dSM_dt_extract(i,:),sminterp) ; 

    dSM_dt_T2M_binning(bin,:,counter(bin,:)) = dSM_dt_extract(i,:); 

    if all( isnan(dSM_dt_T2M_binning(bin,:,counter(bin,:))))

    else
    counter(bin,:) = counter(bin,:) + 1 ;
    end

    %counter(bin,:) 

    if counter(bin,:) == 200000       
    return
    end
    
    end



end


dSM_dt_T2M_binning(dSM_dt_T2M_binning == 0) = NaN ; 
samples_3D_count = sum(~isnan(dSM_dt_T2M_binning),3,'omitnan') ;  

xmap = mean(dSM_dt_T2M_binning,3,'omitnan') ; 
%xmap(samples_3D_count < 100) = NaN ; 
save(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_T2M_spatial_20\dSM_dt_xmap_anomaly5y_spatial_', num2str(row_2_5),'_',num2str(col_2_5))  ,'xmap')
save(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_T2M_spatial_20\dSM_dt_xmap_anomaly5y_samples_', num2str(row_2_5),'_',num2str(col_2_5))  ,'samples_3D_count')





 spatial

end














%% devise some classification scheme .. maybe 4 way color gradient .. or split into pos NDVI and neg NDVI
clear


T2M_binning = linspace(-10,10,41) ; 
dSM_dt_binning = linspace(-0.08,0.08,41) ; 

% cd('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_T2M_spatial')
cd('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_T2M_spatial_20')



fileslist = string(ls('*dSM_dt_xmap_anomaly5y_spatial*')) ; 
fileslist_samples = string(ls('*dSM_dt_xmap_anomaly5y_samples*')) ; 

lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 

lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 

dSM_dt_mean_NDVI_pos = NaN(72,144) ; 
dSM_dt_mean_NDVI_neg = NaN(72,144) ; 
dSM_dt_mean_NDVI_all = NaN(72,144) ; 
dSM_dt_mean_NDVI_slope = NaN(72,144) ; 
dSM_dt_mean_NDVI_FYslope = NaN(72,144) ; 

cd('F:\ESA_CCI\analysis_datasets\analysis_02')

load('DD_SM_sample_count_1990_now_2_5.mat')
mask_spatial = DD_SM_sample_count_1990_now_2_5 < 10000 ; 

cd('F:\ESA_CCI\analysis_datasets')

%   i = 750
dummy_xmap_rows = repmat((1:40)',[1 40]) ; 
dummy_xmap_rows_vec = dummy_xmap_rows(:) ; 

for i = 1:length(fileslist)


% get 2.5 degree rowcol from name
name_files = char(fileslist(i)) ; 
% grep row and col
pat = digitsPattern;
row_col_name = extract(name_files,pat) ; 
row_2_5 = str2double(row_col_name{2}) ; 
col_2_5 = str2double(row_col_name{3}) ;   
% lats_2_5(row_2_5,1) + 1.25
% lons_2_5(1,col_2_5) - 1.25

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new condition to check if too northern
lat_2_5_dummy = lats_2_5(row_2_5,col_2_5) ; 
lon_2_5_dummy = lons_2_5(row_2_5,col_2_5) ; 

if lat_2_5_dummy > 70
    continue
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



dummy_xmap = matfile(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_T2M_spatial\',fileslist(i))) ; 
dummy_samples = matfile(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_T2M_spatial\',fileslist_samples(i))) ; 
dummy_xmap = dummy_xmap.xmap ; 
dummy_samples = dummy_samples.samples_3D_count ; 

% kick out less than 20 samples
 dummy_xmap(dummy_samples < 10) = NaN ; 
% restrict to number of valid conditions.. iguess desert has very few but
% spurious high values
if sum(~isnan(dummy_xmap(:))) < 150
dummy_xmap(:,:) = NaN ; 
end

 % if sum(dummy_samples(:)) < 8000
 % dummy_xmap(:,:) = NaN ; 
 % end

%%%%%%%%%%%% change loss rate to mm/d %%%%%%%%%%%%%
dummy_xmap = dummy_xmap .*100 ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   dummy_xmap = wiener2(dummy_xmap) ; 

 dummy_xmap_mean2 = mean(dummy_xmap,2,'omitnan') ; 
 dummy_xmap_mean2(isnan(dummy_xmap_mean2)) = [] ; 



 [mdl] = fitlm(dummy_xmap_rows_vec,dummy_xmap(:)) ;
 linslope_fit =   mdl.Coefficients(2,1) ; 


%     figure
%     pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (dummy_xmap)) ; colormap(redblue_color) ; caxis([-0.01 0.01]); shading flat

% remove smallest NDVI deviation bins
% dummy_xmap(19:22,:) = NaN ; 

% split 
 dummy_xmap_pos = dummy_xmap(21:end,:) ; 
 dummy_xmap_neg = dummy_xmap(1:20,:) ; 
 [FXxmap, FYxmap] = gradient(dummy_xmap) ; 


% % pos
% pcolor(sminterp(1:40),NDVI_binning(22:end) - diff(NDVI_binning(1:2))/2, (dummy_xmap_pos)) ; colormap(redblue_color) ; caxis([-0.01 0.01]); shading flat
% %  neg
% pcolor(sminterp(1:40),NDVI_binning(1:20) - diff(NDVI_binning(1:2))/2, (dummy_xmap_neg)) ; colormap(redblue_color) ; caxis([-0.01 0.01]); shading flat

 % if mean(dummy_xmap_pos(:),'omitnan') < 0 



dSM_dt_mean_NDVI_pos(row_2_5,col_2_5) = mean(dummy_xmap_pos(:),'omitnan') ; 
dSM_dt_mean_NDVI_neg(row_2_5,col_2_5) = mean(dummy_xmap_neg(:),'omitnan') ; 
dSM_dt_mean_NDVI_all(row_2_5,col_2_5) = mean(dummy_xmap(:),'omitnan') ; 
dSM_dt_mean_NDVI_slope(row_2_5,col_2_5) = table2array(linslope_fit) ;
dSM_dt_mean_NDVI_FYslope(row_2_5,col_2_5) = median(FYxmap(:),'omitnan');


 % end


i

end


% h = pcolor(repmat((1:72)',[1 144]), repmat(1:144,[72 1]), xmap); 
% set(h,'LineStyle','none')
% shading flat
dSM_dt_mean_T2M_FYslope = dSM_dt_mean_NDVI_FYslope ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_mean_T2M_FYslope','dSM_dt_mean_T2M_FYslope')



dSM_dt_mean_NDVI_slope(dSM_dt_mean_NDVI_slope < -1*10^-3 | dSM_dt_mean_NDVI_slope > 1*10^-3) = NaN ; 

 xmap = dSM_dt_mean_NDVI_FYslope ;
% xmap = dSM_dt_mean_NDVI_slope ; 
% xmap(mask_spatial) = NaN ; 

figure('units','centimeters','position',[10 3 39 20])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
colormap(redblue_color) 
 clim([-0.1  0.1])
% caxis([-1*10^-3  1*10^-3])
hcb=colorbar;
ylabel(hcb,'\DeltaSM/\Deltat anomaly [m³/m³/day]','FontSize',18)
title('T2M anomaly negative')
xlabel('longitude')
ylabel('latitude')
plot(CoastlineLon, CoastlineLat,'Color','k');
hold on
% [rowfind colfind] = find(VIP_NDVI_corr_all_p < 0.05) ; 
% crosses1 = plot(lons_2_5(1,colfind),lats_2_5(rowfind,1),'.k','MarkerSize',7) ; 
set(gca,'FontSize',18)
axes('units','centimeters','Position',[ 7, 5, 4,  4])
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.6  0.6])


% add arrows for description
 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[37,  9.5,   0,  9]) ; 
 set(arrow2,'Position',[37 , 10,     0,  -7.85]) ; 
 
% add textbox
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss                               slower SM loss', ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',18,'Units','centimeters');
set(textbox1,'Position',[37.5,  10.25,   0,  9]) ; 
set(gca,'Box','on')





print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSm_dt_anomaly5y_xmap_T2M_neg','-dbmp','-r150')
close













