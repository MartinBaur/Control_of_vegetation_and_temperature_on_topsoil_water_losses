%% MJB do spatial tiles dSM/dt relative to NDVI anomaly. Then do some classification on them


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

% NDVI_binning = linspace(-0.15,0.15,61) ; 
% dSM_dt_binning = linspace(-0.08,0.08,61) ; 




%% get tileplot for each 2.5° bin




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


[Lia Locb] = ismember(NDVI_files_121, NDVI_files_all,'rows') ; 
NDVI_files = NDVI_files_all ; 
NDVI_files(Locb) = [] ; 



ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 

NDVI_binning = linspace(-0.20,0.20,21) ; 
dSM_dt_binning = linspace(-0.08,0.08,41) ; 



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
dSM_dt_NDVI_binning = NaN(20,40,200000) ; 

counter = NaN(40,1) ; 
counter(:,1) = 1 ; 


dummy_dSM_dt = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',dSM_dt_files(spatial))) ; 
NDVI_dummy =  matfile(strcat('F:\ESA_CCI\analysis_datasets\NDVI_sminterp\',NDVI_files_5y(spatial))) ; 



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

 % min dd number
 if size(dummy_dSM_dt,1) < 3000
     continue
 end

% exclude dainage
%dummy_dSM_dt(:,40:end) = NaN ; 


dummy_rowcol = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',rowcol_files(spatial))) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',t_start_end_files(spatial))) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 

% pre filter for time if needed
mask = dummy_t_mean < 4080 | dummy_t_mean > 13210; % 1-Jan-1990
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_NDVI(mask,:) = NaN ; 

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
NDVI_deviations = dummy_NDVI  ; 
% get mean but for each SM bin
dummy_dSM_dt_median = median(dummy_dSM_dt,1,'omitnan') ; 

for bin = 1:20

    cur_bin_min = NDVI_binning(bin) ;
    cur_bin_max = NDVI_binning(bin+1) ;   
    extract_index = NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max ;
    extract_index_row =  any(NDVI_deviations > cur_bin_min & NDVI_deviations < cur_bin_max,2) ; 
    % get dSM/dt anomalies
    extract_index = extract_index(extract_index_row,:) ; 
    dSM_dt_extract = dummy_dSM_dt(extract_index_row,:) - dummy_dSM_dt_median ; 
    dSM_dt_extract(~extract_index) = NaN ; 


    for i =1:size(dSM_dt_extract,1)

    % dSM_dt_extract_40 = interp1(sminterp,dSM_dt_extract(i,:),sminterp) ; 

    dSM_dt_NDVI_binning(bin,:,counter(bin,:)) = dSM_dt_extract(i,:); 

    if all( isnan(dSM_dt_NDVI_binning(bin,:,counter(bin,:))))

    else
    counter(bin,:) = counter(bin,:) + 1 ;
    end

    %counter(bin,:) 

    if counter(bin,:) == 200000       
    return
    end
    
    end



end


dSM_dt_NDVI_binning(dSM_dt_NDVI_binning == 0) = NaN ; 
samples_3D_count = sum(~isnan(dSM_dt_NDVI_binning),3,'omitnan') ;  

xmap = median(dSM_dt_NDVI_binning,3,'omitnan') ; 

%xmap(samples_3D_count < 100) = NaN ; 
save(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_median_spatial_20\dSM_dt_xmap_anomaly5y_spatial_', num2str(row_2_5),'_',num2str(col_2_5))  ,'xmap')
save(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_median_spatial_20\dSM_dt_xmap_anomaly5y_samples_', num2str(row_2_5),'_',num2str(col_2_5))  ,'samples_3D_count')





 spatial

end






load('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_spatial\dSM_dt_xmap_anomaly5y_spatial_52_129.mat')
load('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_spatial\dSM_dt_xmap_anomaly5y_samples_52_129.mat')



dSM_dt_NDVI_binning(dSM_dt_NDVI_binning == 0) = NaN ; 
samples_3D_count = sum(~isnan(dSM_dt_NDVI_binning),3,'omitnan') ;  


xmap = mean(dSM_dt_NDVI_binning,3,'omitnan') ; 
xmap(samples_3D_count < 20) = NaN ; 



Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.01 0.01]) %
 ylim([-0.17 0.17])
% xlim([-0.00 0.6])
colormap(redblue_color) %
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat deviation from mean loss rate [m³/m³/day]";
xticks([0:0.2:0.6])
xlabel('SM [m³/m³]')
ylabel('NDVI deviation from pixel seasonal mean [-]')
%title('2000-2001')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% print('-image','F:\projects\SM_long_term_DDs\figures\02\dSM_dt_dev_from_mean_LAI_dev_season_onepix_20sampl_02','-dbmp','-r150')


close





%% devise some classification scheme .. maybe 4 way color gradient .. or split into pos NDVI and neg NDVI
clear


NDVI_binning = linspace(-0.20,0.20,41) ; 
dSM_dt_binning = linspace(-0.08,0.08,41) ; 

cd('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_spatial')

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


% i = 23 test
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


dummy_xmap = matfile(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_spatial\',fileslist(i))) ; 
dummy_samples = matfile(strcat('F:\ESA_CCI\analysis_datasets\dSm_dt_xmap_spatial\',fileslist_samples(i))) ; 
dummy_xmap = dummy_xmap.xmap ; 
dummy_samples = dummy_samples.samples_3D_count ; 

% kick out less than 20 samples
 dummy_xmap(dummy_samples < 10) = NaN ; 
% restrict to number of valid conditions.. iguess desert has very few but
% spurious high values
  if sum(~isnan(dummy_xmap(:))) < 150
  dummy_xmap(:,:) = NaN ; 
  end

  %%%%%% MJB convert anomalies to mm/day by multiplying by soil layer
  % thickness. Assume 100 mm to start with. 
  dummy_xmap = dummy_xmap .* 100 ; 


 

  % if sum(dummy_samples(:)) < 8000
  % dummy_xmap(:,:) = NaN ; 
  % end


   dummy_xmap_wiener = wiener2(dummy_xmap) ; 

 % dummy_xmap_mean2 = mean(dummy_xmap,2,'omitnan') ; 
 % dummy_xmap_mean2(isnan(dummy_xmap_mean2)) = [] ; 


 % [mdl] = fitlm(dummy_xmap_rows_vec,dummy_xmap(:)) ;
 % linslope_fit =   mdl.Coefficients(2,1) ; 

 % FY gradient is relative to an increase in f_E anomaly
 [FXxmap FYxmap] = gradient(dummy_xmap) ; 


     % figure
     % pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (dummy_xmap)) ; colormap(redblue_color) ; caxis([-0.01 0.01]); shading flat
     % pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, wiener2(dummy_xmap)) ; colormap(redblue_color) ; caxis([-0.01 0.01]); shading flat

% remove smallest NDVI deviation bins
% dummy_xmap(19:22,:) = NaN ; 

% split 
  dummy_xmap_pos = dummy_xmap(21:end,:) ; 
  dummy_xmap_neg = dummy_xmap(1:20,:) ; 

% % pos
% pcolor(sminterp(1:40),NDVI_binning(22:end) - diff(NDVI_binning(1:2))/2, (dummy_xmap_pos)) ; colormap(redblue_color) ; caxis([-0.01 0.01]); shading flat
% %  neg
% pcolor(sminterp(1:40),NDVI_binning(1:20) - diff(NDVI_binning(1:2))/2, (dummy_xmap_neg)) ; colormap(redblue_color) ; caxis([-0.01 0.01]); shading flat

 % if mean(dummy_xmap_pos(:),'omitnan') < 0 



dSM_dt_mean_NDVI_pos(row_2_5,col_2_5) = mean(dummy_xmap_pos(:),'omitnan') ; 
dSM_dt_mean_NDVI_neg(row_2_5,col_2_5) = mean(dummy_xmap_neg(:),'omitnan') ; 
dSM_dt_mean_NDVI_all(row_2_5,col_2_5) = mean(dummy_xmap(:),'omitnan') ; 
% dSM_dt_mean_NDVI_slope(row_2_5,col_2_5) = table2array(linslope_fit) ;
 dSM_dt_mean_NDVI_FYslope(row_2_5,col_2_5) = median(FYxmap(:),'omitnan');


 % end


i

end


% h = pcolor(repmat((1:72)',[1 144]), repmat(1:144,[72 1]), xmap); 
% set(h,'LineStyle','none')
% shading flat


dSM_dt_mean_NDVI_slope(dSM_dt_mean_NDVI_slope < -1*10^-3 | dSM_dt_mean_NDVI_slope > 1*10^-3) = NaN ; 


save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_mean_NDVI_FYslope','dSM_dt_mean_NDVI_FYslope')



xmap = dSM_dt_mean_NDVI_FYslope ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% xmap = DD_dSM_dt_NDVI_gradient_zeppe_ERA ; 

% xmap = dSM_dt_mean_NDVI_slope ; 
% xmap(mask_spatial) = NaN ; 

figure('units','centimeters','position',[10 3 39 20])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
colormap(redblue_color) 
clim([-0.2  0.2])
% caxis([-1*10^-3  1*10^-3])
hcb=colorbar;
ylabel(hcb,'\DeltaSM/\Deltat anomaly [m³/m³/day]','FontSize',18)
title('negative NDVI anomaly')
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
xlim([-1  1])


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





print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSm_dt_anomaly5y_xmap_NDVI_all_03','-dbmp','-r150')
close


% save('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_pos','dSM_dt_mean_NDVI_pos')
% save('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_neg','dSM_dt_mean_NDVI_neg')
% save('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_all','dSM_dt_mean_NDVI_pos')
% 



%% just plot map of median dSM/dt









