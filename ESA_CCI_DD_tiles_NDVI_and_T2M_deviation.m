%% MJB 16.01.2024 new figure Select small bin of SM and then do phase space of T2M and NDVI anomalies


clear



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





%%





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
% dSM_dt_binning = linspace(-0.08,0.08,41) ; 
T2M_binning = linspace(-15,15,41) ; 

% 
% sminterp = linspace(0,1,50) ; 
% f_E_bins = linspace(-0.2,0.2,51) ; 
% T2M_bins = linspace(-15,15,51) ; 




counter = NaN(40,40) ; 
counter(:,:) = 1 ; 


% binning based on NDVI
%dSM_dt_NDVI_binning = NaN(40,40,1000000) ; 
dSM_dt_T2M_NDVI_binning = NaN(40,40,1200000) ; 



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
 dummy_NDVI = NDVI_dummy.NDVI_spatial_sminterp ; 
 dummy_NDVI(:,41:60) = [] ;

 dummy_T2M = dummy_T2M.T2M_spatial_sminterp ; 
 dummy_T2M(:,41:60) = [] ;


dummy_rowcol = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',rowcol_files(spatial))) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',t_start_end_files(spatial))) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 

% pre filter for time if needed
mask = dummy_t_mean < 4080 | dummy_t_start_end(:,1) > 13210; % 1-Jan-1990 until end 2014
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_NDVI(mask,:) = NaN ; 
dummy_T2M(mask,:) = NaN ; 



%% MJB rmeove all data but few bins of SM .. then do phase space of T2M and NDVI anomalies

dummy_dSM_dt = dummy_dSM_dt(:,30:40) ; 
dummy_NDVI = dummy_NDVI(:,30:40) ; 
dummy_T2M = dummy_T2M(:,30:40) ; 


%%



%%%%%%%%%%% convert to mm/day %%%%%%%%%%%
dummy_dSM_dt = dummy_dSM_dt .* 100 ; 

% binning based on NDVI
% dSM_dt_NDVI_binning = NaN(30,60) ; 
%NDVI_mean = mean(dummy_NDVI,'omitnan') ; 
NDVI_deviations = dummy_NDVI  ; 
T2M_deviations = dummy_T2M ; 
% get mean but for each SM bin
dummy_dSM_dt_mean = median(dummy_dSM_dt(:),1,'omitnan') ; 


for bin_T2M = 1:40

    cur_bin_min_T2M = T2M_binning(bin_T2M) ;
    cur_bin_max_T2M = T2M_binning(bin_T2M+1) ; 

for bin_NDVI = 1:40

    cur_bin_min_NDVI = NDVI_binning(bin_NDVI) ;
    cur_bin_max_NDVI = NDVI_binning(bin_NDVI+1) ;   

    % standard bining 
    extract_index = NDVI_deviations > cur_bin_min_NDVI & NDVI_deviations < cur_bin_max_NDVI & ...
                    T2M_deviations > cur_bin_min_T2M & T2M_deviations < cur_bin_max_T2M;
   
    dSM_dt_extract = dummy_dSM_dt(extract_index) ; 

    % get dSM/dt anomalies
    dSM_dt_extract = dSM_dt_extract - dummy_dSM_dt_mean ; 
    dSM_dt_T2M_NDVI_binning(bin_NDVI,bin_T2M,counter(bin_NDVI,bin_T2M):counter(bin_NDVI,bin_T2M)+length(dSM_dt_extract)-1) = dSM_dt_extract; 

    counter(bin_NDVI,bin_T2M) = counter(bin_NDVI,bin_T2M) + length(dSM_dt_extract) ; 

    if counter(bin_NDVI,bin_T2M)  >= 1500000
        return
    end


end

end

 spatial

end





dSM_dt_T2M_NDVI_binning(dSM_dt_T2M_NDVI_binning == 0) = NaN ; 
samples_3D_count = sum(~isnan(dSM_dt_T2M_NDVI_binning),3,'omitnan') ;  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dSM_dt_NDVI_T2M_anomaly_2D_sampling = sum(~isnan(dSM_dt_T2M_NDVI_binning),3) ; 
dSM_dt_NDVI_T2M_anomaly_2D = median(dSM_dt_T2M_NDVI_binning,3,'omitnan');
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_NDVI_T2M_anomaly_2D_sampling','dSM_dt_NDVI_T2M_anomaly_2D_sampling') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_NDVI_T2M_anomaly_2D','dSM_dt_NDVI_T2M_anomaly_2D') ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xmap = dSM_dt_NDVI_T2M_anomaly_2D ; 


xmap = mean(dSM_dt_T2M_NDVI_binning,3,'omitnan') ; 
 % xmap(samples_3D_count < 1000) = NaN ; 
 % xmap(samples_3D_count < 100) = NaN ; 

% imagesc(samples_3D_count)


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(T2M_binning(2:end) - diff(T2M_binning(1:2))/2  , NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.5 0.5]) %
  ylim([-0.2 0.2])
  xlim([-10 10])
colormap(redblue_color) %
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat anomaly [m³/m³/day]";
cbr.Label.FontSize = 17 ; 
xlabel('T2M anomaly [°K]')
ylabel('NDVI anomaly [-]')
 title('SM = 0.3-0.4 [m³/m³]')
 yline(0)
 xline(0)

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



print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSM_dt_anomaly_NDVI_T2M_combined_binning_sm03_04','-dbmp','-r150')
close

































