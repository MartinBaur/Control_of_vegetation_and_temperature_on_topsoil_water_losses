%% MJB 02-05-2023 ..Do plots for soil moisture loss function project

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




%% Do global all data loss function plot both with ESA CCI and ERA cut all analysis to 1990 until now
% think about spatial limits as well
clear
cd('F:\ESA_CCI\2D')
ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 

load('dSM_dt_interpsm_2D_array.mat')
load('t_start_end_2D_array.mat')
load('rowcol_2D_array.mat')


% exclude all data before 1990
dSM_dt_interpsm_2D_array(t_start_end_2D_array(:,1) < 4080,:) = NaN ; 


dSM_dt_median_1990_now =   median(dSM_dt_interpsm_2D_array,1,'omitnan') ;
dSM_dt_mean_1990_now =   mean(dSM_dt_interpsm_2D_array,1,'omitnan') ;

dSM_dt_25prct_1990_now =   prctile(dSM_dt_interpsm_2D_array,25,1) ; 
dSM_dt_75prct_1990_now =   prctile(dSM_dt_interpsm_2D_array,75,1) ; 

samples_across_SM = sum(~isnan(dSM_dt_interpsm_2D_array),1,'omitnan') ; 


test_boxes = vertcat(dSM_dt_25prct_1990_now,dSM_dt_median_1990_now,dSM_dt_75prct_1990_now) ; 
test_boxes(:,samples_across_SM < 1000) = NaN ; 
boxplot(test_boxes,'Whisker',Inf)


plot(dSM_dt_median_1990_now)
hold on
plot(dSM_dt_25prct_1990_now)
plot(dSM_dt_75prct_1990_now)






% median and percentiles of loss function

save('F:\ESA_CCI\analysis_datasets\dSM_dt_median_1990_now','dSM_dt_median_1990_now')
save('F:\ESA_CCI\analysis_datasets\dSM_dt_mean_1990_now','dSM_dt_mean_1990_now')
save('F:\ESA_CCI\analysis_datasets\samples_across_SM','samples_across_SM')



%% plot one global 1990-2020 SM loss function
pre_index = 2:60 ; 


% crop season index?

prct_25 = dSM_dt_25prct_1990_now(pre_index) ;
prct_50 = dSM_dt_median_1990_now(pre_index) ;
prct_75 = dSM_dt_75prct_1990_now(pre_index) ;


prctlies_box = vertcat(prct_25,prct_50,prct_75) ; 

boxplot(prctlies_box,'Whisker',0) ; 



bam_mean = mean(colormap(bam_color(1:128,:)) ,1,'omitnan')

figure('units','centimeters','position',[10 3 39 20])  ;
linepre = plot(sminterp(pre_index),prct_50,'o-','LineWidth',1.5,'Color',col_C) ;
hold on
%
x2 = [sminterp(pre_index)', fliplr(sminterp(pre_index)')];
inBetween = [prct_50, fliplr(prct_25 )];
fillpre = fill(x2, inBetween, col_C,'FaceAlpha',0.25);
x2 = [sminterp(pre_index)', fliplr(sminterp(pre_index)')];
inBetween = [prct_50, fliplr(prct_75)];
fill(x2, inBetween, col_C,'FaceAlpha',0.25);
plot(sminterp(pre_index),prct_50,'-','LineWidth',1.5,'Color',col_C)
% ylim([-0.4 0])
% xlim([0 0.6])
line25 = plot(sminterp(pre_index),prct_25,'-','LineWidth',1.5,'Color',col_C) ;
line75 = plot(sminterp(pre_index),prct_75,'-','LineWidth',1.5,'Color',col_C) ;
%kolsmi_post = plot(sminterp(pre_index),prct_50,'','MarkerSize',15) ;
xlabel('soil moisture [m³/m³]')
ylabel('\DeltaSM/\Deltat [m³/m³/day]')
%title('Soil moisture loss function 1990-2020')
% add indicators
drainage_thresh = xline(0.44,'--','LineWidth',1.7) ; 
annotation('textbox',[.29 .17 .27 .06],'String','Stage II and stage II SM loss','FontSize',22)
annotation('textbox',[.76 .17 .1 .06],'String','Drainage','FontSize',22)
set(gca,'FontSize',18)

print('-image','F:\projects\SM_long_term_DDs\figures\SM_loss_function_1990_2020','-dbmp','-r150')
close




%% use 2.5 degree spatial dataset to find a rule for pixels we don't want to count
clear

cd('F:\ESA_CCI\analysis_datasets\spatial')
dSM_dt_files = string(ls('*dSM*')) ;
t_start_end_files = string(ls('*t_start*')) ;
rowcol_files = string(ls('*rowcol*')) ;

DD_SM_sample_count_1990_now_2_5 = NaN(72,144,1) ; 
DD_dSM_dt_mean_1990_now_2_5= NaN(72,144,1) ; 

for spatial = 1:size(dSM_dt_files,1)

cd('F:\ESA_CCI\analysis_datasets\spatial')   
dummy_dSM_dt = matfile(dSM_dt_files(spatial)) ; 
dummy_dSM_dt = dummy_dSM_dt.dSM_dt_dummy ; 
% is this needed? 
% dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   

% remove drainage
dummy_dSM_dt(:,40:end) = NaN ;

dummy_rowcol = matfile(rowcol_files(spatial)) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(t_start_end_files(spatial)) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 


% pre filter for time if needed
mask = dummy_t_mean < 4080 ; 
mask = dummy_t_mean > 13210 ; 
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 

% convert to positive mm/day loss
dummy_dSM_dt = dummy_dSM_dt .* -100 ;


% get 2.5 degree rowcol from name
name = char(dSM_dt_files(spatial)) ; 
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



DD_sample_count = sum( ~all(isnan(dummy_dSM_dt),2) )  ; 
DD_dSM_dt_mean = median(dummy_dSM_dt(:),'omitnan') ; 

DD_SM_sample_count_1990_now_2_5(row_2_5,col_2_5,:) = DD_sample_count ;
DD_dSM_dt_mean_1990_now_2_5(row_2_5,col_2_5,:) = DD_dSM_dt_mean ;


spatial

end


DD_SM_sample_count_1990_2014_2_5_nodrain = DD_SM_sample_count_1990_now_2_5 ;
DD_dSM_dt_mean_1990_now_2014_nodrain = DD_dSM_dt_mean_1990_now_2_5 ; 


save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\DD_SM_sample_count_1990_2014_2_5_nodrain','DD_SM_sample_count_1990_2014_2_5_nodrain')
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\DD_dSM_dt_mean_1990_now_2014_nodrain','DD_dSM_dt_mean_1990_now_2014_nodrain')





DD_SM_sample_count_1990_now_2_5(DD_SM_sample_count_1990_now_2_5 < 1000) = NaN ; 

F:\projects\SM_long_term_DDs\data_for_figures

%% plot
cd('F:\ESA_CCI\analysis_datasets\analysis_02')
load('dSM_dt_2_5_grow_season_no_mask.mat')
load('dSM_dt_2_5_non_grow_season_no_mask.mat')
load('DD_SM_sample_count_1990_now_2_5.mat')


mask = DD_SM_sample_count_1990_now_2_5 < 10000 ; 
mask = repmat(mask,[1 1 60]) ; 

dSM_dt_2_5_grow_season(mask) = NaN ; 
dSM_dt_2_5_non_grow_season(mask) = NaN ; 


xmap =  dSM_dt_2_5_grow_season - dSM_dt_2_5_non_grow_season ;
xmap = mean(xmap,3,'omitnan') ; 

figure('units','centimeters','position',[10 3 39 20])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
%h = pcolor(lons_5 - 5, lats_5 + 5, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
%colormap((tokyo_color)) 
%colormap(bam_color(1:128,:)) 
%colormap(bam_color) 
colormap(redblue_color)
caxis([-0.03 0.03])  
%caxis([-0.005 0.005])  
hcb=colorbar;
ylabel(hcb,'\DeltaSM/\Deltat [m³/m³/day]','FontSize',17)
xlabel('longitude')
ylabel('latitude')
%title('growing season - non growing season')
plot(CoastlineLon, CoastlineLat,'Color','k');
set(gca,'FontSize',17)
print('-image','F:\projects\SM_long_term_DDs\figures\02\dSM_dt_maps_grow_minus_non_grow','-dbmp','-r150')
close




%% just average loss here form spatial datasets


clear 
ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 
% 4080 is 1st Jan 1990


cd('F:\ESA_CCI\analysis_datasets\spatial')
dSM_dt_files = string(ls('*dSM*')) ;
t_start_end_files = string(ls('*t_start*')) ;
rowcol_files = string(ls('*rowcol*')) ;


dSM_dt_2_5 = NaN(72,144,60) ; 
sample_count_array = NaN(72,144,60) ; 

for spatial = 1:size(dSM_dt_files,1)

cd('F:\ESA_CCI\analysis_datasets\spatial')   
dummy_dSM_dt = matfile(dSM_dt_files(spatial)) ; 
dummy_dSM_dt = dummy_dSM_dt.dSM_dt_dummy ; 
dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   


dummy_rowcol = matfile(rowcol_files(spatial)) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(t_start_end_files(spatial)) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 


% pre filter for time if needed later than 1990
mask = dummy_t_mean < 4080 ; 
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 


% get 2.5 degree rowcol from name
name = char(dSM_dt_files(spatial)) ; 
% grep row and col
pat = digitsPattern;
row_col_name = extract(name,pat) ; 
row_2_5 = str2double(row_col_name{1}) ; 
col_2_5 = str2double(row_col_name{2}) ;     
   

dummy_dSM_dt(:,sum(~isnan(dummy_dSM_dt),1) < 1) = NaN ; 
sample_count_array(row_2_5,col_2_5,:) = sum(~isnan(dummy_dSM_dt),1) ; 
dSM_dt_2_5(row_2_5,col_2_5,:) = mean(dummy_dSM_dt,1,'omitnan') ;


   spatial
   
end



% dSM_dt_mean_1990_now_02 = NaN(1,1,60) ;
% dSM_dt_mean_1990_now_02(1,1,:) = dSM_dt_mean_1990_now ;
% Mean_full = repmat(dSM_dt_mean_1990_now_02(15:25),[72 144 1]) ; 




xmap =  dSM_dt_2_5 ; 
xmap(sample_count_array < 1000) = NaN ; 
xmap = mean(xmap,3,'omitnan') ; 
xmap(DD_SM_sample_count_1990_now_2_5 < 10000) = NaN ; 


xmap = DD_SM_sample_count_1990_now_2_5 ; 
xmap(DD_SM_sample_count_1990_now_2_5 < 10000) = NaN ; 


figure('units','centimeters','position',[10 3 39 20])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
%h = pcolor(lons_5 - 5, lats_5 + 5, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
%colormap((tokyo_color)) 
colormap(bam_color(129:end,:)) 
%colormap(bam_color) 
%colormap(redblue_color)
% caxis([-0.05 0.00])  
 caxis([10000 40000])  
% caxis([-0.07 0.000])  
hcb=colorbar;
% ylabel(hcb,'\DeltaSM/\Deltat [m³/m³/day]','FontSize',17)
ylabel(hcb,'Number of observed drydowns [-]','FontSize',17)
xlabel('longitude')
ylabel('latitude')
% title('growing season - non growing season')
plot(CoastlineLon, CoastlineLat,'Color','k');
set(gca,'FontSize',17)
print('-image','F:\projects\SM_long_term_DDs\figures\02\sampling_number','-dbmp','-r150')
close





%% visualize individual dd functions
xmap =  dSM_dt_2_5 ; 
xmap(sample_count_array < 1000) = NaN ; 
xmap_mean = mean(xmap,3,'omitnan') ; 
xmap_mean(DD_SM_sample_count_1990_now_2_5 < 10000) = NaN ; 
imagesc(xmap_mean)

[xs ys] = getpts() ; xs = round(xs) ; ys = round(ys) ; 
close

for i = 1:length(xs)
 subplot(1,2,1)
 plot(sminterp,squeeze(xmap(ys(i),xs(i),:)))
 hold on
 ylim([-0.3 0])
 xlim([0 0.6])
 subplot(1,2,2)
 plot(sminterp,squeeze(sample_count_array(ys(i),xs(i),:)))
 hold on

end


%% visualize individual pixel nongrow vs grow


xmap =  dSM_dt_2_5_grow_season - dSM_dt_2_5_non_grow_season ;
xmap = mean(xmap,3,'omitnan') ; 
imagesc(xmap) ;
colormap(redblue_color) ;
clim([-0.03 0.03]) ; 

[xs ys] = getpts() ; xs = round(xs) ; ys = round(ys) ; 
close

for i = 1:length(xs)

 plot(sminterp,squeeze(dSM_dt_2_5_grow_season(ys(i),xs(i),:)))
 hold on
 ylim([-0.3 0])
 xlim([0 0.6])
 plot(sminterp,squeeze(dSM_dt_2_5_non_grow_season(ys(i),xs(i),:)))
 legend('grow','nongrow')


end


%% broken stick results 
clear

cd('F:\ESA_CCI\analysis_datasets\analysis_02')
load('dSM_dt_2_5_grow_season_no_mask.mat')
load('dSM_dt_2_5_non_grow_season_no_mask.mat')
load('DD_SM_sample_count_1990_now_2_5.mat')

cd('F:\ESA_CCI\broken_stick')
load('SM_loss_regime_02.mat')
load('SM_loss_regime_slope_02.mat')
load('SM_loss_regime_theta_start_02.mat')


mask = DD_SM_sample_count_1990_now_2_5 < 10000 ; 


xmap =  SM_loss_regime_theta_start;
xmap(mask) = NaN ; 
%xmap = categorical(xmap) ; 


figure('units','centimeters','position',[10 3 39 20])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
%h = pcolor(lons_5 - 5, lats_5 + 5, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
%colormap((batlow_color)) 
%colormap(jet(3)) 
%colormap(bam_color(1:128,:))
colormap(bam_color(129:end,:))
%caxis([-0.5 0])  
%caxis([-0.005 0.005])  
caxis([0 0.4])
hcb=colorbar;
%ylabel(hcb,'linear SM slope in stage I  [-]','FontSize',17)
ylabel(hcb,'\theta* (stage I to II transition)  [m³/m³]','FontSize',17)
xlabel('longitude')
ylabel('latitude')
%title('growing season - non growing season')
plot(CoastlineLon, CoastlineLat,'Color','k');
set(gca,'FontSize',17)

print('-image','F:\projects\SM_long_term_DDs\figures\02\theta_star_02','-dbmp','-r150')
close


%% trends in ESA CCI SM raw .. from yearly averages
clear
cd('F:\ESA_CCI\analysis_datasets')
load('ESA_CCI_median_yearly_2_5.mat')
% 1978:2020 is the extent of years. cut to 1990
ESA_CCI_median_yearly_2_5 = ESA_CCI_median_yearly_2_5(:,:,13:43) ; 

cd('F:\ESA_CCI\analysis_datasets\analysis_02')
load('DD_SM_sample_count_1990_now_2_5.mat')
load('dSM_dt_2_5_grow_season_no_mask.mat')

mask = DD_SM_sample_count_1990_now_2_5 < 5000 |  isnan(DD_SM_sample_count_1990_now_2_5) ; 



ESA_CCI_Sm_TS_slope_yearly_2_5 = NaN(72,144) ; 
ESA_CCI_Sm_MK_p_yearly_2_5 = NaN(72,144) ; 

% r was at 44 and c was at 144

for r =1:72
    for c = 1:144
  
  dummy = squeeze(ESA_CCI_median_yearly_2_5(r,c,13:end)) ; 
  cd('E:\Daten Baur\Matlab code')
  [m_cur, b_cur] = TheilSen([(1:length(dummy))' , dummy]) ; 
   ESA_CCI_Sm_TS_slope_yearly_2_5(r,c) = m_cur ;    
  [H,p_value] = Mann_Kendall(dummy,0.05)  ;   
  ESA_CCI_Sm_MK_p_yearly_2_5(r,c) = p_value ; 
    end
    r
end


 ESA_CCI_Sm_MK_p_yearly_2_5(mask) = NaN ; 
 ESA_CCI_Sm_TS_slope_yearly_2_5(mask) = NaN ; 



xmap = ESA_CCI_Sm_TS_slope_yearly_2_5 ;
figure('units','centimeters','position',[10 3 39 20])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
%colormap(redblue_color) 
colormap(bam_color)
clim([-3*10^-3  3*10^-3])
hcb=colorbar;
ylabel(hcb,'SM trend [m³/m³/yr]','FontSize',17)
xlabel('longitude')
ylabel('latitude')
title('ESA CCI SM 1990-2020 trend')
plot(CoastlineLon, CoastlineLat,'Color','k');
hold on
 [rowfind colfind] = find(ESA_CCI_Sm_MK_p_yearly_2_5 < 0.05) ; 
 % crosses1 = plot(lons_2_5(1,colfind),lats_2_5(rowfind,1),'.k','MarkerSize',5) ; 
set(gca,'FontSize',17)
% axes('units','centimeters','Position',[ 7, 5, 4,  4])
% box on
% histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
% set(gca, 'FontSize',12)
% xlim([-3*10^-3  3*10^-3])
print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\ESA_CCI_SM_trend_yearly_1990_2020_02','-dbmp','-r150')
close



%% dSM/dt TS slope based on annual mean loss functions ans 2.5° 
% masks for 10k total observations per pixel and 100 per loss function bin


clear

ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 
% 4080 is 1st Jan 1990

cd('F:\ESA_CCI\analysis_datasets\spatial_02')
dSM_dt_files = string(ls('*dSM_dt_interpsm*')) ;
t_start_end_files = string(ls('*t_start*')) ;
rowcol_files = string(ls('*rowcol*')) ;

cd('F:\ESA_CCI\analysis_datasets\analysis_02')
load('DD_SM_sample_count_1990_now_2_5.mat')

mask_mean = DD_SM_sample_count_1990_now_2_5 < 10000 | isnan(DD_SM_sample_count_1990_now_2_5) ; 

cd('F:\AVHRR_phenology')
load('Phenology_2D_EOS.mat')
load('Phenology_2D_SOS.mat')


% find spatial locations of pos neg 
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_pos')
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_neg')
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_all')
dSM_dt_mean_NDVI_pos_mask = dSM_dt_mean_NDVI_pos > 0 ; 


dSM_dt_2_5_grow_season = NaN(72,144,60) ; 
dSM_dt_2_5_non_grow_season = NaN(72,144,60) ; 
dSM_dt_TS_slope_all = NaN(72,144,1) ;
dSM_dt_mean_all = NaN(72,144,31) ;
dSM_dt_TS_MK_trend_p = NaN(72,144,1) ;
i = 1990:2020 ; 

% spatial = 2450
for spatial = 1:size(dSM_dt_files,1)

cd('F:\ESA_CCI\analysis_datasets\spatial_02')   
dummy_dSM_dt = matfile(dSM_dt_files(spatial)) ; 
dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy ; 
dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% exclude dainage
dummy_dSM_dt(:,40:end) = NaN ; 


dummy_rowcol = matfile(rowcol_files(spatial)) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(t_start_end_files(spatial)) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 


% pre filter for time if needed
mask = dummy_t_mean < 4080 ; % 1-Jan-1990
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 


% get 2.5 degree rowcol from name
name = char(dSM_dt_files(spatial)) ; 
% grep row and col
pat = digitsPattern;
row_col_name = extract(name,pat) ; 
row_2_5 = str2double(row_col_name{1}) ; 
col_2_5 = str2double(row_col_name{2}) ;     
   
% if (dSM_dt_mean_NDVI_pos_mask(row_2_5,col_2_5) > 0) 
%     continue
% end


% indexing in growing and non growing season
EOS_dummy = squeeze(Phenology_2D_EOS(row_2_5,col_2_5,:)) ; 
SOS_dummy = squeeze(Phenology_2D_SOS(row_2_5,col_2_5,:)) ; 

if (all(isnan(EOS_dummy)) || all(isnan(SOS_dummy)))
  continue  
end

% build logical indices based on histcounts                
cd('E:\Daten Baur\Matlab code')       
[left, right] = IntervalUnion(SOS_dummy, EOS_dummy); % FEX
edges = [left,right]';
edges = edges(:);
[~,~,loc] = histcounts(dummy_t_mean,edges);
L = mod(loc,2)==1;          
not_L = ~L ;                 

% get growing season
dummy_dSM_dt_gseason =    dummy_dSM_dt(L,:) ;             
dummy_dSM_dt_non_gseason =    dummy_dSM_dt(not_L,:) ;  
dummy_dSM_dt_all = dummy_dSM_dt ; 
% get timing
dummy_t_mean_gseason = dummy_t_mean(L,:) ; 
dummy_t_mean_non_gseason = dummy_t_mean(not_L,:) ; 
dummy_t_mean_all = dummy_t_mean ; 


% get all data


 % get average annual and then trend   
 dSM_dt_year_dummy = NaN(31,60) ; 
% year_vector = 1990:2:2019 ; 
 year_vector = 1990:2020 ; 
 
for i =  1:length(year_vector)% 1990:2:2019
    
    year_cur = year_vector(i) ;
   % years_overlap = find(year(ESA_CCI_datetime) == year_cur | year(ESA_CCI_datetime) == (year_cur + 1) ) ; 
    years_overlap = find(year(ESA_CCI_datetime) == year_cur) ; 
   
    % gseason
    % select_year = dummy_t_mean_gseason > min(years_overlap) & dummy_t_mean_gseason < max(years_overlap) ; 
    % dummy_select = dummy_dSM_dt_gseason(select_year,:) ; 
    % dummy_select(:,sum(~isnan(dummy_select)) < 100) = NaN ; 

    % non gseason
    % select_year = dummy_t_mean_non_gseason > min(years_overlap) & dummy_t_mean_non_gseason < max(years_overlap) ; 
    % dummy_select = dummy_dSM_dt_non_gseason(select_year,:) ; 
    % dummy_select(:,sum(~isnan(dummy_select)) < 100) = NaN ; 

    % all
     select_year = dummy_t_mean_all > min(years_overlap) & dummy_t_mean_all < max(years_overlap) ; 
     dummy_select = dummy_dSM_dt_all(select_year,:) ; 
     dummy_select(:,sum(~isnan(dummy_select)) < 100) = NaN ; 

    %imagesc(dummy_select)
    dSM_dt_year_dummy(i,:) = mean( dummy_select,1,'omitnan') ; 

    
end

% imagesc(dSM_dt_year_dummy)
dSM_dt_year_dummy(dSM_dt_year_dummy == 0) = NaN ; 

dSM_dt_year_dummy = mean(dSM_dt_year_dummy,2,'omitnan') ; 
cd('E:\Daten Baur\Matlab code')
% [m_cur b_cur] = TheilSen([(1:15)' , dSM_dt_year_dummy]) ; 
  [m_cur b_cur] = TheilSen([(1:31)' , dSM_dt_year_dummy]) ; 

   [H,p_value] = Mann_Kendall(dSM_dt_year_dummy,0.05)  ;   
   dSM_dt_TS_MK_trend_p(row_2_5,col_2_5) = p_value ; 

% figure ; plot(year_vector,dSM_dt_year_dummy,'-o') ; hold on ; 
% plot(year_vector,m_cur*(1:31) + b_cur ) ; yyaxis right ;plot(year_vector,squeeze(VIP_NDVI_yearly_1990_2020_2_5(row_2_5,col_2_5,:)),'-og' )


dSM_dt_mean_all(row_2_5,col_2_5,:) = dSM_dt_year_dummy ; 
dSM_dt_TS_slope_all(row_2_5,col_2_5,:) = m_cur ; 


   spatial
   
end





%%


mask_spatial = DD_SM_sample_count_1990_now_2_5 < 5000 |  isnan(DD_SM_sample_count_1990_now_2_5) ; 
dSM_dt_TS_slope_all(mask_spatial) = NaN ; 
dSM_dt_TS_MK_trend_p(mask_spatial) = NaN ; 

xmap = dSM_dt_TS_slope_all ;
figure('units','centimeters','position',[10 3 39 20])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
%colormap(redblue_color) 
colormap(bam_color)
clim([-8*10^-4  8*10^-4])
hcb=colorbar;
ylabel(hcb,'dSM/dt trend [m³/m³/day/yr]','FontSize',17)
% title('non growing season','FontSize',17)
xlabel('longitude')
ylabel('latitude')
%title('ESA CCI SM 1990-2020 trend')
plot(CoastlineLon, CoastlineLat,'Color','k');
hold on
[rowfind colfind] = find(dSM_dt_TS_MK_trend_p < 0.05) ; 
% crosses1 = plot(lons_2_5(1,colfind),lats_2_5(rowfind,1),'.k','MarkerSize',7) ; 
set(gca,'FontSize',17)


print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSM_dt_TS_slope_1990_2010_02','-dbmp','-r150')
close


imagesc(mean(dSM_dt_1990_2020_TS_slope,3,'omitnan'))



%% boxplots of median loss rate
% 
% 
% clear
% cd('F:\ESA_CCI\2D')
% ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
% tt = 1:15402 ; 
% 
% load('dSM_dt_interpsm_2D_array.mat')
% load('t_start_end_2D_array.mat')
% load('rowcol_2D_array.mat')
% 
% 
% 
% 
% for i = 1:60
% 
% 
% Post_cur = boxplot(dSM_dt_interpsm_2D_array(:,i)  ,'Labels',num2str(sminterp(i)), ...
%     'Whisker', 1.5, 'Jitter', 0.0001, 'Positions',sminterp(i),'Width',0.005) ;
% boxes1 = findobj(gcf,'tag','Box','-and','DisplayName','') ; 
% set(Post_cur(7,:),'Visible','off')
% set(findobj(gcf,'tag','Box','-and','DisplayName',''), 'Color', [0.25 0.25 0.25],'DisplayName','xx','LineWidth',1.3);
% set(findobj(gcf,'tag','Upper Whisker','-and','DisplayName',''), 'Color', [0 0 0],'DisplayName','xx');
% set(findobj(gcf,'tag','Lower Whisker','-and','DisplayName',''), 'Color', [0 0 0],'DisplayName','xx');
% h_post = findobj(gca,'Tag','Box','-and','DisplayName','');
% set(h_post,'DisplayName','MB')
% hold on
%  for j=1:length(boxes1)
%    patches1 =  patch(get(boxes1(j),'XData'),get(boxes1(j),'YData'),col_L,'FaceAlpha',.5);
%  end
%  hold on
% 
% Post_cur = boxplot(dSM_dt_interpsm_2D_array(:,i)  ,'Labels',num2str(sminterp(i)), ...
%     'Whisker', 1.5, 'Jitter', 0.0001, 'Positions',sminterp(i),'Width',0.005) ;
% set(Post_cur(7,:),'Visible','off')
% set(findobj(gcf,'tag','Box','-and','DisplayName',''), 'Color', [0.25 0.25 0.25],'DisplayName','xx','LineWidth',1.3);
% set(findobj(gcf,'tag','Upper Whisker','-and','DisplayName',''), 'Color', [0 0 0],'DisplayName','xx');
% set(findobj(gcf,'tag','Lower Whisker','-and','DisplayName',''), 'Color', [0 0 0],'DisplayName','xx');
% h_post = findobj(gca,'Tag','Box','-and','DisplayName','');
% set(h_post,'DisplayName','MB')
% 
% 
% end
% 


%% correlate SM loss rate and NDVI over time

clear


ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 
% 4080 is 1st Jan 1990

cd('F:\ESA_CCI\analysis_datasets\spatial_02')
dSM_dt_files = string(ls('*dSM_dt_interpsm*')) ;
t_start_end_files = string(ls('*t_start*')) ;
rowcol_files = string(ls('*rowcol*')) ;


cd('F:\ESA_CCI\analysis_datasets\NDVI_sminterp_noanomaly')
VIP_NDVI_files = string(ls('*VIP_NDVI*')) ;


cd('F:\ESA_CCI\analysis_datasets\analysis_02')
load('DD_SM_sample_count_1990_now_2_5.mat')

mask_mean = DD_SM_sample_count_1990_now_2_5 < 6000 | isnan(DD_SM_sample_count_1990_now_2_5) ; 


dSM_dt_2_5_grow_season = NaN(72,144,60) ; 
dSM_dt_2_5_non_grow_season = NaN(72,144,60) ; 
dSM_dt_TS_slope_all = NaN(72,144,1) ;
dSM_dt_mean_all = NaN(72,144,25) ;
dSM_dt_TS_MK_trend_p = NaN(72,144,1) ;
i = 1990:2020 ; 

dSM_dt_NDVI_corr_all = NaN(72,144,1) ;
dSM_dt_NDVI_corr_all_p = NaN(72,144,1) ;


% spatial = 2450
for spatial = 1:size(dSM_dt_files,1)

cd('F:\ESA_CCI\analysis_datasets\spatial_02')   
dummy_dSM_dt = matfile(dSM_dt_files(spatial)) ; 
dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy ; 
% dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;   
% exclude dainage
dummy_dSM_dt(:,40:end) = NaN ; 


% convert to mm/day
% dummy_dSM_dt = dummy_dSM_dt .* -100 ; 


% NDVI 
dummy_NDVI = matfile(strcat('F:\ESA_CCI\analysis_datasets\NDVI_sminterp_noanomaly\',VIP_NDVI_files(spatial))) ; 
dummy_NDVI = dummy_NDVI.NDVI_spatial_sminterp ; 

dummy_rowcol = matfile(rowcol_files(spatial)) ; 
dummy_rowcol = dummy_rowcol.rowcol_dummy ;     
    
dummy_t_start_end = matfile(t_start_end_files(spatial)) ; 
dummy_t_start_end = dummy_t_start_end.t_start_end_dummy ;    
dummy_t_mean = round(mean(dummy_t_start_end,2,'omitnan')) ; 


% pre filter for time if needed. Shall we allow al data until 2020? Or also
% cut 2014 cause no NDVI? Probably cut

 mask = dummy_t_mean < 4080 & dummy_t_start_end(:,1) > 13210; % 1-Jan-1990 until end 2014
 % mask = dummy_t_mean < 4080 ; % 1-Jan-1990 until end 2014
dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_NDVI(mask,:) = NaN ; 


% get 2.5 degree rowcol from name
name = char(dSM_dt_files(spatial)) ; 
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



% convert to positive mm/day loss
dummy_dSM_dt = dummy_dSM_dt .* -100 ; 



 % get average annual and then trend   
 dSM_dt_year_dummy = NaN(25,60) ; 
 NDVI_year_dummy = NaN(25,60) ;  
% year_vector = 1990:2:2019 ; 
 year_vector = 1990:2014 ; 
 
for i =  1:length(year_vector)% 1990:2:2019
    
    year_cur = year_vector(i) ;
   % years_overlap = find(year(ESA_CCI_datetime) == year_cur | year(ESA_CCI_datetime) == (year_cur + 1) ) ; 
    years_overlap = find(year(ESA_CCI_datetime) == year_cur) ; 
   
     select_year = dummy_t_mean > min(years_overlap) & dummy_t_mean < max(years_overlap) ; 
     dummy_select = dummy_dSM_dt(select_year,:) ; 
     dummy_select(:,sum(~isnan(dummy_select)) < 100) = NaN ; 

     dummy_select_NDVI = dummy_NDVI(select_year,:) ; 
     dummy_select_NDVI(:,sum(~isnan(dummy_select)) < 100) = NaN ; 


    %imagesc(dummy_select)
    dSM_dt_year_dummy(i,:) = median( dummy_select,1,'omitnan') ; 
    NDVI_year_dummy(i,:) =   median( dummy_select_NDVI,1,'omitnan') ;     

    
end


% figure
% for i = 1:20
% 
%     plot(dSM_dt_year_dummy(i,:),'o-')
%     hold on
%     pause(0.5)
% 
% end



% imagesc(dSM_dt_year_dummy)
dSM_dt_year_dummy(dSM_dt_year_dummy == 0) = NaN ; 
NDVI_year_dummy  (dSM_dt_year_dummy == 0) = NaN ; 
length_SM_loss_function = sum(~isnan(dSM_dt_year_dummy),2) ; 
all_valid_SM_loss_function = all(~isnan(dSM_dt_year_dummy),1) ; 


% check new condition, remove SM loss functions shortert than 10. Arbitrary
% but should get rid of mini loss functions that are annomalous.
% dSM_dt_year_dummy(length_SM_loss_function < 10,:) = NaN ; 


dSM_dt_year_dummy = median(dSM_dt_year_dummy,2,'omitnan') ; 
NDVI_year_dummy = median(NDVI_year_dummy,2,'omitnan') ; 



cd('E:\Daten Baur\Matlab code')
% correlation between loss rates and NDVI could do correlation across SM
% conditions and then mean or distirbution?? 
  [m_cur, b_cur] = TheilSen([(1:25)' , dSM_dt_year_dummy]) ; 



mask_2 = isnan(NDVI_year_dummy) | isnan(dSM_dt_year_dummy) ; 

dSM_dt_year_dummy_short = dSM_dt_year_dummy ; 
dSM_dt_year_dummy_short(mask_2) = [] ; 
NDVI_year_dummy(mask_2) = [] ; 

 if all(isnan(dSM_dt_year_dummy))
     continue
 end

  if length(NDVI_year_dummy) > 5
  % [R, P] = corr(NDVI_year_dummy , dSM_dt_year_dummy_short) ; 
          % all
        [R,P] = corrcoef([NDVI_year_dummy dSM_dt_year_dummy_short])    ;  
        dSM_dt_NDVI_corr_all(row_2_5,col_2_5,:) =  R(1,2); 
        dSM_dt_NDVI_corr_all_p(row_2_5,col_2_5,:) =  P(1,2) ;   
  end




   [H,p_value] = Mann_Kendall(dSM_dt_year_dummy,0.05)  ;   
   dSM_dt_TS_MK_trend_p(row_2_5,col_2_5) = p_value ; 

% figure ; plot(year_vector,dSM_dt_year_dummy,'-o') ; hold on ; 
% plot(year_vector,m_cur*(1:31) + b_cur ) ; yyaxis right ;plot(year_vector,squeeze(VIP_NDVI_yearly_1990_2020_2_5(row_2_5,col_2_5,:)),'-og' )


dSM_dt_mean_all(row_2_5,col_2_5,:) = dSM_dt_year_dummy ; 
dSM_dt_TS_slope_all(row_2_5,col_2_5,:) = m_cur ; 



   spatial
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dSM_dt_slope_1990_2014 = dSM_dt_TS_slope_all ; 
dSM_dt_MK_trend_p_1990_2014 = dSM_dt_TS_MK_trend_p ; 

save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_slope_1990_2014','dSM_dt_slope_1990_2014') ; 
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_MK_trend_p_1990_2014','dSM_dt_MK_trend_p_1990_2014') ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for test = 1:25
% 
%     plot(dSM_dt_year_dummy(test,:),'-o','Color',bam_color(test.*8,:))
%     hold on
% 
% end






mask_spatial = DD_SM_sample_count_1990_now_2_5 < 3000 |  isnan(DD_SM_sample_count_1990_now_2_5) ; 
% xmap(mask_spatial) = NaN ; 
% xmap(mask_spatial) = NaN ; 



xmap = dSM_dt_TS_slope_all ;
xmap(mask_spatial) = NaN ; 
xmap(mask_spatial) = NaN ; 
figure('units','centimeters','position',[10 3 39 20])  ;
h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
set(h,'LineStyle','none')
shading flat
hold on
%colormap(redblue_color) 
colormap(bam_color)
   clim([-0.1  0.1])
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


print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSM_dt_TS_slope_1990_2010_02','-dbmp','-r150')
close













%% do standard tileplot, based on xmap positive and negativ epixels
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


% find spatial locations of pos neg 
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_pos')
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_neg')
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_all')

load('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_mean_NDVI_FYslope.mat') 

% dSM_dt_mean_NDVI_pos_mask = dSM_dt_mean_NDVI_pos > 0 ; 




ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 

NDVI_binning = linspace(-0.20,0.20,41) ; 
dSM_dt_binning = linspace(-0.08,0.08,41) ; 

sxity_to_thirty = 0:2:60 ; sxity_to_thirty(1) = 1 ; 

counter = NaN(40,40) ; 
counter(:,:) = 1 ; 


% binning based on NDVI
dSM_dt_NDVI_binning = NaN(40,40,1000000) ; 
%dSM_dt_NDVI_binning = NaN(40,40,2000000) ; 


lons_2_5 = (-180+1.25):2.5:(180-1.25)   ; 
lons_2_5 = repmat(lons_2_5,[72, 1]) ; 
lats_2_5 = fliplr((-90+1.25):2.5:(90-1.25))   ; 
lats_2_5 = repmat(lats_2_5',[1, 144]) ; 



% spatial = 2450
% spatial = 2550
cd('F:\ESA_CCI\analysis_datasets')

for spatial = 1:length(dSM_dt_files)


% get drydown data

dummy_dSM_dt = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',dSM_dt_files(spatial))) ; 
NDVI_dummy =  matfile(strcat('F:\ESA_CCI\analysis_datasets\NDVI_sminterp\',NDVI_files_5y(spatial))) ; 
% rowcol_dummy = matfile(strcat('F:\ESA_CCI\analysis_datasets\spatial_02\',rowcol_files(spatial))) ; 


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



% logical condition to only analyze 2,5 degree boxes with extraction of
% cover effect
 if (dSM_dt_mean_NDVI_FYslope(row_2_5,col_2_5) > 0)
     continue
 end




dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy   ; 
% dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;

% less than 3k drydowns remove
if size(dummy_dSM_dt,1) < 3000
    continue
end



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

% pre filter for time if needed pre 1990 and post 2014
mask = dummy_t_start_end(:,2) < 4080 |  dummy_t_start_end(:,1) > 13210; % 1-Jan-1990


dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_NDVI(mask,:) = NaN ; 

% convert to mm/days loss
dummy_dSM_dt = dummy_dSM_dt .* -100  ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporal condition to 1 year
% time_start = find(ESA_CCI_datetime == datetime('01-Jan-2003')) ; 
% time_end = find(ESA_CCI_datetime == datetime('01-Jan-2004')) ; 
% % time_start = find(ESA_CCI_datetime == year_vector_start) ; 
% % time_end = find(ESA_CCI_datetime == year_vector_end) ; 



% mask = dummy_t_mean < time_start |  dummy_t_mean > time_end ; % 
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
dummy_dSM_dt_mean = median(dummy_dSM_dt,1,'omitnan') ; 

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




save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_NDVI_binning_median_extraction','dSM_dt_NDVI_binning_median_extraction')
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\dSM_dt_NDVI_binning_median_extraction_samples','dSM_dt_NDVI_binning_median_extraction_samples')



dSM_dt_NDVI_binning_median_extraction = median(dSM_dt_NDVI_binning,3,'omitnan') ; 
dSM_dt_NDVI_binning_median_extraction_samples = samples_3D_count ; 

dSM_dt_NDVI_binning(dSM_dt_NDVI_binning == 0) = NaN ; 
samples_3D_count = sum(~isnan(dSM_dt_NDVI_binning),3,'omitnan') ;  



xmap = median(dSM_dt_NDVI_binning,3,'omitnan') ; 
 xmap(samples_3D_count < 1000) = NaN ; 
% xmap(samples_3D_count < 100) = NaN ; 

% imagesc(samples_3D_count)


Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
set(h1,'LineStyle','none')
shading flat
clim([-0.5 0.5]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(redblue_color) %
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat anomaly [m³/m³/day]";
cbr.Label.FontSize = 17 ; 
xticks([0:0.2:0.6])
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
title('Only 2.5 degree boxes with cover effect')
pbaspect([0.8 0.8 0.8])
set(gca,'FontSize',17)


print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSM_dt_anomaly_cover_pix_only','-dbmp','-r150')
close



%% do trend analysis for extraction and cover effect only



clear

NDVI_binning = linspace(-0.20,0.20,41) ; 
% dSM_dt_binning = linspace(-0.08,0.08,41) ; 


% cd('F:\ESA_CCI\analysis_datasets\xmap_dSm_dt_02')
cd('F:\ESA_CCI\analysis_datasets\xmap_dSm_dt_03')


filenames = string(ls('*dSM_dt_dev*')) ; 



xmap_array = NaN(40,40,10) ; 

for i = 1:24

    dummy = matfile(filenames(i)) ; 
    dummy_samples = dummy.samples_3D_count ;     
    dummy = dummy.xmap ;
    dummy(dummy_samples < 10) = NaN ; 
    xmap_array(:,:,i) = dummy ; 

i
end

% just plot average?


%





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

SM_loss_mm_d_ESACCI_trend = Xmap_space_slope ;
SM_loss_mm_d_ESACCI_trend_p =  Xmap_space_slope_p ; 

% h1 = pcolor(1:40,1:40,Xmap_space_slope) ;
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\SM_loss_mm_d_ESACCI_trend','SM_loss_mm_d_ESACCI_trend')
save('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\SM_loss_mm_d_ESACCI_trend_p','SM_loss_mm_d_ESACCI_trend_p')



Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = imagesc(sminterp(1:40) ,(NDVI_binning(2:end) - diff(NDVI_binning(1:end))./2)  ,Xmap_space_slope) ;
hold on
% set all nan as white
set(h1, 'AlphaData', ~isnan(Xmap_space_slope))
set(gca,'YDir','normal') 
clim([-0.05 0.05]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.6])
colormap(redblue_color) %
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
 points1 = plot(sminterp(colfind),NDVI_binning(rowfind+1) - 0.0050 ,'.k','MarkerSize',5) ; 

print('-image','F:\projects\SM_long_term_DDs\figures\04_5y_anomaly\dSM_dt_anomaly_NDVI_tile_trend_extraction','-dbmp','-r150')
close





Fig_Panel = figure('units','centimeters','position',[10 2 30 21])  ;
h1 = imagesc(sminterp(1:40) ,(NDVI_binning(2:end) - diff(NDVI_binning(1:end))./2)  ,xmap_array(:,:,1)) ;
hold on
% set all nan as white
set(h1, 'AlphaData', ~isnan(xmap_array(:,:,1)))
set(gca,'YDir','normal') 
clim([-0.05 0.05]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.6])
% colormap(redblue_color) %
cbr = colorbar ;
cbr.Label.String = "\DeltaSM/\Deltat anomaly trend [m³/m³/days/year]";
xticks([0:0.2:0.6])
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
title('1990-2014 Theil Sen slope')
set(gca,'FontSize',17)
pbaspect([1 1 1])





%% redo all xmap SM and NDVI bins for each year. Then use to calculate a trend through the 40x40 surface



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


% find spatial locations of pos neg 
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_pos')
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_neg')
load('F:\ESA_CCI\analysis_datasets\xmap_spatial_means\dSM_dt_mean_NDVI_all')

load('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_mean_NDVI_FYslope.mat') 

% dSM_dt_mean_NDVI_pos_mask = dSM_dt_mean_NDVI_pos > 0 ; 




ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 

NDVI_binning = linspace(-0.20,0.20,41) ; 
dSM_dt_binning = linspace(-0.08,0.08,41) ; 

sxity_to_thirty = 0:2:60 ; sxity_to_thirty(1) = 1 ; 

counter = NaN(40,40) ; 
counter(:,:) = 1 ; 


% ESA_CCI_datetime(4080)
% ESA_CCI_datetime(end)
years_vec = 1990:2014 ; 
for i = 1:length(years_vec)
    datetime_dummy_years(i) = datetime(strcat('01-Jan-',num2str(years_vec(i)))) ;
end
% datetime_dummy_years(end) = ESA_CCI_datetime(end) ; 
 datetime_dummy_years(end) = datetime('31-Dec-2014') ; 



% spatial = 2450
% spatial = 2550
cd('F:\ESA_CCI\analysis_datasets')





for years_index = 10:length(datetime_dummy_years)


counter = NaN(40,40) ; 
counter(:,:) = 1 ; 
% binning based on NDVI
% NDVI_dSM_dt_binning = NaN(60,60,500000) ; 


year_vector_start = datetime_dummy_years(years_index) ; 
year_vector_end = datetime_dummy_years(years_index+1) ; 


% binning based on NDVI
dSM_dt_NDVI_binning = NaN(40,40,300000) ; 





% spatial = 2450
% spatial = 2550



for spatial = 1:length(dSM_dt_files)


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


% logical condition to only analyze 2,5 degree boxes with extraction of
% cover effect
 % if (dSM_dt_mean_NDVI_FYslope(row_2_5,col_2_5) > 0)
 %     continue
 % end




dummy_dSM_dt = dummy_dSM_dt.dSM_dt_interpsm_dummy   ; 
% dummy_dSM_dt(dummy_dSM_dt < -0.4) = NaN ;

% less than 3k drydowns remove
if size(dummy_dSM_dt,1) < 3000
    continue
end



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

% pre filter for time if needed pre 1990 and post 2014
mask = dummy_t_start_end(:,2) < 4080 |  dummy_t_start_end(:,1) > 13210; % 1-Jan-1990


dummy_t_mean(mask) = NaN ; 
dummy_dSM_dt(mask,:) = NaN ; 
dummy_rowcol(mask,:) = NaN ; 
dummy_NDVI(mask,:) = NaN ; 

% convert to mm/days loss
dummy_dSM_dt = dummy_dSM_dt .* -100  ; 


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



% binning based on NDVI
% dSM_dt_NDVI_binning = NaN(30,60) ; 
%NDVI_mean = mean(dummy_NDVI,'omitnan') ; 
NDVI_deviations = dummy_NDVI  ; 
% get mean but for each SM bin
dummy_dSM_dt_mean = median(dummy_dSM_dt,1,'omitnan') ; 

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


    if counter(bin,j) == 300000
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
% xmap(samples_3D_count < 100) = NaN ; 
  % figure
  % pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2,  (xmap)) ; colormap(redblue_color) ; clim([-0.5 0.5])  

save(strcat('F:\ESA_CCI\analysis_datasets\xmap_dSm_dt_03\dSM_dt_dev_',num2str(year(year_vector_start))),'xmap','samples_3D_count')

year_vector_start



end







