%% MJB  15-5-2024 Script to do all plots for ESA CCI content Egg paper.

clear




%% load default plotting variables

set(0, 'DefaultAxesFontSize',16)
% set(0,'defaultlinelinewidth',2)
set(0, 'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesTitleFontWeight','normal')



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
sminterp_zeppe = (linspace(0,1,50))' ; 


NDVI_binning = linspace(-0.20,0.20,41) ; 
T2M_binning  = linspace(-10,10,41) ; 

NDVI_binning_Zeppe = linspace(-0.2,0.2,41) ; 
T2M_binning_Zeppe  = linspace(-10,10,41) ; 





%% Fig. 1 Theoretic loss curve and estimated loss curves maybe for regions
% do pre analysis here
clear
cd('F:\ESA_CCI\2D_02')
ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:15402 ; 

load('dSM_dt_interpsm_2D_array.mat')
load('t_start_end_2D_array.mat')
load('rowcol_2D_array.mat')



% exclude all data before 1990
dSM_dt_interpsm_2D_array(t_start_end_2D_array(:,1) < 4080,:) = NaN ; 

[Lia, finddate] = find(ESA_CCI_datetime == datetime('31-Dec-2014')) ; 

% exclude all data after 2014 .. NDVI only goes until then
dSM_dt_interpsm_2D_array(t_start_end_2D_array(:,2) > finddate,:) = NaN ; 

% convert to mm/day
dSM_dt_interpsm_2D_array = dSM_dt_interpsm_2D_array .* -100 ; 

dSM_dt_median_1990_2014 =   median(dSM_dt_interpsm_2D_array,1,'omitnan') ;
dSM_dt_mean_1990_2014 =   mean(dSM_dt_interpsm_2D_array,1,'omitnan') ;
dSM_dt_25prct_1990_2014 =   prctile(dSM_dt_interpsm_2D_array,25,1) ; 
dSM_dt_75prct_1990_2014 =   prctile(dSM_dt_interpsm_2D_array,75,1) ; 
samples_across_SM = sum(~isnan(dSM_dt_interpsm_2D_array),1,'omitnan') ; 
clear dSM_dt_interpsm_2D_array 

% add zeppe model

% load('F:\projects\SM_long_term_DDs\Zeppe_model_output_06_mm_d\DD_dSM_dt_shallow_interpsm.mat')
load('F:\projects\SM_long_term_DDs\Zeppe_model_output_08\DD_dSM_dt_shallow_interpsm.mat')

% convert to mm/d loss
DD_dSM_dt_shallow_interpsm = DD_dSM_dt_shallow_interpsm .* -100 ; 
dSM_dt_median_zeppe =   median(DD_dSM_dt_shallow_interpsm,1,'omitnan') ;
dSM_dt_mean_zeppe =   mean(DD_dSM_dt_shallow_interpsm,1,'omitnan') ;
dSM_dt_25prct_zeppe =   prctile(DD_dSM_dt_shallow_interpsm,25,1) ; 
dSM_dt_75prct_zeppe =   prctile(DD_dSM_dt_shallow_interpsm,75,1) ; 
clear DD_dSM_dt_shallow_interpsm
pre_index_zeppe = 2:49 ; 

% pre_index_zeppe = 2:59 ; 
% test_boxes = vertcat(dSM_dt_25prct_1990_2014,dSM_dt_median_1990_2014,dSM_dt_75prct_1990_2014) ; 
% test_boxes(:,samples_across_SM < 1000) = NaN ; 
% boxplot(test_boxes,'Whisker',Inf)


% here theoretical SM loss curve
SM = 0:0.001:0.6 ; 
loss_rate = NaN(1,601) ; 
% stage 2
loss_rate(1,1:40) = 0 ; 
loss_rate(1,41:250) = linspace(0,1.01,210) ; 
% stage 1
loss_rate(1,251:400) = 1.01 ; 
% drainage
loss_rate(1,401:601) = 1 + exp(linspace(0,5,201)) ./ 100 ; 

pre_index = 2:60 ; 
prct_25 = dSM_dt_25prct_1990_2014(pre_index) ;
prct_50 = dSM_dt_median_1990_2014(pre_index) ;
prct_75 = dSM_dt_75prct_1990_2014(pre_index) ;


% save('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_median_1990_2014','dSM_dt_median_1990_2014')
% save('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_mean_1990_2014','dSM_dt_mean_1990_2014')
% save('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_25prct_1990_2014','dSM_dt_25prct_1990_2014')
% save('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_75prct_1990_2014','dSM_dt_75prct_1990_2014')
% save('F:\projects\SM_long_term_DDs\data_for_figures\samples_across_SM','samples_across_SM')
% save('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_median_zeppe','dSM_dt_median_zeppe')




load('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_median_1990_2014')
load('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_mean_1990_2014')
load('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_25prct_1990_2014')
load('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_75prct_1990_2014')
load('F:\projects\SM_long_term_DDs\data_for_figures\samples_across_SM')
load('F:\projects\SM_long_term_DDs\data_for_figures\dSM_dt_median_zeppe')



dSM_dt_median_1990_2014 = dSM_dt_median_1990_2014 .* -100 ;
dSM_dt_mean_1990_2014 = dSM_dt_mean_1990_2014 .* -100 ;
dSM_dt_25prct_1990_2014 = dSM_dt_25prct_1990_2014 .* -100 ;
dSM_dt_75prct_1990_2014 = dSM_dt_75prct_1990_2014 .* -100 ;




%% start figure
Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

%----------------- Panel 1 theoretical loss function ----------------
sub1 = subplot(2,2,1) ; 
plot(SM,loss_rate,'k','LineWidth',1.5)
xlabel('SM [m³/m³]','FontSize',16)
ylabel('SM loss [mm/day]','FontSize',16)
yticklabels('')
xticklabels('')
yline(1.01,'k--','LineWidth',1)
xline(SM(40),'k--','LineWidth',1)
xline(SM(251),'k--','LineWidth',1)
xline(SM(401),'k--','LineWidth',1)
text(SM(45),2.25,'\theta_{\it w}','FontSize',16)
text(SM(256),2.25,'\theta_*','FontSize',16)
text(SM(406),2.25,'\theta_{\it fc}','FontSize',16)
fontsize(16,'points')


%----------------- Panel 2 global loss function ----------------
sub2 = subplot(2,2,2) ; 
linepre = plot(sminterp(pre_index),prct_50,'-','LineWidth',1.5,'Color','k') ;
hold on
%
x2 = [sminterp(pre_index)', fliplr(sminterp(pre_index)')];
inBetween = [prct_50, fliplr(prct_25 )];
fillpre = fill(x2, inBetween, col_L,'FaceAlpha',0.25);
x2 = [sminterp(pre_index)', fliplr(sminterp(pre_index)')];
inBetween = [prct_50, fliplr(prct_75)];
fill(x2, inBetween, col_L,'FaceAlpha',0.25);
plot(sminterp(pre_index),prct_50,'-','LineWidth',1.5,'Color','k')
% yticks(-0.3:0.075:0)
% ylim([-0.4 0])
% xlim([0 0.6])
line25 = plot(sminterp(pre_index),prct_25,'-','LineWidth',1.5,'Color','k') ;
line75 = plot(sminterp(pre_index),prct_75,'-','LineWidth',1.5,'Color','k') ;
%kolsmi_post = plot(sminterp(pre_index),prct_50,'','MarkerSize',15) ;
xlabel('SM [m³/m³]','FontSize',16)
ylabel('SM loss [mm/day]','FontSize',16)
%title('Soil moisture loss function 1990-2020')
% add indicators
drainage_thresh = xline(0.437,'k--','LineWidth',1) ; 
drrainage_y = yline(5,'k--','LineWidth',1) ;
text(0.45,26,'\theta_{\it fc}','FontSize',16)
% annotation('textbox',[.29 .17 .27 .06],'String','Stage II and stage II SM loss','FontSize',16)
% annotation('textbox',[.76 .17 .1 .06],'String','Drainage','FontSize',22)
set(gca,'FontSize',16)

%----------------- Panel Add on model loss function ----------------

% linepre = plot(sminterp(pre_index_zeppe),dSM_dt_median_zeppe(pre_index_zeppe)./2,'-','LineWidth',0.5,'Color','k') ;
% hold on
% %
% x2 = [sminterp(pre_index_zeppe)', fliplr(sminterp(pre_index_zeppe)')];
% inBetween = [dSM_dt_median_zeppe(pre_index_zeppe), fliplr(dSM_dt_25prct_zeppe(pre_index_zeppe) )];
% fillpre = fill(x2, inBetween, col_C,'FaceAlpha',0.25);
% x2 = [sminterp(pre_index_zeppe)', fliplr(sminterp(pre_index_zeppe)')];
% inBetween = [dSM_dt_median_zeppe(pre_index_zeppe), fliplr(dSM_dt_75prct_zeppe(pre_index_zeppe))];
% fill(x2, inBetween, col_C,'FaceAlpha',0.25);
% model_median_ref = plot(sminterp(pre_index_zeppe),dSM_dt_median_zeppe(pre_index_zeppe),'-','LineWidth',1.5,'Color','r') ; 

model_median_ref = plot(sminterp(1:44),dSM_dt_median_zeppe(1:44),'-','LineWidth',1.5,'Color','r') ; 
% yticks(-0.3:0.075:0)
% ylim([-0.4 0])
% xlim([0 0.6])
% line25 = plot(sminterp(pre_index_zeppe),dSM_dt_25prct_zeppe(pre_index_zeppe),'-','LineWidth',1.5,'Color','k') ;
% line75 = plot(sminterp(pre_index_zeppe),dSM_dt_75prct_zeppe(pre_index_zeppe),'-','LineWidth',1.5,'Color','k') ;
% 
legend('ESA CCI','','','','','','','','median model','Location','northwest')



%----------------- Panel 3 map of average SM loss ----------------
% thought is to show and quantify differences in stage 1 +2 loss .. maybe
% make dependent on NDVI, so it fits with the next analysis steps. Could
% then later compare with model? 
% load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\DD_dSM_dt_mean_1990_now_2014.mat')
% load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\DD_SM_sample_count_1990_2014_2_5.mat')
load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\DD_dSM_dt_median_1990_now_2014_nodrain.mat')
load('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\DD_SM_sample_count_1990_2014_2_5_nodrain.mat')


DD_dSM_dt_mean_1990_now_2014_nodrain(lats_2_5 > 70) = NaN ; 

xmap = median(DD_dSM_dt_mean_1990_now_2014_nodrain,3,'omitnan');
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% MJB comment: !! set 3000 total samples as threshold for this study. Gets
% rid of coastal pixels etc. Shall we carry this forward for binning plots
% as well? Probably smart to do

sub3 = subplot(2,2,3) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat
h = imagesc(lons_2_5(1,:), lats_2_5(:,1), xmap) ; 
set(h, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
% colormap(sub3,flipud(bam_color(1:128,:))) 
 colormap(sub3,flipud(batlow_color)) 
% colormap(sub3,(batlow_color)) 
clim([0 3])% dSM/dt diurnal
xticks(-160:80:160)
yticks(-80:40:80)
hcb2=colorbar;
set(hcb2, 'FontSize',16)
set(hcb2, 'FontSize',16,'YTick',0:1:4)
xlabel('longitude','FontSize',16)
ylabel('latitude','FontSize',16)
ylabel(hcb2,' median SM loss [mm/day]','FontSize',16)
%ylabel(hcb2,'\DeltaSM/\Deltat change [%]','FontSize',16)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
% inset with pdf
histoaxes = axes('units','centimeters','Position',[ 4.5, 16.5-8-4.25+1, 3,  3]) ; 
box on
histo2 = histogram(xmap,50,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(histoaxes, 'FontSize',12)
xlim([0 3])

% set positions relative to FIG = [10 2 40 24]

set(sub1,'units','centimeter','position',  [ 2.5,2+11.5+3, 11, 7])
set(sub2,'units','centimeter','position',  [ 2.5+12+2.5, 2+11.5+3, 11, 7])
set(sub3,'units','centimeter','position',  [ 2.5, 1.7, 11+11+3.5,(11+11+3.5) / 2])
% see if this works for everything


textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 2.5,2+11.6+3+7,   1,1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [ 2.5+12+2.5, 2+11.6+3+7, 1, 1], 'EdgeColor', 'none')
textbox3_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'c)', 'FontSize', 20) ; 
set(textbox3_label,'Units','centimeters', 'Position', [ 2.5, 1.7+(11+11+3.6) / 2, 1,1], 'EdgeColor', 'none')

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',2,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[ 2.5+11/4,2+13.4, 5.5, 0]) ; 


fontsize(16,'points')
set(histoaxes, 'FontSize',12)



 set(Fig_Panel, 'PaperUnits', 'inches');        % use inches for consistency
 set(Fig_Panel, 'PaperSize', [15 11]);     % A4 size
% set(Fig_Panel, 'PaperPosition', [0 0 8.27 11.69]);  % fill the page


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_1_model_median_04','svg')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_1_model_median_04','pdf')

% set(Fig_Panel, 'PaperPositionMode', 'auto');  % ensures the figure fills the page
% print(Fig_Panel, 'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_1_model_median_04.pdf', '-dpdf', '-bestfit');


close 







%% Fig Panel 2 Global effect of NDVI anomalies on SM loss rate and map of average slope relativ eto NDVI
% prepare data
clear

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision')

load('dSM_dt_NDVI_binning_median.mat')
load('dSM_dt_NDVI_binning_median_samples.mat')
load('DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
load('dSM_dt_mean_NDVI_FYslope.mat')




%% start figure




%----------------- subplot 1 SM loss function relative to NDVI phase space global ----------------

xmap = dSM_dt_NDVI_binning_median ; 
xmap(dSM_dt_NDVI_binning_median_samples < 1000) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(1,2,1) ;
% h1 = pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat

h1 = imagesc(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

clim([-0.5 0.5]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-0.25:0.125:0.25)
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
pbaspect([1 1 1])


%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------


xmap = -dSM_dt_mean_NDVI_FYslope ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% xmap(xmap > 0) = NaN ; 

sub2 = subplot(1,2,2) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat
h1 = imagesc(lons_2_5(1,:), lats_2_5(:,1), xmap) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(redblue_color) 
clim([-0.1  0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaNDVI [mm/day]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes = axes('units','centimeters','Position',[ 3+8+4+3.2, 7.5, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.1  0.1])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow1,'Position',[15+20+2 (6+8)-3   0  3]) ; 
 set(arrow2,'Position',[15+20+2 (6+8)-5  0 -3]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+20+2.6 (6+8)-5   3  0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+20+2.6 (6+8)-3  -3 0]) ; 
set(gca,'Box','on');

textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+0.1, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+6,6+8+0.1, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,6, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+4, 6, 20, 8])

fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)




saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures\Panel_2','svg')
close 




%% -------------------   Panel 3 --------------------------
% display areas with negative slope relative to NDVI and show that
% extraction effect is limited by available soil moisture. Do same as panel
% 2 basically


clear

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d')

load('dSM_dt_NDVI_binning_median_extraction_samples.mat')
load('dSM_dt_NDVI_binning_median_extraction.mat')

load('DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
load('dSM_dt_mean_NDVI_FYslope.mat')






%----------------- Subplot 1 SM loss function relative to NDVI phase space global ----------------

xmap = dSM_dt_NDVI_binning_median_extraction ; 
xmap(dSM_dt_NDVI_binning_median_extraction_samples < 1000) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(1,2,1) ;
% h1 = pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat

h1 = imagesc(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

clim([-0.5 0.5]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-0.25:0.125:0.25)
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
pbaspect([1 1 1])
fontsize(16,'points')


%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------


xmap = -dSM_dt_mean_NDVI_FYslope ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
xmap(xmap < 0) = NaN ; 

sub2 = subplot(1,2,2) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat

h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(sub2,redblue_color(51:100,:)) 
clim([0  0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaNDVI [mm/day]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes = axes('units','centimeters','Position',[ 3+8+4+3.2, 7.5, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([0  0.1])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map


 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow2,'Position',[15+20+2 (6)  0 8]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+20+2.6 (6+8)-2   3  0]) ; 
set(gca,'Box','on');
 


textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+0.1, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+6,6+8+0.1, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,6, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+4, 6, 20, 8])

fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)




saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures\Panel_3','svg')
close 




%% -------------------   Panel 4 --------------------------
% Same as standard SM loss but for T2M 


clear

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d')

load('dSM_dt_T2M_binning_median.mat')
 load('dSM_dt_T2M_binning_median_samples.mat')
load('dSM_dt_mean_T2M_FYslope.mat')
load('DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
load('DD_drydown_sampling_zeppe.mat')


% load SM constant NDVI and T2; 
load('dSM_dt_NDVI_T2M_anomaly_Zeppe_05_2D.mat')
load('dSM_dt_NDVI_T2M_anomaly_Zeppe_05_2D_sampling.mat')




%% start figure




%----------------- subplot 1 SM loss function relative to NDVI phase space global ----------------

xmap = -dSM_dt_T2M_binning_median ; 
xmap(dSM_dt_T2M_binning_median_samples < 1000) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(1,2,1) ;
% h1 = pcolor(sminterp(1:40),T2M_binning(2:end) - diff(T2M_binning(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat

h1 = imagesc(sminterp(1:40),T2M_binning(2:end) - diff(T2M_binning(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

clim([-0.5 0.5]) %
 ylim([-10 10])
 xlim([-0.00 0.5])
colormap(redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-10:5:10)
xlabel('SM [m³/m³]')
ylabel('T2M anomaly [°K]')
pbaspect([1 1 1])
fontsize(16,'points')


%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------


xmap = -dSM_dt_mean_T2M_FYslope ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% xmap(xmap > 0) = NaN ; 

sub2 = subplot(1,2,2) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat

h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(redblue_color) 
clim([-0.1  0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaT2M [mm/day/°K]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes = axes('units','centimeters','Position',[ 3+8+4+3.2, 7.5, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.1  0.1])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow1,'Position',[15+20+2 (6+8)-3   0  3]) ; 
 set(arrow2,'Position',[15+20+2 (6+8)-5  0 -3]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+20+2.6 (6+8)-5   3  0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+20+2.6 (6+8)-3  -3 0]) ; 
set(gca,'Box','on');

textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+0.1, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+6,6+8+0.1, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,6, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+4, 6, 20, 8])

fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)

saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures\Panel_4','svg')
close 








%% ------------------ Panel 5  -----------------------------------
% make model phase space and map. shall we use 2.5 or 5 degree map? and
% what about anomalies vs seasonality? Anomalies probably better
clear

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d')


load('dSM_dt_NDVI_anomaly_2D_zeppe_ERA.mat')
load('dSM_dt_NDVI_anomaly_2D_zeppe_ERA_sampling.mat')

load('DD_drydown_sampling_zeppe_ERA.mat')
load('DD_dSM_dt_NDVI_gradient_zeppe_ERA.mat')

load('DD_SM_sample_count_1990_2014_2_5_nodrain')



sminterp_zeppe = linspace(0,0.4,40)' ; 
%----------------- subplot 1 SM loss function relative to NDVI phase space global ----------------

xmap = dSM_dt_NDVI_anomaly_2D_zeppe_ERA ; 
 xmap(dSM_dt_NDVI_anomaly_2D_zeppe_ERA_sampling < 1000) = NaN ; 
xmap(:,sminterp_zeppe > 0.4) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(1,2,1) ;
% h1 = pcolor(sminterp_zeppe,NDVI_binning_Zeppe(2:end) - diff(NDVI_binning_Zeppe(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat

h1 = imagesc(sminterp(1:40),NDVI_binning_Zeppe(2:end) - diff(NDVI_binning_Zeppe(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

clim([-0.5 0.5]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-0.25:0.125:0.25)
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
pbaspect([1 1 1])
fontsize(16,'points')


%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = DD_dSM_dt_NDVI_gradient_zeppe_ERA ;
% xmap(DD_drydown_sampling_zeppe_05 < 3000) = NaN ; 
 xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% xmap(xmap > 0) = NaN ; 

sub2 = subplot(1,2,2) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat

h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(redblue_color) 
clim([-0.1 0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaNDVI [mm/day]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes = axes('units','centimeters','Position',[ 3+8+4+3.2, 7.5, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.1  0.1])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow1,'Position',[15+20+2 (6+8)-3   0  3]) ; 
 set(arrow2,'Position',[15+20+2 (6+8)-5  0 -3]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+20+2.6 (6+8)-5   3  0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+20+2.6 (6+8)-3  -3 0]) ; 
set(gca,'Box','on');

textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+0.1, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+6,6+8+0.1, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,6, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+4, 6, 20, 8])

fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)



saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures\Panel_5_model','svg')
close 





%% ------------------------ Panel 6 --------------------------------------
% repeat Panel 1 but all based on Model results. Maps of median Sm loss
% rate as well as sampling and map
cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision')

load('DD_dSM_dt_loss_rate_2_5_zeppe_ERA.mat')
load('DD_drydown_sampling_zeppe_ERA.mat')
load('F:\projects\SM_long_term_DDs\Zeppe_model_output_06_mm_d\DD_dSM_dt_shallow_interpsm.mat')


DD_dSM_dt_shallow_interpsm = DD_dSM_dt_shallow_interpsm .* -100 ; 
% DD_dSM_dt_shallow_interpsm = DD_dSM_dt_shallow_interpsm .*2 ; 
dSM_dt_median_zeppe =   median(DD_dSM_dt_shallow_interpsm,1,'omitnan') ;
dSM_dt_mean_zeppe =   mean(DD_dSM_dt_shallow_interpsm,1,'omitnan') ;
dSM_dt_25prct_zeppe =   prctile(DD_dSM_dt_shallow_interpsm,25,1) ; 
dSM_dt_75prct_zeppe =   prctile(DD_dSM_dt_shallow_interpsm,75,1) ; 

sminterp_zeppe = linspace(0,0.5,50)' ; 



% here theoretical SM loss curve
SM = 0:0.001:0.6 ; 
loss_rate = NaN(1,601) ; 
% stage 2
loss_rate(1,1:40) = 0 ; 
loss_rate(1,41:250) = linspace(0,1.01,210) ; 
% stage 1
loss_rate(1,251:400) = 1.01 ; 
% drainage
loss_rate(1,401:601) = 1+ exp(linspace(0,5,201)) ./ 100 ; 




Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

%----------------- Panel 1 theoretical loss function ----------------
sub1 = subplot(2,2,1) ; 
plot(SM,loss_rate,'k','LineWidth',1.5)
xlabel('SM [m³/m³]','FontSize',16)
ylabel('SM loss [mm/day]','FontSize',16)
yticklabels('')
xticklabels('')
yline(1.01,'k--','LineWidth',1)
xline(SM(40),'k--','LineWidth',1)
xline(SM(251),'k--','LineWidth',1)
xline(SM(401),'k--','LineWidth',1)
text(SM(45),2.25,'\theta_{\it w}','FontSize',16)
text(SM(256),2.25,'\theta_*','FontSize',16)
text(SM(406),2.25,'\theta_{\it fc}','FontSize',16)
fontsize(16,'points')

pre_index = 2:49 ; 
%----------------- Panel 2 global loss function ----------------
sub2 = subplot(2,2,2) ; 
linepre = plot(sminterp_zeppe(pre_index),dSM_dt_median_zeppe(pre_index),'-','LineWidth',1.5,'Color','k') ;
hold on
%
x2 = [sminterp_zeppe(pre_index)', fliplr(sminterp_zeppe(pre_index)')];
inBetween = [dSM_dt_median_zeppe(pre_index), fliplr(dSM_dt_25prct_zeppe(pre_index) )];
fillpre = fill(x2, inBetween, col_L,'FaceAlpha',0.25);
x2 = [sminterp_zeppe(pre_index)', fliplr(sminterp_zeppe(pre_index)')];
inBetween = [dSM_dt_median_zeppe(pre_index), fliplr(dSM_dt_75prct_zeppe(pre_index))];
fill(x2, inBetween, col_L,'FaceAlpha',0.25);
plot(sminterp_zeppe(pre_index),dSM_dt_median_zeppe(pre_index),'-','LineWidth',1.5,'Color','k')
% yticks(-0.3:0.075:0)
 ylim([0 30])
 xlim([0 0.6])
line25 = plot(sminterp_zeppe(pre_index),dSM_dt_25prct_zeppe(pre_index),'-','LineWidth',1.5,'Color','k') ;
line75 = plot(sminterp_zeppe(pre_index),dSM_dt_75prct_zeppe(pre_index),'-','LineWidth',1.5,'Color','k') ;
%kolsmi_post = plot(sminterp(pre_index),prct_50,'','MarkerSize',15) ;
xlabel('SM [m³/m³]','FontSize',16)
ylabel('SM loss model [mm/day]','FontSize',16)
%title('Soil moisture loss function 1990-2020')
% add indicators
% drainage_thresh = xline(0.437,'k--','LineWidth',1) ; 
% drrainage_y = yline(-0.05,'k--','LineWidth',1) ;
% text(0.45,-0.26,'\theta_{\it fc}','FontSize',16)
% annotation('textbox',[.29 .17 .27 .06],'String','Stage II and stage II SM loss','FontSize',16)
% annotation('textbox',[.76 .17 .1 .06],'String','Drainage','FontSize',22)
set(gca,'FontSize',16)


%----------------- Panel 3 map of average SM loss ----------------
% thought is to show and quantify differences in stage 1 +2 loss .. maybe
% make dependent on NDVI, so it fits with the next analysis steps. Could
% then later compare with model? 


xmap = median(DD_dSM_dt_loss_rate_2_5_zeppe_ERA,3,'omitnan') ;
% xmap(DD_drydown_sampling_zeppe_ERA < 3000) = NaN ; 
% MJB comment: !! set 3000 total samples as threshold for this study. Gets
% rid of coastal pixels etc. Shall we carry this forward for binning plots
% as well? Probably smart to do

sub3 = subplot(2,2,3) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat

h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(sub3,flipud(bam_color(1:128,:))) 
clim([0 5])% dSM/dt diurnal
xticks(-160:80:160)
yticks(-80:40:80)
hcb2=colorbar;
set(hcb2, 'FontSize',16)
set(hcb2, 'FontSize',16,'YTick',0:1.25:5)
xlabel('longitude','FontSize',16)
ylabel('latitude','FontSize',16)
ylabel(hcb2,' median SM loss model [mm/day]','FontSize',16)
%ylabel(hcb2,'\DeltaSM/\Deltat change [%]','FontSize',16)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
% inset with pdf
histoaxes = axes('units','centimeters','Position',[ 4.5, 16.5-8-4.25+1, 3,  3]) ; 
box on
histo2 = histogram(xmap,50,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(histoaxes, 'FontSize',12)
xlim([0 5])

% set positions relative to FIG = [10 2 40 24]

set(sub1,'units','centimeter','position',  [ 2.5,2+11.5+3, 11, 7])
set(sub2,'units','centimeter','position',  [ 2.5+12+2.5, 2+11.5+3, 11, 7])
set(sub3,'units','centimeter','position',  [ 2.5, 1.7, 11+11+3.5,(11+11+3.5) / 2])
% see if this works for everything


textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 2.5,2+11.6+3+7,   1,1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [ 2.5+12+2.5, 2+11.6+3+7, 1, 1], 'EdgeColor', 'none')
textbox3_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'c)', 'FontSize', 20) ; 
set(textbox3_label,'Units','centimeters', 'Position', [ 2.5, 1.7+(11+11+3.6) / 2, 1,1], 'EdgeColor', 'none')

fontsize(16,'points')
set(histoaxes, 'FontSize',12)


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures\Panel_1_model','svg')
close 



%% -----------------------   Panel 7 --------------------------------------
% dSM/dt relative to T2M global phase space and map
clear

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d\')

load('DD_SM_sample_count_1990_2014_2_5_nodrain.mat')

load('dSM_dt_T2M_anomaly_2D_zeppe_ERA_sampling.mat')
load('dSM_dt_T2M_anomaly_2D_zeppe_ERA.mat')
load('DD_dSM_dt_T2M_gradient_zeppe_ERA.mat')



sminterp_zeppe = linspace(0,0.6,60)' ; 

%----------------- subplot 1 SM loss function relative to NDVI phase space global ----------------

xmap = dSM_dt_T2M_anomaly_2D_zeppe_ERA ; 
xmap(dSM_dt_T2M_anomaly_2D_zeppe_ERA_sampling < 1000) = NaN ; 
% xmap(:,sminterp_zeppe > 0.4) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(2,2,1) ;
% h1 = pcolor(sminterp_zeppe(1:40),T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
h1 = imagesc(sminterp(1:40),T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 


clim([-0.5 0.5]) %
 ylim([-10 10])
 xlim([-0.00 0.5])
colormap(redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-10:5:10)
xlabel('SM [m³/m³]')
ylabel('T2M anomaly [°K]')
pbaspect([1 1 1])
fontsize(16,'points')


% ----------------  subplot 2 Add T2M and NDVI phase space with fixed SM ------------








%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = DD_dSM_dt_T2M_gradient_zeppe_ERA ;
  xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% xmap(xmap > 0) = NaN ; 

sub2 = subplot(2,2,3) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat

h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(redblue_color) 
clim([-0.1  0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaT2M [m³/m³/day/°K]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes = axes('units','centimeters','Position',[ 3+8+4+3.2, 7.5, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.1  0.1])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow1,'Position',[15+20+2 (6+8)-3   0  3]) ; 
 set(arrow2,'Position',[15+20+2 (6+8)-5  0 -3]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+20+2.6 (6+8)-5   3  0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+20+2.6 (6+8)-3  -3 0]) ; 
set(gca,'Box','on');

textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+0.1, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+6,6+8+0.1, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,6, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+4, 6, 20, 8])

fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures\Panel_8_fixed_SM_model','svg')
close 





%% -----------------------   Panel 8 --------------------------------------
% dSM/dt relative to T2M global phase space and map
clear


cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\')

load('DD_SM_sample_count_1990_2014_2_5_nodrain.mat')

load('dSM_dt_T2M_anomaly_2D_zeppe_ERA_sampling.mat')
load('dSM_dt_T2M_anomaly_2D_zeppe_ERA.mat')
load('DD_dSM_dt_T2M_gradient_zeppe_ERA.mat')

load('dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA.mat')
load('dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA_sampling.mat')


sminterp_zeppe = linspace(0,0.6,60)' ; 

%----------------- subplot 1 SM loss function relative to NDVI phase space global ----------------

xmap = dSM_dt_T2M_anomaly_2D_zeppe_ERA ; 
xmap(dSM_dt_T2M_anomaly_2D_zeppe_ERA_sampling < 1000) = NaN ; 
% xmap(:,sminterp_zeppe > 0.4) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(2,2,1) ;
% h1 = pcolor(sminterp_zeppe(1:40),T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
h1 = imagesc(sminterp(1:40),T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

clim([-0.5 0.5]) %
 ylim([-10 10])
 xlim([-0.00 0.5])
colormap(redblue_color) %
% hcb = colorbar ;
% hcb.Label.String = "SM loss anomaly [mm/day]";
% set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-10:5:10)
xlabel('SM [m³/m³]')
ylabel('T2M anomaly [°K]')
pbaspect([1 1 1])
fontsize(16,'points')
yline(0,'--','LineWidth',0.5)


% ----------------  subplot 2 Add T2M and NDVI phase space with fixed SM ------------


sub2 = subplot(2,2,2) ;
xmap = dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA ; 
xmap(dSM_dt_NDVI_T2M_anomaly_2D_zeppe_ERA_sampling  < 1000) = NaN ; 
% h2 = pcolor(T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2,   NDVI_binning_Zeppe(2:end) - diff(NDVI_binning_Zeppe(1:2))/2, (xmap)) ;
% set(h2,'LineStyle','none')
% shading flat
h1 = imagesc(T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2,NDVI_binning_Zeppe(2:end) - diff(NDVI_binning_Zeppe(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 


clim([-0.5 0.5]) %
% clim([-6e-4 0.0]) %
 ylim([-0.2 0.2])
 xlim([-10 10])
 colormap(redblue_color) %
% colormap parula
 cbr1 = colorbar ;
 cbr1.Label.String = "Sm loss anomaly [mm/day]";
 cbr1.Units = 'centimeters' ; 
 set(cbr1, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xlabel('T2M anomaly [°K]')
ylabel('NDVI anomaly [-]')
% title('Europe')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% fontsize(16,'points')
yline(0,'--','LineWidth',0.5)
xline(0,'--','LineWidth',0.5)


%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = DD_dSM_dt_T2M_gradient_zeppe_ERA ;
  xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
 xmap(lats_2_5 > 70) = NaN ; 

sub3 = subplot(2,2,3) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat

h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(redblue_color) 
clim([-0.15  0.15])
hcb=colorbar('XTick',-0.15:0.075:0.15);
ylabel(hcb,'\DeltaSM loss / \DeltaT2M [mm day^{-1} K^{-1}]','FontSize',18)
hcb.Units = 'centimeters' ; 
cbr1_pos = cbr1.Position ; 
hcb.Position =[23.0893, 2.5, 0.5644,(10+8+2) / 2];       
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes = axes('units','centimeters','Position',[  4, 4.5, 2.5 ,2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.2  0.2])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow1,'Position',[ 2.5+12+8+4.2, 2+11+2+5,   0,  3]) ; 
 set(arrow2,'Position',[ 2.5+12+8+4.2, 2+11+2+3,  0, -3]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[2.5+12+8+4.7, 2+11+2+3,   3,  0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[2.5+12+8+4.7, 2+11+2+5, -3, 0]) ; 
set(gca,'Box','on');




 arrow3 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow4 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow3,'Position',[ 2.5+12+8+4.2, 2+1.5+5,   0,  4]) ; 
 set(arrow4,'Position',[ 2.5+12+8+4.2, 2+1.5+3,  0, -4]) ; 
 
textbox3 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox3,'Position',[2.5+12+8+4.7, 2+1.5+3,   3,  0]) ; 
set(gca,'Box','on');
 
textbox4 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox4,'Position',[2.5+12+8+4.7, 2+1.5+5, -3, 0]) ; 
set(gca,'Box','on')




textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 2.5,2+11+2+8, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [ 2.5+12, 2+11+2+8, 1, 1], 'EdgeColor', 'none')
textbox3_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'c)', 'FontSize', 20) ; 
set(textbox3_label,'Units','centimeters', 'Position', [ 2.5, 2.5+10, 1,1], 'EdgeColor', 'none')



set(sub1,'units','centimeter','position',  [ 2.5,2+11+2, 8, 8])
set(sub2,'units','centimeter','position',  [ 2.5+12, 2+11+2, 8, 8])
set(sub3,'units','centimeter','position',  [ 2.5, 2.5, 10+8+2,(10+8+2) / 2])



fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_4_T2M_fixed_SM_model','svg')
saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_4_T2M_fixed_SM_model','pdf')

exportgraphics(Fig_Panel, 'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_4_T2M_fixed_SM_model.pdf', ...
    'ContentType', 'vector');

close 




%% -----------------------   Panel 9 --------------------------------------
% dSM/dt relative to T2M global phase space and map and SM held constant
% for observation


clear


cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision\')

load('DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
load('dSM_dt_mean_T2M_FYslope.mat')

load('dSM_dt_T2M_binning_median.mat')
load('dSM_dt_T2M_binning_median_samples.mat')

load('dSM_dt_NDVI_T2M_anomaly_2D.mat')
load('dSM_dt_NDVI_T2M_anomaly_2D_sampling.mat')

sminterp_zeppe = linspace(0,0.6,60)' ; 

%----------------- subplot 1 SM loss function relative to NDVI phase space global ----------------


xmap = -dSM_dt_T2M_binning_median ; 
xmap(dSM_dt_T2M_binning_median_samples < 1000) = NaN ; 
% xmap(:,sminterp_zeppe > 0.4) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(2,2,1) ;
% h1 = pcolor(sminterp_zeppe(1:40),T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
h1 = imagesc(sminterp(1:40),T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

clim([-0.5 0.5]) %
 ylim([-10 10])
 xlim([-0.00 0.5])
colormap(redblue_color) %
% hcb = colorbar ;
% hcb.Label.String = "SM loss anomaly [mm/day]";
% set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-10:5:10)
xlabel('SM [m³/m³]')
ylabel('T2M anomaly [°K]')
pbaspect([1 1 1])
fontsize(16,'points')
yline(0,'--','LineWidth',0.5)


% ----------------  subplot 2 Add T2M and NDVI phase space with fixed SM ------------


sub2 = subplot(2,2,2) ;
xmap = -dSM_dt_NDVI_T2M_anomaly_2D ; 
xmap(dSM_dt_NDVI_T2M_anomaly_2D_sampling  < 1000) = NaN ; 
% h2 = pcolor(T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2,   NDVI_binning_Zeppe(2:end) - diff(NDVI_binning_Zeppe(1:2))/2, (xmap)) ;
% set(h2,'LineStyle','none')
% shading flat
h1 = imagesc(T2M_binning_Zeppe(2:end) - diff(T2M_binning_Zeppe(1:2))/2,   NDVI_binning_Zeppe(2:end) - diff(NDVI_binning_Zeppe(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 


clim([-0.5 0.5]) %
% clim([-6e-4 0.0]) %
 ylim([-0.2 0.2])
 xlim([-10 10])
 colormap(redblue_color) %
% colormap parula
 cbr1 = colorbar ;
 cbr1.Label.String = "Sm loss anomaly [mm/day]";
 cbr1.Units = 'centimeters' ; 
 set(cbr1, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xlabel('T2M anomaly [°K]')
ylabel('NDVI anomaly [-]')
% title('Europe')
set(gca,'FontSize',17)
pbaspect([1 1 1])
% fontsize(16,'points')
yline(0,'--','LineWidth',0.5)
xline(0,'--','LineWidth',0.5)


%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = -dSM_dt_mean_T2M_FYslope ;
  xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
xmap(lats_2_5 > 70) = NaN ; 

sub3 = subplot(2,2,3) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat

h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(redblue_color) 
clim([-0.15  0.15])
hcb=colorbar('XTick',-0.15:0.075:0.15);
ylabel(hcb,'\DeltaSM loss / \DeltaT2M [mm day^{-1} K^{-1}]','FontSize',18)
hcb.Units = 'centimeters' ; 
cbr1_pos = cbr1.Position ; 
hcb.Position =[23.0893, 2.5, 0.5644,(10+8+2) / 2];
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes = axes('units','centimeters','Position',[  4, 4.5, 2.5 ,2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.2  0.2])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow1,'Position',[ 2.5+12+8+4.2, 2+11+2+5,   0,  3]) ; 
 set(arrow2,'Position',[ 2.5+12+8+4.2, 2+11+2+3,  0, -3]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[2.5+12+8+4.7, 2+11+2+3,   3,  0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[2.5+12+8+4.7, 2+11+2+5, -3, 0]) ; 
set(gca,'Box','on');




 arrow3 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow4 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow3,'Position',[ 2.5+12+8+4.2, 2+1.5+5,   0,  4]) ; 
 set(arrow4,'Position',[ 2.5+12+8+4.2, 2+1.5+3,  0, -4]) ; 
 
textbox3 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox3,'Position',[2.5+12+8+4.7, 2+1.5+3,   3,  0]) ; 
set(gca,'Box','on');
 
textbox4 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox4,'Position',[2.5+12+8+4.7, 2+1.5+5, -3, 0]) ; 
set(gca,'Box','on')




textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 2.5,2+11+2+8, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [ 2.5+12, 2+11+2+8, 1, 1], 'EdgeColor', 'none')
textbox3_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'c)', 'FontSize', 20) ; 
set(textbox3_label,'Units','centimeters', 'Position', [ 2.5, 2.5+10, 1,1], 'EdgeColor', 'none')



set(sub1,'units','centimeter','position',  [ 2.5,2+11+2, 8, 8])
set(sub2,'units','centimeter','position',  [ 2.5+12, 2+11+2, 8, 8])
set(sub3,'units','centimeter','position',  [ 2.5, 2.5, 10+8+2,(10+8+2) / 2])



fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_3_T2M_obs_fixed_SM_model','svg')
close 




%% ------------------------ Panel 10 -----------------------------
% redo basic NDVI phase space for Zeppe model extraction only


clear

cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision')


load('dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction.mat')
load('dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction_sampling.mat')


load('DD_drydown_sampling_zeppe_ERA.mat')
load('DD_dSM_dt_NDVI_gradient_zeppe_ERA.mat')
load('DD_SM_sample_count_1990_2014_2_5_nodrain')



sminterp_zeppe = linspace(0,0.4,40)' ; 
%----------------- subplot 1 SM loss function relative to NDVI phase space global ----------------

xmap = dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction ; 
xmap(dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction_sampling < 1000) = NaN ; 
xmap(:,sminterp_zeppe > 0.4) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(1,2,1) ;
% h1 = pcolor(sminterp_zeppe,NDVI_binning_Zeppe(2:end) - diff(NDVI_binning_Zeppe(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
h1 = imagesc(sminterp(1:40),NDVI_binning_Zeppe(2:end) - diff(NDVI_binning_Zeppe(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

clim([-0.5 0.5]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-0.25:0.125:0.25)
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
pbaspect([1 1 1])
fontsize(16,'points')


%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = DD_dSM_dt_NDVI_gradient_zeppe_ERA ;
% xmap(DD_drydown_sampling_zeppe_05 < 3000) = NaN ; 
 xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
 xmap(xmap < 0) = NaN ; 

sub2 = subplot(1,2,2) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat
h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(sub2,redblue_color(51:100,:)) 
clim([0  0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaNDVI [mm/day]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes = axes('units','centimeters','Position',[ 3+8+4+3.2, 7.5, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([0  0.1])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[15+20+2, 6,   0,  8]) ; 
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+20+2.6, 6+2,  -3, 0]) ; 
set(gca,'Box','on');

textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+0.1, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+6,6+8+0.1, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,6, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+4, 6, 20, 8])

fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)




saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures\Panel_10_SM_loss_NDVI_extract_model','svg')
close 










%% ------------------------ Panel 10 -----------------------------
% Panel 10 a) Do standard cover effect but double for obs and model


clear

% CCI
cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision')

load('dSM_dt_NDVI_binning_median.mat')
load('dSM_dt_NDVI_binning_median_samples.mat')
load('DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
load('dSM_dt_mean_NDVI_FYslope.mat')

% Zeppe
cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision')

load('dSM_dt_NDVI_anomaly_2D_zeppe_ERA.mat')
load('dSM_dt_NDVI_anomaly_2D_zeppe_ERA_sampling.mat')
load('DD_drydown_sampling_zeppe_ERA.mat')
load('DD_dSM_dt_NDVI_gradient_zeppe_ERA.mat')
load('DD_SM_sample_count_1990_2014_2_5_nodrain')






%----------------- Subplot 1 SM loss function relative to NDVI phase space global ----------------



xmap = dSM_dt_NDVI_binning_median ; 
xmap(dSM_dt_NDVI_binning_median_samples < 1000) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(2,2,1) ;
h1 = imagesc(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
hold on
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 
clim([-0.5 0.5]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-0.25:0.125:0.25)
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
pbaspect([1 1 1])
fontsize(16,'points')
xtickangle(0)
yline(0,'--','LineWidth',0.5)

%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------


xmap = -dSM_dt_mean_NDVI_FYslope ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
% xmap(xmap < 0) = NaN ; 
xmap(lats_2_5 > 70) = NaN ; 

sub2 = subplot(2,2,2) ; 
h = imagesc(lons_2_5(1,:) , lats_2_5(:,1) , xmap); 
set(h, 'AlphaData', ~isnan(xmap))
set(sub2,'YDir','normal') 
hold on
colormap(sub2,redblue_color) 
clim([-0.1  0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaNDVI [mm/day]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes1 = axes('units','centimeters','Position',[ 3+8+4+3.2+2, 7.5+8, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.1  0.1])
set(histoaxes1, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

 arrow1 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[15+0.3, 6+13,  0, 3]) ; 

 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow2,'Position',[15+0.3, 6+11,  0, -3]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+0.3+0.7, 6+13,   -3,  0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+0.3+0.7, 6+11,   3,  0]) ; 
set(gca,'Box','on');


textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+0.1+8, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+6+2,6+8+0.1+8, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,6+8, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+4+2, 6+8, 20, 8])




% ------------------- here same but zeppe model -----------------------------


xmap = dSM_dt_NDVI_anomaly_2D_zeppe_ERA ; 
 xmap(dSM_dt_NDVI_anomaly_2D_zeppe_ERA_sampling < 1000) = NaN ; 
% xmap(:,sminterp_zeppe > 0.4) = NaN ; 


sub3 = subplot(2,2,3) ;
h2 = imagesc(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
hold on
set(h2, 'AlphaData', ~isnan(xmap))
set(sub3,'YDir','normal') 
clim([-0.5 0.5]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(sub3,redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-0.25:0.125:0.25)
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
pbaspect([1 1 1])
fontsize(16,'points')
yline(0,'--','LineWidth',0.5)


%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = DD_dSM_dt_NDVI_gradient_zeppe_ERA ;
% xmap(DD_drydown_sampling_zeppe_05 < 3000) = NaN ; 
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
 % xmap(xmap < 0) = NaN ; 
xmap(lats_2_5 > 70) = NaN ; 

sub4 = subplot(2,2,4) ; 
h3 = imagesc(lons_2_5(1,:) , lats_2_5(:,1) , xmap); 
set(h3, 'AlphaData', ~isnan(xmap))
set(sub4,'YDir','normal') 
hold on
colormap(sub4,redblue_color) 
clim([-0.1  0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaNDVI [mm/day]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes2 = axes('units','centimeters','Position',[ 3+8+4+3.2+2, 4.5, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.1  0.1])
set(histoaxes2, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map


 arrow3 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow3,'Position',[15+0.3, 6,   0,  -3]) ; 

 arrow4 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow4,'Position',[15+0.3, 8,   0,  3]) ; 
 

textbox3 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox3,'Position',[15+0.3+0.7, 8,  -3, 0]) ; 
set(gca,'Box','on');

textbox4 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox4,'Position',[15+0.3+0.7, 6,  3, 0]) ; 
set(gca,'Box','on');


textbox3_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'c)', 'FontSize', 20) ; 
set(textbox3_label,'Units','centimeters', 'Position', [ 3,3+8+0.1, 1, 1], 'EdgeColor', 'none')
textbox4_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'd)', 'FontSize', 20) ; 
set(textbox4_label,'Units','centimeters', 'Position', [3+8+6+2,3+8+0.1, 1, 1], 'EdgeColor', 'none')


set(sub3,'units','centimeter','position',  [ 3,3, 8, 8])
set(sub4,'units','centimeter','position',  [ 3+8+4+2, 3, 20, 8])

fontsize(16,'points')
xtickangle(sub3, 0 )
% final only set histoaxes smaller
set(histoaxes1, 'FontSize',12)
set(histoaxes2, 'FontSize',12)




saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_1_SM_loss_NDVI_ObsModel_double','svg')
close 













%% ------------------------ Panel 10 -----------------------------
% Do 4 panel extraction effect limit for both model and observations

clear

% CCI
cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision')

load('dSM_dt_NDVI_binning_median_extraction_samples.mat')
load('dSM_dt_NDVI_binning_median_extraction.mat')

load('DD_SM_sample_count_1990_2014_2_5_nodrain.mat')
load('dSM_dt_mean_NDVI_FYslope.mat')

% Zeppe
cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision')

load('dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction.mat')
load('dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction_sampling.mat')

load('DD_drydown_sampling_zeppe_ERA.mat')
load('DD_dSM_dt_NDVI_gradient_zeppe_ERA.mat')
load('DD_SM_sample_count_1990_2014_2_5_nodrain')









%----------------- Subplot 1 SM loss function relative to NDVI phase space global ----------------


xmap = dSM_dt_NDVI_binning_median_extraction ; 
xmap(dSM_dt_NDVI_binning_median_extraction_samples < 1000) = NaN ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;

sub1 = subplot(2,2,1) ;
% h1 = pcolor(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
h1 = imagesc(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 



clim([-0.5 0.5]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-0.25:0.125:0.25)
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
pbaspect([1 1 1])
fontsize(16,'points')
yline(0,'--','LineWidth',0.5)

%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------


xmap = -dSM_dt_mean_NDVI_FYslope ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
xmap(xmap < 0) = NaN ; 

sub2 = subplot(2,2,2) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat
h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 


hold on
colormap(sub2,redblue_color(51:100,:)) 
clim([0  0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaNDVI [mm/day]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes1 = axes('units','centimeters','Position',[ 3+8+4+3.2+2, 7.5+8, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([0  0.1])
set(histoaxes1, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

  arrow1 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[15+0.3, 6+13,  0, 3]) ; 

 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow2,'Position',[15+0.3, 6+11,  0, -3]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+0.3+0.7, 6+13,   -3,  0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+0.3+0.7, 6+11,   3,  0]) ; 
set(gca,'Box','on');

 
textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+0.1+8, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+6+2,6+8+0.1+8, 1, 1], 'EdgeColor', 'none')


set(sub1,'units','centimeter','position',  [ 3,6+8, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+4+2, 6+8, 20, 8])

fontsize(16,'points')
xtickangle(sub1, 0 )



% ------------------- here same but zeppe model -----------------------------


xmap = dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction ; 
 xmap(dSM_dt_NDVI_anomaly_2D_zeppe_ERA_extraction_sampling < 1000) = NaN ; 
% xmap(:,sminterp_zeppe > 0.4) = NaN ; 
% xmap(:,sminterp > 0.4) = NaN ; 


sub3 = subplot(2,2,3) ;
% h1 = pcolor(sminterp_zeppe,NDVI_binning_Zeppe(2:end) - diff(NDVI_binning_Zeppe(1:2))/2, (xmap)) ;
% set(h1,'LineStyle','none')
% shading flat
h1 = imagesc(sminterp(1:40),NDVI_binning(2:end) - diff(NDVI_binning(1:2))/2, (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

clim([-0.5 0.5]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])
colormap(sub3,redblue_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly [mm/day]";
set(hcb, 'FontSize',16,'YTick',-0.5:0.25:0.5)
xticks(0:0.125:0.5)
yticks(-0.25:0.125:0.25)
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
pbaspect([1 1 1])
fontsize(16,'points')
yline(0,'--','LineWidth',0.5)


%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = DD_dSM_dt_NDVI_gradient_zeppe_ERA ;
% xmap(DD_drydown_sampling_zeppe_05 < 3000) = NaN ; 
 xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
 xmap(xmap < 0) = NaN ; 

sub4 = subplot(2,2,4) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat
h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 

hold on
colormap(sub4,redblue_color(51:100,:)) 
clim([0  0.1])
hcb=colorbar;
ylabel(hcb,'\DeltaSM loss / \DeltaNDVI [mm/day]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on
histoaxes2 = axes('units','centimeters','Position',[ 3+8+4+3.2+2, 4.5, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([0  0.1])
set(histoaxes2, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

  arrow3 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow3,'Position',[15+0.3, 6,   0,  -3]) ; 

 arrow4 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow4,'Position',[15+0.3, 8,   0,  3]) ; 
 

textbox3 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox3,'Position',[15+0.3+0.7, 8,  -3, 0]) ; 
set(gca,'Box','on');

textbox4 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox4,'Position',[15+0.3+0.7, 6,  3, 0]) ; 
set(gca,'Box','on');

textbox3_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'c)', 'FontSize', 20) ; 
set(textbox3_label,'Units','centimeters', 'Position', [ 3,3+8+0.1, 1, 1], 'EdgeColor', 'none')
textbox4_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'd)', 'FontSize', 20) ; 
set(textbox4_label,'Units','centimeters', 'Position', [3+8+6+2,3+8+0.1, 1, 1], 'EdgeColor', 'none')


set(sub3,'units','centimeter','position',  [ 3,3, 8, 8])
set(sub4,'units','centimeter','position',  [ 3+8+4+2, 3, 20, 8])

fontsize(16,'points')
xtickangle(sub3, 0 )

set(histoaxes1, 'FontSize',12)
set(histoaxes2, 'FontSize',12)


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_2_SM_loss_NDVI_extract_double','svg')
close 





%% do trends plot from oservations




clear

% CCI
cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d')

load('SM_loss_mm_d_ESACCI_trend.mat')
load('SM_loss_mm_d_ESACCI_trend_p.mat')
load('dSM_dt_MK_trend_p_1990_2014.mat')
load('dSM_dt_slope_1990_2014.mat')

load('DD_SM_sample_count_1990_2014_2_5_nodrain.mat')



% Model .. probably do the same






%----------------- subplot 1 SM loss function trend in phase space ----------------

xmap = SM_loss_mm_d_ESACCI_trend ; 

Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;
sub1 = subplot(1,2,1) ;
h1 = imagesc(sminterp(1:40) ,(NDVI_binning(2:end) - diff(NDVI_binning(1:end))./2)  ,xmap) ;
hold on
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 
clim([-0.025 0.025]) %
 ylim([-0.25 0.25])
 xlim([-0.00 0.5])

% plot significance
 [rowfind colfind] = find(SM_loss_mm_d_ESACCI_trend_p < 0.05) ; 
 points1 = plot(sminterp(colfind),NDVI_binning(rowfind+1) - 0.0050 ,'.k','MarkerSize',4) ; 

colormap(sub1,bam_color) %
hcb = colorbar ;
hcb.Label.String = "SM loss anomaly trend [mm/day/year]";
set(hcb, 'FontSize',16,'YTick',-0.025:0.0125:0.025)
xticks(0:0.125:0.5)
yticks(-0.25:0.125:0.25)
xlabel('SM [m³/m³]')
ylabel('NDVI anomaly [-]')
pbaspect([1 1 1])
fontsize(16,'points')


%----------------- subplot 2 SM loss function trend to NDVI 2.5 degree map ----------------
xmap = dSM_dt_slope_1990_2014 ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
xmap(xmap < -0.4 | xmap > 0.4) = NaN ; 

sub2 = subplot(1,2,2) ; 
h = imagesc(lons_2_5(1,:), lats_2_5(:,1) , xmap); 
set(h, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 
hold on
colormap(sub2,bam_color) 
clim([-0.1  0.1])
hcb=colorbar;
ylabel(hcb,'SM loss trend [mm/day/year]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');

% plot significance
 [rowfind colfind] = find(dSM_dt_MK_trend_p_1990_2014 < 0.05) ; 
 points1 = plot(lons_2_5(1,colfind),lats_2_5(rowfind,:),'.k','MarkerSize',4) ; 


pbaspect([144 72 1])
hold on
histoaxes = axes('units','centimeters','Position',[ 3+8+4+3.2, 7.5, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.1  0.1])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');
% calc positions
% 15     6    20     8  are the pos of the map

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 arrow2 = annotation('arrow',[0.955 0.955],[0.8  0.6],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 
 set(arrow1,'Position',[15+20+2 (6+8)-3   0  3]) ; 
 set(arrow2,'Position',[15+20+2 (6+8)-5  0 -3]) ; 
 
textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+20+2.6 (6+8)-5   3  0]) ; 
set(gca,'Box','on');
 
textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+20+2.6 (6+8)-3  -3 0]) ; 
set(gca,'Box','on');

textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+0.1, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+6,6+8+0.1, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,6, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+4, 6, 20, 8])

fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures\Panel_12_SM_loss_CCI_trends','svg')
close 



% set renderer
     % set(0, 'DefaultFigureRenderer', 'painters');
     % set(0, 'DefaultFigureRenderer', 'opengl');








%% Panel 13 do  trend figure. Left panels will be line plots of Sm loss trends, right panels trend maps
% this is probably easier than understanding the SM NDVI phase space trend
% figure
% ---------------- for now only do Obs ---------------------------
clear



cd('F:\projects\SM_long_term_DDs\data_for_figures_mm_d_revision')

load('dSM_dt_MK_trend_p_1990_2014.mat')
load('Global_SM_loss_trend_ESCACCI_p.mat')
load('Global_SM_loss_trend_ESACCI.mat')
load('dSM_dt_slope_1990_2014.mat')
load('DD_SM_sample_count_1990_2014_2_5_nodrain')

% model
load('Global_SM_loss_trend_p.mat')
load('Global_SM_loss_trend.mat')
% model
load('dSM_dt_slope_1990_2014_zeppe')
load('dSM_dt_MK_trend_p_1990_2014_zeppe')


sminterp_zeppe = linspace(0,0.4,40)' ; 

%----------------- subplot 1 SM loss function trend line  ----------------

 Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;
% Fig_Panel = figure('units','centimeters','position',[10 2 40 12])  ;


sub1 = subplot(1,2,1) ;
linepre = plot(sminterp_zeppe,Global_SM_loss_trend_ESACCI,'bo-','LineWidth',2.0,'Color',bam_color(240,:)) ;
hold on
 ylim([-0.01 0.05])
xlim([-0.00 0.4])
% xticks(0:0.125:0.5)
yticks(-0.01:0.01:0.05)
xlabel('SM [m³/m³]')
ylabel('SM loss trend [mm day^{-1} year^{-1}]')
pbaspect([1 1 1])
fontsize(16,'points')
plot(sminterp_zeppe(Global_SM_loss_trend_ESCACCI_p < 0.05),Global_SM_loss_trend_ESACCI(:,Global_SM_loss_trend_ESCACCI_p < 0.05),'k.','MarkerSize',15) ;
legend('','p < 0.05')

%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = dSM_dt_slope_1990_2014 ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 
xmap(lats_2_5 > 70) = NaN ; 

sub2 = subplot(1,2,2) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat
h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 
hold on
colormap(sub2,bam_color) 
% clim([-0.2  0.2])
clim([-0.1 0.1])
hcb=colorbar;
ylabel(hcb,'SM loss trend [mm day^{-1} year^{-1}]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on

[rowfind colfind] = find(dSM_dt_MK_trend_p_1990_2014 < 0.05) ; 
crosses1 = plot(sub2,lons_2_5(1,colfind),lats_2_5(rowfind,1),'.k','MarkerSize',4) ; 


histoaxes = axes('units','centimeters','Position',[ 3+8+4+0.2, 7.5+9, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.2  0.2])
% xlim([-0.1  0.1])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[15+18.5, 6+5+9,   0,  3]) ; 

  arrow2 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow2,'Position',[15+18.5, 6+3+9,   0,  -3]) ; 


textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+19+0.1, 6+5+9,  -3, 0]) ; 
set(gca,'Box','on');

textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+19+0.1, 6+3+9,  3, 0]) ; 
set(gca,'Box','on');


textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+9+0.1, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+3,6+8+9+0.1, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,6+9, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+1, 6+9, 20, 8])






% PANELs for model
Global_SM_loss_trend = Global_SM_loss_trend(1:40)' ; 
Global_SM_loss_trend_p = Global_SM_loss_trend_p(1:40)' ; 

% Global_SM_loss_trend = Global_SM_loss_trend.*10 ; 


sub3 = subplot(2,2,3) ;
linepre = plot(sminterp_zeppe,Global_SM_loss_trend,'bo-','LineWidth',2.0,'Color',bam_color(240,:)) ;
hold on
 ylim([-0.01 0.05])
xlim([-0.00 0.4])
% xticks(0:0.125:0.5)
yticks(-0.01:0.01:0.05)
xlabel('SM [m³/m³]')
ylabel('SM loss trend model [mm day^{-1} year^{-1}]')
pbaspect([1 1 1])
fontsize(16,'points')
plot(sminterp_zeppe(Global_SM_loss_trend_p < 0.05),Global_SM_loss_trend(Global_SM_loss_trend_p < 0.05),'k.','MarkerSize',15) ;
legend('','p < 0.05')



% map
xmap = dSM_dt_slope_1990_2014_zeppe ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 

% xmap = xmap.*10 ; 


sub4 = subplot(2,2,4) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat
h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 
hold on
colormap(sub4,bam_color) 
% clim([-0.2  0.2])
clim([-0.1 0.1])
hcb=colorbar;
ylabel(hcb,'SM loss trend model [mm day^{-1} year^{-1}]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on

[rowfind colfind] = find(dSM_dt_MK_trend_p_1990_2014_zeppe < 0.05) ; 
crosses1 = plot(sub4,lons_2_5(1,colfind),lats_2_5(rowfind,1),'.k','MarkerSize',4) ; 


histoaxes2 = axes('units','centimeters','Position',[ 3+8+4+0.2, 7.5-2, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.2  0.2])
% xlim([-0.1  0.1])
set(histoaxes2, 'FontSize',12)

axes_diff_Position = get(sub4, 'Position');

 arrow3 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow3,'Position',[15+18.5, 6+5-2,   0,  3]) ; 

  arrow4 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow4,'Position',[15+18.5, 6+3-2,   0,  -3]) ; 


textbox3 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox3,'Position',[15+19+0.1, 6+5-2,  -3, 0]) ; 
set(gca,'Box','on');

textbox4 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox4,'Position',[15+19+0.1, 6+3-2,  3, 0]) ; 
set(gca,'Box','on');


textbox3_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'c)', 'FontSize', 20) ; 
set(textbox3_label,'Units','centimeters', 'Position', [ 3,6+8-2+0.1, 1, 1], 'EdgeColor', 'none')
textbox4_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'd)', 'FontSize', 20) ; 
set(textbox4_label,'Units','centimeters', 'Position', [3+8+3,6+8-2+0.1, 1, 1], 'EdgeColor', 'none')



set(sub3,'units','centimeter','position',  [ 3,4, 8, 8])
set(sub4,'units','centimeter','position',  [ 3+8+1, 4, 20, 8])


fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)
set(histoaxes2, 'FontSize',12)


saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_figures_02\Panel_13_SM_loss_trends_ESACCI_single','svg')
close 



%%




%----------------- subplot 1 SM loss function trend line  ----------------

 Fig_Panel = figure('units','centimeters','position',[10 2 40 24])  ;


sub1 = subplot(1,2,1) ;
linepre = plot(sminterp_zeppe,Global_SM_loss_trend_ESACCI,'bo-','LineWidth',2.0,'Color',bam_color(240,:)) ;
hold on
 ylim([-0.01 0.05])
xlim([-0.00 0.4])
% xticks(0:0.125:0.5)
yticks(-0.01:0.01:0.05)
xlabel('SM [m³/m³]')
ylabel('SM loss trend [mm day^{-1} year^{-1}]')
pbaspect([1 1 1])
fontsize(16,'points')
plot(sminterp_zeppe(Global_SM_loss_trend_ESCACCI_p < 0.05),Global_SM_loss_trend_ESACCI(:,Global_SM_loss_trend_ESCACCI_p < 0.05),'k.','MarkerSize',15) ;
legend('','p < 0.05')
yline(0,'--','LineWidth',0.5)
%----------------- subplot 2 SM loss function relative to NDVI 2.5 degree map ----------------
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, DD_drydown_sampling_zeppe_03); 

xmap = dSM_dt_slope_1990_2014 ;
xmap(DD_SM_sample_count_1990_2014_2_5_nodrain < 3000) = NaN ; 


sub2 = subplot(1,2,2) ; 
% h = pcolor(lons_2_5 - 1.25, lats_2_5 + 1.25, xmap); 
% set(h,'LineStyle','none')
% shading flat
h1 = imagesc(lons_2_5(1,:),lats_2_5(:,1), (xmap)) ; 
set(h1, 'AlphaData', ~isnan(xmap))
set(gca,'YDir','normal') 
hold on
colormap(sub2,bam_color) 
% clim([-0.2  0.2])
clim([-0.1 0.1])
hcb=colorbar;
ylabel(hcb,'SM loss trend [mm day^{-1} year^{-1}]','FontSize',18)
xlabel('longitude')
ylabel('latitude')
xticks(-160:80:160)
yticks(-80:40:80)
plot(CoastlineLon, CoastlineLat,'Color','k');
pbaspect([144 72 1])
hold on

[rowfind colfind] = find(dSM_dt_MK_trend_p_1990_2014 < 0.05) ; 
crosses1 = plot(sub2,lons_2_5(1,colfind),lats_2_5(rowfind,1),'.k','MarkerSize',4) ; 


histoaxes = axes('units','centimeters','Position',[ 3+8+4+0.2, 7.5+9-12, 2.5,  2.5]) ; 
box on
histo2 = histogram(xmap,'FaceColor',[0.402, 0.402, 0.402],'EdgeColor','none') ;
set(gca, 'FontSize',12)
xlim([-0.2  0.2])
% xlim([-0.1  0.1])
set(histoaxes, 'FontSize',12)

axes_diff_Position = get(sub2, 'Position');

 arrow1 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow1,'Position',[15+18.5, 6+5+9-12,   0,  3]) ; 

  arrow2 = annotation('arrow',[0.955 0.955],[0.7 0.9],'LineWidth',5,'HeadLength',15,'HeadWidth',15,'Units','centimeters') ;
 set(arrow2,'Position',[15+18.5, 6+3+9-12,   0,  -3]) ; 


textbox1 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','faster SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox1,'Position',[15+19+0.1, 6+5+9-12,  -3, 0]) ; 
set(gca,'Box','on');

textbox2 =  annotation('textarrow',[0.5 0.5],[0.5 0.5],'string','slower SM loss' , ...
'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.02 .6 0 0],'FontSize',16,'Units','centimeters');
set(textbox2,'Position',[15+19+0.1, 6+3+9-12,  3, 0]) ; 
set(gca,'Box','on');


textbox1_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'a)', 'FontSize', 20) ; 
set(textbox1_label,'Units','centimeters', 'Position', [ 3,6+8+9+0.1-12, 1, 1], 'EdgeColor', 'none')
textbox2_label = annotation('textbox', [0, 0.2, 0, 0],  'string', 'b)', 'FontSize', 20) ; 
set(textbox2_label,'Units','centimeters', 'Position', [3+8+3,6+8+9+0.1-12, 1, 1], 'EdgeColor', 'none')

set(sub1,'units','centimeter','position',  [ 3,3, 8, 8])
set(sub2,'units','centimeter','position',  [ 3+8+1, 3, 20, 8])


fontsize(16,'points')
xtickangle(sub1, 0 )
set(histoaxes, 'FontSize',12)

saveas(Fig_Panel,'F:\projects\SM_long_term_DDs\figures\zeppe_model_02\Panel_Figures_revision\Panel_13_SM_loss_trends_ESACCI_single','svg')
close 




