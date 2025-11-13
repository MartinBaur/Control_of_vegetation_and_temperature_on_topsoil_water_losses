%% MJB 07.03.2023 extract DDs and seasonally correct from ESA CCI

clear
cd('F:\ESA_CCI\processed')

% get datetime array
ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
tt = 1:14301 ; 
dLength = 4 ; 

%% build data arrys .. maybe prepare for parfor 
load('F:\ESA_CCI\means\SM_CCI_Passive_mean.mat')
load('SM_CCI_Combined_mean.mat')
load('F:\ESA_CCI\CCI_lon.mat')
load('F:\ESA_CCI\CCI_lat.mat')

CCI_lat = repmat(CCI_lat,[1 1440]) ; 
CCI_lon = repmat(CCI_lon',[720 1]) ; 
% filter snow and ice areas
SM_CCI_Passive_mean(CCI_lat > 70) = NaN ; 

%get row and col list for extracting data 
[CCI_rows, CCI_cols] = find(~isnan(SM_CCI_Passive_mean)) ; 
CCI_rowcol_list = [CCI_rows CCI_cols] ; 


[CCI_rows, CCI_cols] = find(~isnan(SM_CCI_Combined_mean)) ; 
CCI_rowcol_list = [CCI_rows CCI_cols] ; 


%% extract data in equal size chunks .. then run parallel drydown on chunks.
% check  memory to find a decent chunksize
cd('F:\ESA_CCI\processed')
SM_CCI_01 = matfile('SM_CCI_1.mat') ; 
SM_CCI_02 = matfile('SM_CCI_2.mat') ; 
SM_CCI_03 = matfile('SM_CCI_3.mat') ; 
SM_CCI_04 = matfile('SM_CCI_4.mat') ; 
SM_CCI_05 = matfile('SM_CCI_5.mat') ; 
SM_CCI_06 = matfile('SM_CCI_6.mat') ; 
SM_CCI_07 = matfile('SM_CCI_7.mat') ; 
SM_CCI_08 = matfile('SM_CCI_8.mat') ; 
SM_CCI_09 = matfile('SM_CCI_9.mat') ; 





total_length = size(CCI_rowcol_list,1) ; 
sequence = (1:5000:total_length)' ; 
sminterp = (0.01:0.01:0.6)' ; 




%% out loop going through spatial chunks
row_chunks = round(linspace(1,720,6)) ; 
col_chunks = round(linspace(1,1440,6)) ; 

% rchunk = 3 ; 
% cchunk = 3 ; 

for rchunk = 1:size(row_chunks,1)-1
for cchunk = 1:size(col_chunks,1)-1

 rchunk_start = row_chunks(rchunk) ;    
 rchunk_end   = row_chunks(rchunk+1) ;     
 cchunk_start = col_chunks(cchunk) ;    
 cchunk_end   = col_chunks(cchunk+1) ;      
 
 rows_cur = rchunk_start:rchunk_end ;
 cols_cur = cchunk_start:cchunk_end ; 
 rows_cur = repmat(rows_cur',[1 length(cols_cur)]) ; 
 cols_cur = repmat(cols_cur,[size(rows_cur,1) 1]) ; 
 
 
    SM_CCI_Passive_mean_subset = SM_CCI_Passive_mean(rchunk_start:rchunk_end,cchunk_start:cchunk_end) ; 
   
    SM_CCI_01_subset = squeeze(SM_CCI_01.SM_CCI(rchunk_start:rchunk_end,cchunk_start:cchunk_end,:)) ; 
    SM_CCI_02_subset = squeeze(SM_CCI_02.SM_CCI(rchunk_start:rchunk_end,cchunk_start:cchunk_end,:)) ; 
    SM_CCI_03_subset = squeeze(SM_CCI_03.SM_CCI(rchunk_start:rchunk_end,cchunk_start:cchunk_end,:)) ; 
    SM_CCI_04_subset = squeeze(SM_CCI_04.SM_CCI(rchunk_start:rchunk_end,cchunk_start:cchunk_end,:)) ; 
    SM_CCI_05_subset = squeeze(SM_CCI_05.SM_CCI(rchunk_start:rchunk_end,cchunk_start:cchunk_end,:)) ;     
    SM_CCI_06_subset = squeeze(SM_CCI_06.SM_CCI(rchunk_start:rchunk_end,cchunk_start:cchunk_end,:)) ;     
    SM_CCI_07_subset = squeeze(SM_CCI_07.SM_CCI(rchunk_start:rchunk_end,cchunk_start:cchunk_end,:)) ;     
    SM_CCI_08_subset = squeeze(SM_CCI_08.SM_CCI(rchunk_start:rchunk_end,cchunk_start:cchunk_end,:)) ;    
    SM_CCI_09_subset = squeeze(SM_CCI_09.SM_CCI(rchunk_start:rchunk_end,cchunk_start:cchunk_end,:)) ;    
    
    SM_CCI_full_dummy = cat(3,SM_CCI_01_subset,SM_CCI_02_subset,SM_CCI_03_subset,SM_CCI_04_subset,...
        SM_CCI_05_subset,SM_CCI_06_subset,SM_CCI_07_subset,SM_CCI_08_subset,SM_CCI_09_subset) ; 
     
    clear SM_CCI_01_subset SM_CCI_02_subset SM_CCI_03_subset SM_CCI_04_subset SM_CCI_05_subset
    clear SM_CCI_06_subset SM_CCI_07_subset SM_CCI_08_subset SM_CCI_09_subset
  
    
   % now get into 2D and only keep pixels that are in rowcollist (valid
   % land)
   SM_CCI_full_par = NaN(size(SM_CCI_full_dummy,1)*size(SM_CCI_full_dummy,2),size(SM_CCI_full_dummy,3) ) ;
   for i = 1 : size(SM_CCI_full_dummy,3)     
           SM_CCI_full_par(:,i) = reshape(  SM_CCI_full_dummy(:,:,i),[size(SM_CCI_full_dummy,1)*...
                                                              size(SM_CCI_full_dummy,2),1] )            ;   
   end
   
   
   rows_cur_par = reshape(  rows_cur,[size(rows_cur,1)* size(rows_cur,2)],1 )  ;     
   cols_cur_par = reshape(  cols_cur,[size(cols_cur,1)* size(cols_cur,2)],1 )  ;        
   %SM_CCI_Passive_mean_par =  reshape(  SM_CCI_Passive_mean_subset,[size(SM_CCI_Passive_mean_subset,1)* size(SM_CCI_Passive_mean_subset,2)],1 )  ;     
   mask_invalid =  all(isnan(SM_CCI_full_par),2) ; 
   SM_CCI_full_par(mask_invalid,:) = [] ; 
   rows_cur_par(mask_invalid,:) = [] ;    
   cols_cur_par(mask_invalid,:) = [] ;  
 
   
    %% see if array for chunk size is doable ... use less than 60 sm steps
   dSM_dt_interpsm_chunk_array = NaN(size(SM_CCI_full_par,1),500,60) ; 
   t_chunk = cell(size(SM_CCI_full_par,1),500,1) ;
   row_chunk =  NaN(size(SM_CCI_full_par,1),1) ; 
   col_chunk =  NaN(size(SM_CCI_full_par,1),1) ;    
   
   %% start with parfor loop to detect drydowns
   % if save time row col we can just merge chunks later   
   parfor DD = 1:size(SM_CCI_full_par,1)
   
   sm_dummy = SM_CCI_full_par(DD,:)' ; 
   
   cd('E:\Daten Baur\Matlab code\Project IGARSS multi frq tau\noodles_L_C_X')
   [NDryL,timevOL,smvOL] = DryDowns_SM(tt,sm_dummy,dLength); 
   % convert to arrays
    % for time
    cols = cellfun(@numel,timevOL);
    rows = NDryL ;
    timevOL_array = NaN(rows,max(cols));
    for k = 1:rows
    timevOL_array(k,1:cols(k)) = timevOL{k};
    end   
    
    % for sm
    cols = cellfun(@numel,smvOL);
    rows = NDryL ;
    smvOL_array = NaN(rows,max(cols));
    for k = 1:rows
    smvOL_array(k,1:cols(k)) = smvOL{k};
    end
    
    % calc dSM/dt 
    dSM_dt = diff(smvOL_array,1,2) ; 
    dt     = diff(timevOL_array,1,2) ; 
    dSM_dt = dSM_dt ./ dt ; 
    %dSM_dt(dSM_dt > 0) = NaN ; 
    % get short sm 
    smvOL_array_short = NaN(size(smvOL_array)) ; 
    for k = 1:rows
    dummy =  smvOL_array(k,:) ; 
    dummy(isnan(dummy)) = [] ; 
    smvOL_array_short(k,1:length(dummy)-1) = dummy(1:end-1) ; 
    end
    
    % interpolate on SM conditions
    dSm_dt_intersm_array = NaN(rows,60) ; 
    
    for k = 1:rows
    dummy_sm = smvOL_array_short(k,:)' ; 
    dummy_sm(isnan(dummy_sm)) = [] ; 
    dummy_dSM_dt = dSM_dt(k,:)' ; 
    dummy_dSM_dt(isnan(dummy_dSM_dt)) = [] ;  
    % check for fastest 1d interp methods .. multiple options maybe use
    % profiler to check difference mayb euse !interp1q!
    if ( issorted(dummy_sm,'descend') && (length(dummy_sm) == length(unique(dummy_sm)))  )
        
    dSm_dt_intersm_array(k,:) = interp1((dummy_sm),(dummy_dSM_dt),sminterp)  ;
   % dSm_dt_intersm_array(k,:) = interp1q((dummy_sm),(dummy_dSM_dt),sminterp)  ;
    end
    end
    % plot(sminterp, mean(dSm_dt_intersm_array,1,'omitnan'))
   
   for j = 1:rows
   dSM_dt_interpsm_chunk_array(DD,j,:) = dSm_dt_intersm_array(j,:) ; 
   t_chunk(DD,j) = {timevOL_array(j,:)} ; 
   row_chunk = rows_cur_par(DD) ; 
   col_chunk = cols_cur_par(DD) ; 
   end
   
   % show current parfor iteration
   DD
   end
   
  % save chunk output maybe possible to aggregate later
  
   
   
   
end
end

%








%% passive
clear
filelist = string(ls('SM_CCI*')) ; 
filelist = filelist(1:9) ; 
ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
load('F:\ESA_CCI\means\SM_CCI_Passive_mean.mat')


SM_CCI_Passive_mean_vec = SM_CCI_Passive_mean(:) ; 

[CCI_rows, CCI_cols] = find(~isnan(SM_CCI_Passive_mean)) ; 
CCI_rowcol_list = [CCI_rows CCI_cols] ; 
total_length = size(CCI_rowcol_list,1) ; 

SM_CCI_Passive_2D = NaN(total_length,15402) ; 

index = 1 ; 

for i = 1:length(filelist)

  dummy = load(filelist(i)) ; 
  dummy = dummy.SM_CCI ;
  dummy = reshape(dummy, [size(dummy,2)*size(dummy,1) , size(dummy,3)]) ; 
  dummy(isnan(SM_CCI_Passive_mean_vec),:) = [] ; 

  SM_CCI_Passive_2D(:,index:index + size(dummy,2) - 1) = dummy ; 

  index = index + size(dummy,2)  ; 
  i
  
end


save('F:\ESA_CCI\processed\SM_CCI_Passive_2D','SM_CCI_Passive_2D','-v7.3')




%% combined
clear
cd('F:\ESA_CCI\processed_combined')

ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
load('F:\ESA_CCI\means\SM_CCI_Combined_mean.mat')
% SM_CCI_Combined_mean(CCI_lat > 70) = NaN ; 

[CCI_rows, CCI_cols] = find(~isnan(SM_CCI_Combined_mean)) ; 
CCI_rowcol_list = [CCI_rows CCI_cols] ; 
total_length = size(CCI_rowcol_list,1) ; 
sequence = (1:5000:total_length)' ; 
sminterp = (0.01:0.01:0.6)' ; 


filelist = string(ls('SM_CCI*')) ; 
filelist = filelist(2:end) ; 

SM_CCI_Passive_mean_vec = SM_CCI_Combined_mean(:) ; 

sum(~isnan(SM_CCI_Passive_mean_vec)) ; 

SM_CCI_Combined_2D = NaN(total_length,15402) ; 

index = 1 ; 

for i = 1:length(filelist)

  dummy = load(filelist(i)) ; 
  dummy = dummy.SM_CCI_combined ;
  dummy = reshape(dummy, [size(dummy,2)*size(dummy,1) , size(dummy,3)]) ; 
  dummy(isnan(SM_CCI_Passive_mean_vec),:) = [] ; 

  SM_CCI_Combined_2D(:,index:index + size(dummy,2) - 1) = dummy ; 

  index = index + size(dummy,2)  ; 
  i
  
end



save('F:\ESA_CCI\processed_combined\SM_CCI_Combined_2D','SM_CCI_Combined_2D','-v7.3')









