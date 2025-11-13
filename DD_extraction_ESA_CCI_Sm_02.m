%% MJB 07.03.2023 extract DDs and seasonally correct from ESA CCI Passive and Combined


clear

% get datetime array
ESA_CCI_datetime = datetime('01-Nov-1978'):days(1):datetime('31-Dec-2020') ; 
dLength = 4 ; 
tt = 1:15402 ; 
sminterp = (0.01:0.01:0.6)' ; 
timeinterp = 1:1:50 ; 


%% define drydown function to not ot cd

% Isolate Soil Moisture DryDowns in SM and VOD Time-Series
% Place Individual DryDowns in SMv and VODv Cell Vector With Length Ndry
% DryDowns of Minimum DDLength Days Duration
cd('E:\Daten Baur\Matlab code\Project IGARSS multi frq tau\noodles_L_C_X')
DD_find_function = @DryDowns_SM ; 


%% build data arrys .. maybe prepare for parfor 
load('F:\ESA_CCI\means\SM_CCI_Combined_mean.mat')
load('F:\ESA_CCI\CCI_lon.mat')
load('F:\ESA_CCI\CCI_lat.mat')

CCI_lat = repmat(CCI_lat,[1 1440]) ; 
CCI_lon = repmat(CCI_lon',[720 1]) ; 

% filter snow and ice areas optional
%SM_CCI_Passive_mean(CCI_lat > 70) = NaN ; 

%get row and col list for extracting data 
[CCI_rows, CCI_cols] = find(~isnan(SM_CCI_Combined_mean)) ; 
CCI_rowcol_list = [CCI_rows CCI_cols] ; 



%% load 2D data
cd('F:\ESA_CCI\processed_combined')
%load('SM_CCI_Combined_2D.mat')

% matfile and parallelize in chunks
SM_CCI_Combined_2D = matfile('SM_CCI_Combined_2D.mat') ; 

% maybe smaller chunks so we can do it computationally
row_chunks = round(linspace(1,213841,12)) ; 



%% start parfor loops over chunks



for rchunk = 8:size(row_chunks,2)-1
    
    curchunk_start = row_chunks(rchunk) ; 
    curchunk_end =   row_chunks(rchunk+1) ; 
    % load chunk
    SM_CCI_chunk = SM_CCI_Combined_2D.SM_CCI_Combined_2D(curchunk_start:curchunk_end,:) ; 
    CCI_rowcol_list_chunk = CCI_rowcol_list(curchunk_start:curchunk_end,:) ; 

    % output arrays
    dSM_dt_interpsm_chunk_array = NaN(size(SM_CCI_chunk,1),900,60) ; 
    dSM_dt_interpt_chunk_array = NaN(size(SM_CCI_chunk,1),900,50) ;     
    Sm_chunk_interpt_array = NaN(size(SM_CCI_chunk,1),900,50) ; 

    % maybe tchunk only saves min and max time stamp. should be enough
    t_chunk_start_array = NaN(size(SM_CCI_chunk,1),900) ;
    t_chunk_end_array = NaN(size(SM_CCI_chunk,1),900) ;
    row_chunk_array =   NaN(size(SM_CCI_chunk,1),900) ; 
    col_chunk_array =  NaN(size(SM_CCI_chunk,1),900) ; 
    
    
    % profile on
    % mpiprofile on
    % parfor over dry downs
   % tic
     for DD = 1:size(SM_CCI_chunk,1)
         
         DD       
         sm_dummy = SM_CCI_chunk(DD,:)' ; 
%          rowcol_dummy = CCI_rowcol_list(DD,:) ;         
%          row_dummy = rowcol_dummy(1) ;          
%          col_dummy = rowcol_dummy(2) ; 
          row_dummy = CCI_rowcol_list_chunk(DD,1) ;           
          col_dummy = CCI_rowcol_list_chunk(DD,2) ;  

         
         % add a tinning step and read McColl for that issue
         % do we need to include the deseason step? Is this necessary for
         % this type of study? Maybe not, deseason AF steps did very little
         % anyways
         [NDryL,timevOL,smvOL] = DD_find_function(tt,sm_dummy,dLength); 

         cols = cellfun(@numel,timevOL);
         rows = NDryL ;
         timevOL_array = NaN(rows,max(cols));
         timevOL_array_max = NaN(900,1);
         timevOL_array_min = NaN(900,1);
          
         row_DD_array = NaN(900,1) ; 
         col_DD_array = NaN(900,1) ; 

        cols = cellfun(@numel,smvOL);
        smvOL_array = NaN(rows,max(cols));

         for k = 1:rows
         timevOL_array(k,1:cols(k)) = timevOL{k};
         timevOL_array_max(k) = max(timevOL{k});
         timevOL_array_min(k) = min(timevOL{k});     
         row_DD_array(k) = row_dummy ;
         col_DD_array(k) = col_dummy ;
         smvOL_array(k,1:cols(k)) = smvOL{k};
         end

         % calc dSM/dt 
         dSM_dt = diff(smvOL_array,1,2) ; 
         dt     = diff(timevOL_array,1,2) ; 
         % this might be key .. only allow dt of 2 weeks. minimum
         % difference between observations maybe be even more restrictive
         mask_dt = any(dt > 12,2) ; 

         dt(mask_dt,:) = NaN ;    
         timevOL_array(mask_dt,:) = NaN ; 
         smvOL_array(mask_dt,:) = NaN ; 
         col_DD_array(mask_dt) = NaN ; 
         row_DD_array(mask_dt) = NaN ; 
         timevOL_array_max(mask_dt) = NaN ; 
         timevOL_array_min(mask_dt) = NaN ; 

         % calc change rate
         dSM_dt = dSM_dt ./ dt ; 

         smvOL_array_short = NaN(size(smvOL_array)) ; 
         timevOL_array_short = NaN(size(smvOL_array)) ;   

         % interpolate on SM conditions
         dSm_dt_intersm_array = NaN(900,60) ; 
         Sm_intert_array = NaN(900,50) ; 
         dSm_dt_intert_array = NaN(900,50) ; 

         for k = 1:rows
         dummy_sm =  smvOL_array(k,:) ; 
         dummy_t = timevOL_array(k,:) ; 

         dummy_t(isnan(dummy_sm)) = [] ;         
         dummy_sm(isnan(dummy_sm)) = [] ; 

         dummy_dSM_dt = dSM_dt(k,:)' ; 
         dummy_dSM_dt(isnan(dummy_dSM_dt)) = [] ;  


         if ( issorted(dummy_sm,'descend') && (length(dummy_sm) == length(unique(dummy_sm))) && ~isempty(dummy_dSM_dt)  )

         dummy_reltime = dummy_t(1:end-1) - dummy_t(1) + 1;

         dSm_dt_intersm_array(k,:) = interp1((dummy_sm(1:end-1)),(dummy_dSM_dt),sminterp)  ;
         dSm_dt_intert_array(k,:) =  interp1((dummy_reltime),(dummy_dSM_dt),timeinterp)  ;    
         Sm_intert_array(k,:) = interp1((dummy_reltime),(dummy_sm(1:end-1)),timeinterp)  ; 

         end
         end

      % plot(sminterp, mean(dSm_dt_intersm_array,1,'omitnan'))
      % imagesc(dSm_dt_intersm_array)
       %  dSm_dt_intersm_array(dSm_dt_intersm_array == 0) = NaN ; 
         
             
         dSM_dt_interpsm_chunk_array(DD,:,:) = dSm_dt_intersm_array ; 
         dSM_dt_interpt_chunk_array(DD,:,:) = dSm_dt_intert_array ;          
         Sm_chunk_interpt_array(DD,:,:) = Sm_intert_array ; 

         % save only min and max time for DD
         t_chunk_start_array(DD,:) =timevOL_array_min ; 
         t_chunk_end_array(DD,:) = timevOL_array_max ;  
         row_chunk_array(DD,:) = row_DD_array ; 
         col_chunk_array(DD,:) = col_DD_array ; 

         
                               
     end
 % toc

 
          % save chunk output
         save(strcat('F:\ESA_CCI\DD_chunk_02\dSM_dt_interpsm_chunk_array_',num2str(rchunk)),'dSM_dt_interpsm_chunk_array','-v7.3')
         save(strcat('F:\ESA_CCI\DD_chunk_02\dSM_dt_interpt_chunk_array',num2str(rchunk)),'dSM_dt_interpt_chunk_array','-v7.3')
         save(strcat('F:\ESA_CCI\DD_chunk_02\Sm_chunk_interpt_array',num2str(rchunk)),'Sm_chunk_interpt_array','-v7.3')


         save(strcat('F:\ESA_CCI\DD_chunk_02\t_chunk_start_array',num2str(rchunk)),'t_chunk_start_array','-v7.3') 
         save(strcat('F:\ESA_CCI\DD_chunk_02\t_chunk_end_array',num2str(rchunk)),'t_chunk_end_array','-v7.3') 
         save(strcat('F:\ESA_CCI\DD_chunk_02\row_chunk_array',num2str(rchunk)),'row_chunk_array','-v7.3')
         save(strcat('F:\ESA_CCI\DD_chunk_02\col_chunk_array',num2str(rchunk)),'col_chunk_array','-v7.3')   

end










%% start parfor loops over chunks



for rchunk = 1:size(row_chunks,2)-1
    
    curchunk_start = row_chunks(rchunk) ; 
    curchunk_end =   row_chunks(rchunk+1) ; 
    % load chunk
    SM_CCI_chunk = SM_CCI_Combined_2D.SM_CCI_Combined_2D(curchunk_start:curchunk_end,:) ; 
    CCI_rowcol_list_chunk = CCI_rowcol_list(curchunk_start:curchunk_end,:) ; 
    % output arrays
    dSM_dt_interpt_chunk_array = NaN(size(SM_CCI_chunk,1),900,60) ; 
    % maybe tchunk only saves min and max time stamp. should be enough
    t_chunk_start_array = NaN(size(SM_CCI_chunk,1),900) ;
    t_chunk_end_array = NaN(size(SM_CCI_chunk,1),900) ;
    row_chunk_array =   NaN(size(SM_CCI_chunk,1),900) ; 
    col_chunk_array =  NaN(size(SM_CCI_chunk,1),900) ; 
    
    
    % profile on
    % mpiprofile on
    % parfor over dry downs
   % tic
     parfor DD = 1:size(SM_CCI_chunk,1)
         
         DD       
         sm_dummy = SM_CCI_chunk(DD,:)' ; 
%          rowcol_dummy = CCI_rowcol_list(DD,:) ;         
%          row_dummy = rowcol_dummy(1) ;          
%          col_dummy = rowcol_dummy(2) ; 
          row_dummy = CCI_rowcol_list_chunk(DD,1) ;           
          col_dummy = CCI_rowcol_list_chunk(DD,2) ;  

         
       
         [NDryL,timevOL,smvOL] = DD_find_function(tt,sm_dummy,dLength); 

         cols = cellfun(@numel,timevOL);
         rows = NDryL ;
         timevOL_array = NaN(rows,max(cols));
         timevOL_array_max = NaN(900,1);
         timevOL_array_min = NaN(900,1);
          
         row_DD_array = NaN(900,1) ; 
         col_DD_array = NaN(900,1) ; 
         
         for k = 1:rows
         timevOL_array(k,1:cols(k)) = timevOL{k};
         timevOL_array_max(k) = max(timevOL{k});
         timevOL_array_min(k) = min(timevOL{k});     
         row_DD_array(k) = row_dummy ;
         col_DD_array(k) = col_dummy ;
         end   
    

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
         timevOL_array_short = NaN(size(timevOL_array)) ; 
         for k = 1:rows
         dummy =  timevOL_array(k,:) -  timevOL_array(k,1) +1; 
         dummy(isnan(dummy)) = [] ; 
         timevOL_array_short(k,1:length(dummy)-1) = dummy(1:end-1) ; 
         end
    
         % interpolate on SM conditions
         dSm_dt_interpt_array = NaN(900,60) ; 
    
         for k = 1:rows
         dummy_t = timevOL_array_short(k,:)' ; 
         dummy_t(isnan(dummy_t)) = [] ; 
         dummy_dSM_dt = dSM_dt(k,:)' ; 
         dummy_dSM_dt(isnan(dummy_dSM_dt)) = [] ;  
         
         if ( issorted(dummy_t,'ascend') && (length(dummy_t) == length(unique(dummy_t)))  )
        
         dSm_dt_interpt_array(k,:) = interp1((dummy_t),(dummy_dSM_dt),1:60)  ;
       % dSm_dt_intersm_array(k,:) = interp1q((dummy_sm),(dummy_dSM_dt),sminterp)  ;
         end
         end
      % plot(sminterp, mean(dSm_dt_intersm_array,1,'omitnan'))
      % imagesc(dSm_dt_intersm_array)
         dSm_dt_interpt_array(dSm_dt_interpt_array == 0) = NaN ; 
         
         
     
         dSM_dt_interpt_chunk_array(DD,:,:) = dSm_dt_interpt_array ; 
         % save only min and max time for DD
%          t_chunk_start_array(DD,:) =timevOL_array_min ; 
%          t_chunk_end_array(DD,:) = timevOL_array_max ;  
%          row_chunk_array(DD,:) = row_DD_array ; 
%          col_chunk_array(DD,:) = col_DD_array ; 
                               
     end
 % toc

 
          % save chunk output
         save(strcat('F:\ESA_CCI\DD_chunk\dSM_dt_interpsm_chunk_array_',num2str(rchunk)),'dSM_dt_interpsm_chunk_array','-v7.3')
%          save(strcat('F:\ESA_CCI\DD_chunk\t_chunk_start_array',num2str(rchunk)),'t_chunk_start_array','-v7.3') 
%          save(strcat('F:\ESA_CCI\DD_chunk\t_chunk_end_array',num2str(rchunk)),'t_chunk_end_array','-v7.3') 
%          save(strcat('F:\ESA_CCI\DD_chunk\row_chunk_array',num2str(rchunk)),'row_chunk_array','-v7.3')
%          save(strcat('F:\ESA_CCI\DD_chunk\col_chunk_array',num2str(rchunk)),'col_chunk_array','-v7.3')         
end














