function met_normalzation_diff_samples()

%fname = '/storage/home/h/hxc62/work/';%define the home folder
%addpath(fname);
%fname = 'E:/d/Dropbox/';%set the home folder
fname = 'C:\Users\bailab\Dropbox\';

list1 = [{'17D'};{'17R'};{'ID'}; {'IR'};{'ND'};{'NR'};{'TD'};{'TR'}];% set the list of sample names
%list2 = [{'WR2C'};{'NR2C'};{'IR2C'};{'T2D2'};{'T2R2'}];
%list2 = [{'N1R3'};{'N3R3'};{'O2R3'};{'O3R3'};{'T2R3'};{'W4R3'}] ;
list2 = [{'NR3C'};{'OR3C'};{'T2RC'}];
motifs = load([fname, 'Pacbio_analysis/ref/motif_pos_v2.mat']);%read motif positions
mkdir([fname, 'Pacbio_analysis/202006_comb/matrix/normalized']);%make folder to save nucleosome prediction data
%%%%%%%%%%%%%%%%%%%%%
con_num = 11;% set the first con_num constantly methylated GC as control

for k = 3
    %% read unnormalized reads and control reads; create empty matrix for saving data
    name1 = list1{2,1};% select the control sample name from the list
    name2 = list2{k,1};% select the sample which needs to be normalized
    data = load([fname, 'Pacbio_analysis/re_start/matrix/matrix_comb_',name1,'_trimmed_more.mat']);%load GC matrix of WT sample, C in GC is 1, T in GC is -1, other nucleotides are 0.
    data2 = load([fname, 'Pacbio_analysis/202006_comb/matrix/unnormalized/matrix_',name2,'.mat']);
    tic;
    met_av_all = zeros(169,1);
    %met_av_long = zeros(169,11);
    met_av_all2 = zeros(169,1);
    met_av_long2 = zeros(169,11);
    count2 = 0;
    len_all = 0;
    len_all2 = 0;
    met_sum1 = zeros(169,11);
    met_sum2 = zeros(169,11);
    %% calculate the methylation level in the control region of all reads
    for x = 1:169
        position = xlsread([fname,'Pacbio_analysis\ref\GC_positions.xlsx'], x);% read the position of GCs
        %mpos = motifs.data(x).pos;%read the position file
        matrix = data.data(x).C_T_sum_trim;%read the trimmed sequence matrix
        matrix2 = data2.data(x).C_T_sum_trim;%read the trimmed sequence matrix
        [len,~] = size(matrix);%count the number of reads in the matrix
        [len2,~] = size(matrix2);
        GC_control = position(1,1:con_num);%use the first 11(con_num) GCs to normalize the rest of GCs
        if len~=0
        met_av_all(x,1) = sum(sum(matrix(:,GC_control)==1))/len;% methylation level in the control region of one matrix
        %met_av_long(x,:) = sum(matrix(:,GC_control)==1)/len;% methylation level at each GC position in the control region of one matrix
        len_all = len_all+len;
        met_sum1(x,:) = sum(matrix(:,GC_control)==1);
        end
        if len2~=0
            met_av_all2(x,1) = sum(sum(matrix2(:,GC_control)==1))/len2;
            met_av_long2(x,:) = sum(matrix2(:,GC_control)==1)/len2;
            met_sum2(x,:) = sum(matrix2(:,GC_control)==1);
            count2 = count2+1;% count2 is the number of non-empty matrix
            len_all2 = len_all2+len2;
        end       
    end
    %% normalized reads by the methylation levels in the constant region; 
    %GCs was stochastically added/removed to adjust the methylation level
    for i = 1:169
        position = xlsread([fname,'Pacbio_analysis\ref\GC_positions.xlsx'], i);
        matrix2 = data2.data(i).C_T_sum_trim;%read the trimmed sequence matrix
        [len2,~] = size(matrix2);%count the number of reads in the matrix2
             if len2==0
                continue
             end
        met_av_control = sum(met_av_all)/169/11;% the methylation level of the control
        met_av_real = sum(met_av_all2)/count2/11;% the methylation level of the target sample
            met_change = met_av_control-met_av_real;% absolute change of methylation
            met_change_ratio = met_av_control/met_av_real-1;% ratio change of methylation
     %% change the number of GCs in each position
            for n = 1:length(position) 
                if len2 == 0
                    continue
                end
                GC = position(1,n);
                met_level2 = sum(matrix2(:,GC)==1)/len2;
                %% calculate how much GCs need to add/remove
                if (met_level2+met_change)<0.95 
                    met_add = round(met_change_ratio*met_level2*len2);
                else % if the methylation level after adjustment is too high, keep the methylation level at 95%
                    met_add = round((0.95-met_level2)*len2);
                end
                
                if met_add >0% add GCs
                    if met_add<=sum(matrix2(:,GC)==-1)
                        add = randperm(sum(matrix2(:,GC)==-1),met_add)';% randomly choose reads and add GCs
                        unmetC = find(matrix2(:,GC)==-1);
                        unmetC_add = unmetC(add,1);
                        matrix2(unmetC_add,GC) = 1;
                    end
                elseif met_add <0% remove GCs
                    if -met_add<=sum(matrix2(:,GC)==1)
                        reduce = randperm(sum(matrix2(:,GC)==1),-met_add)';
                        metC = find(matrix2(:,GC)==1);
                        metC_reduce = metC(reduce,1);
                        matrix2(metC_reduce,GC) = -1;
                    end
                end
            end
        C_T_all_sum(i).C_T_sum_trim = matrix2;
        %C_T_all_sum(i).header = headers_trim;
        C_T_all_sum(i).C_T_pos = position(1,:);%write C_T_sum, headers, and position to a structure
    end
    save_data([fname, 'Pacbio_analysis/202006_comb/matrix/normalized/matrix_',name1,'_',name2,'_ratio.mat'],C_T_all_sum);
end