 %% this code was used to combine technical replicates or normalized biological replicates. 
 
 %name1 = 'O2R3';
 %name2 = 'O3R3';
% name1 = 'N1R3';
% name2 = 'N3R3';
% name1 = 'T2R2';
% name2 = 'T2R3';
name1 = 'TR';
name2 = 'T2RCn';
%fname = 'C:\Users\bailab\Dropbox\Pacbio_analysis\202006_comb\';
fname = 'E:\d\Dropbox\Pacbio_analysis\all_final\norm_';
% data1 = load([fname, 'matrix/normalized/matrix_norm_',name1,'.mat']);
% data2 = load([fname, 'matrix/normalized/matrix_norm_',name1,'.mat']);
data1 = load([fname, 'matrix/matrix_',name1,'.mat']);
data2 = load([fname, 'matrix/matrix_',name2,'.mat']);

for i = 1:169
    matrix2 = [data1.data(i).C_T_sum_trim;data2.data(i).C_T_sum_trim];
    pos = data1.data(i).C_T_pos;
        C_T_all_sum(i).C_T_sum_trim = matrix2;
        %C_T_all_sum(i).header = headers_trim;
        C_T_all_sum(i).C_T_pos = pos;%write C_T_sum, headers, and position to a structure
end

save_data([fname, 'matrix/matrix_TRall.mat'],C_T_all_sum);
