%function align_met_all()
%% this function align CCS to reference sequences

fname = '/storage/home/h/hxc62/work/';%define the home folder
addpath(fname);
%fname = 'E:/Hengye/';
%fname = 'C:\Users\BAI Lab\Dropbox\';
%list = [{'17D'};{'17R'};{'ID'}; {'IR'};{'ND'};{'NR'};{'TD'};{'TR'}];% set the name of samples
%list = [{'I2R2'};{'I3R2'};{'N2R2'}; {'N3R2'};{'T2D2'};{'T2R2'};{'W1R2'};{'W4R2'};{'test'}];
list = [{'N1R3'};{'N3R3'};{'O2R3'};{'O3R3'};{'T2R3'};{'W4R3'}] ;

mkdir([fname, 'Pacbio_analysis/202006_comb']);
mkdir([fname, 'Pacbio_analysis/202006_comb/aligned']);% make the main folder for all samples
%mkdir([fname, 'Pacbio_analysis/202006_comb/aligned']);% make a subfolder for the selected sample
for l = 6
    mkdir([fname, 'Pacbio_analysis/202006_comb/aligned/',num2str(l)])
    name = list{l,1};% select the file of one sample
    tic;
        fileID = fopen([fname, 'Pacbio_analysis/ref/refs_Y.txt']);%open the reference file(orignial Cs are converted to Ys) 
        refs = textscan(fileID,'%s');%read reference sequences
        refs = refs{1,1};%convert refs from a cell format to a char array
        fclose(fileID);
        
        data = load([fname, 'Pacbio_analysis/202006_comb/seqs/', name,'.mat']);%load the unmatched reads from file
        seqs_headers = data.data;
        seqs = seqs_headers(:,1);
        headers = seqs_headers(:,2);
        alignment = struct('reads', cell(1,169));%make a 1X169 structure named alignment to store the rematched reads 
        mot_pos = load([fname,'Pacbio_analysis/ref/motif_pos_v2.mat']);% read positions and sequences of motifs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        parfor r = 1:169
                ref=refs{r,1};%read a reference sequence
                aligned = [];
                mheader = cell(length(seqs),1);
                for i = 1:length(seqs)% align each seq to the selected ref
                    header = headers{i,1};
                    try
                        %% align each read to the reference sequence
                        [~, local_align] = swalign(ref,seqs{i,1},...%locally align unmatched seqs to the ref seq using Smith-Waterman algorithm
                           'Alphabet', 'NT',...%set seq type as nucleotide('NT')
                           'ScoringMatrix','NUC44');%use matrix NUC44 to align sequence, only the well-matched part will be aligned

                        %% In this module the indel mutations in the sequencing reads are removed, and the output data aseq has the same length as the ref seq
                        local_match = strrep(local_align(3,:),'-','N');%replace all the gaps in the alignment to N
                        [~,align] = nwalign(ref,local_match,...%globally align realigned seqs to the ref using Needleman-Wunsch algorithm
                                            'Alphabet', 'NT',...%set seq type as nucleotide('NT')
                                            'ScoringMatrix','NUC44');%use matrix NUC44 to align sequence, the whole sequence will be aligned
                                            
                        alignr = [align(1,:);align(3,:)];%get the aligned sequences from the alignment output data
                        [~,w] = size(alignr);% get the length of the aligned sequence    
                        for width = 1:w
                            if alignr(2,width)=='-'%replace the gaps in the sequence by N
                                alignr(2,width)=strrep(alignr(2,width),'-','N');
                            end
                            if alignr(1,width)=='-'%replace the insert mutations in the sequence by x
                                alignr(2,width)='x';
                            end
                        end

                        aseq=regexprep(alignr(2,:),'x','');% delete all x in the sequence, so that the length of aseq is same as the reference seq
                        %% determine whether the sequence is well aligned or not.                    
                        [~, align2] = nwalign(ref,aseq,...% globally align the trimmed aseq to ref to get the similarity
                                              'Alphabet', 'NT',...
                                              'ScoringMatrix','NUC44');
                                          
                         sim_all = sum(align2(2,:)== '|' | align2(2,:)==':');%count how many nt are matched in the whole seq
                         sim_v = sum(align2(2,430:540)== '|' | align2(2,430:540)==':');%count how many nt are matched in the variable region
                         %if the overall similarity is above 90% and similarity at variable region is higher than 99%, add this sequence to 'aligned'
                         if sim_all > 0.90*length(ref) && sim_v > 0.99*(540-430)
                             real = motif_align(aseq,mot_pos,r);
                             if real == 1
                                aligned = [aligned; aseq];
                                [len,~] = size(aligned);
                                mheader{len,1} = header;
                             end
                         end   
                    catch
                    end
                end
                alignment(1,r).reads = aligned; %input the aligned seq to the struct, r is the index of the ref
                alignment(1,r).header = mheader(~cellfun('isempty',mheader));%remove empty variables in the cell array mheader
                save_data([fname, 'Pacbio_analysis/202006_comb/aligned/',num2str(l),'/aligned_', name,'_',num2str(r), '.mat'], aligned)% save this set of remathced seq
        end
        save([fname, 'Pacbio_analysis/202006_comb/aligned/aligned_', name, '.mat'], 'alignment');%save all rematch seqs
    toc;        
end
