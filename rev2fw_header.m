function rev2fw_header()
%% this function generate reverse complimentary sequence of CCS for alignment


%fname = '/storage/home/h/hxc62/work/';%define the home folder
%addpath(fname);
%fname = 'E:/Hengye/';
fname = 'C:\Users\bailab\Dropbox';

%list = [{'17D'};{'17R'};{'ID'}; {'IR'};{'ND'};{'NR'};{'TD'};{'TR'}];%set the names of samples
%list = [{'I2R2'};{'I3R2'};{'N2R2'}; {'N3R2'};{'T2D2'};{'T2R2'};{'W1R2'};{'W4R2'}];%set the names of samples
list = [{'17RX'};{'N1R3'};{'N3R3'};{'O2R3'};{'O3R3'};{'T2D2X'};{'T2R2X'};{'T2R3'};{'TDX'};{'TRX'};{'W4R3'};{'W4R3X'}];
%mkdir('Y:\Hengye\Pacbio_analysis\re_start');
mkdir([fname,'\Pacbio_analysis\202006_comb']);
mkdir([fname,'\Pacbio_analysis\202006_comb\seqs']);
for k = 1:12
        name = list{k,1};% get the name from the list
        File1 = [fname, '\Pacbio_analysis\20200610\barcode_parsed_ccs\',name,'_ccs.fa'];% find the file corresponding to the sample
        [header1,seqs1] = fastaread(File1);% read the headers and sequences from the fastq file
        File2 = [fname, '\Pacbio_analysis\20200612\barcode_parsed_ccs\',name,'_ccs.fa'];% find the file corresponding to the sample
        %% only need to combine reads for duplicated samples
        [header2,seqs2] = fastaread(File2);% read the headers and sequences from the fastq file
        header = [header1,header2];%combine headers from two sequel runs
        seqs = [seqs1,seqs2];%combine reads from two sequel runs
        %% generate a new file containing reads and reversed reads
        [~,len] = size(seqs);
        seqs_c = cell(len*2,2);% create a new cell from headers and sequences
        for i = 1:len%input headers, seqs, and the reverse complimentary seqs to seqs_c
            seq = seqs{1,i};% get one read from the seqs
            seqs_c{i,2} = header{1,i};% set the second column of seqs_c as headers
            seqs_c{len+i,2} = header{1,i};
            seqs_c{i,1}= seq;
            seqs_c{len+i,1}= seqrcomplement(seq);
        end
    save_data([fname, '\Pacbio_analysis\202006_comb\seqs\',name,'.mat'],seqs_c);
end
