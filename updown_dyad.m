function updown_dyad()
%% this function calculate the positions of the nucleosomes in the upstream and downstream of motifs

%fname = '/storage/home/h/hxc62/work/';%define the home folder
%addpath(fname);
%fname = 'E:/d/Dropbox/';%set the home folder
fname = 'C:\Users\bailab\Dropbox\';
list = [{'IR2C'};{'NR3Cn'};{'OR3Cn'};{'TRall'};{'WTall'};{'WTcomb'};{'WR2C'}];% read sample list
motifs = load([fname '\Pacbio_analysis\ref\motif_pos_v2.mat']);%read motif positions
for k = 1:5
    name = list{k,1};
    data = load([fname 'Pacbio_analysis\all_final4\nuc_prediction\tt\nuc_pred_TTxx_',name,'_73_20.mat']);% read inferred nuc position
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all = cell(169,1);
    allp = [];
    for i = 1:169
        mot_pos = motifs.data(i).pos;
        peak_all = data.data(i).pred;
        peak_all(sum(peak_all,2)==0,:)=[];
        [len,~] = size(peak_all);
        %%%%%%%%%%%%%%%%%%%%%%%
        motu = mot_pos(1,1);%find the upstream edge and downstream edge of motifs
        motd = mot_pos(end,end);
        dist = zeros(len,3);
        dnuc = zeros(len,3);
        %% find the dyad of the up and downstream nucleosome in the presence or absence of NDR
        for j = 1:len
            peak = peak_all(j,:);
            if sum(peak(1,motu:motd)<=0)>0 % part of the motif is not protected
                dist(j,1) = 1;  % ndr exist
                dnuc(j,1) = 1;
                op = find(peak(1,motu:motd)==0);
                start = op(1,1)+motu-1;
                dnuc(j,2:3) = udnuc(peak,start)-1208;
            else % motif is protected
                opu = find(peak(1,1:motu)==0);
                opd = find(peak(1,motu:end)==0);
                startu = opu(end,end);
                startd = opd(1,1)+motu-1;
                uupnuc = udnuc(peak,startu);
                ddnuc = udnuc(peak,startd);
                dnuc(j,2:3) = [uupnuc(1) ddnuc(2)]-1208;
                
            end
        end
        all(i,1) = {dnuc};
        allp = [allp;dnuc];
    end
    dyad_all.(name) = all;
end
    save_data([fname 'Pacbio_analysis\all_final4\ud_dyad_pos.mat'],dyad_all);
end




function dyad = udnuc(peak,start)
%% find the dyad of the up and downstream nucleosome
dyad = [0,0];
try 
    upn = find(peak(1,1:start)>0);
    updown = upn(1,end);
    upm = find(peak(1,1:updown)==0);
    upup = upm(1,end);
    dyad(1,1) = round((upup+updown)/2);
catch
end
try
    downn = find(peak(1,start:end)>0);
    downup = downn(1,1)+start-1;
    downm = find(peak(1,downup:end)==0);
    downdown = downm(1,1)+downup-1;
    dyad(1,2) = round((downup+downdown)/2);
catch
end
    


end
