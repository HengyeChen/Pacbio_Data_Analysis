function updown_dyad_oc_change()
%% this function calculated the change of position and occupancy of nucleosome -3 and -5.


dbstop if error
%fname = 'E:/d/Dropbox/';%set the home folder
fname = 'C:\Users\bailab\Dropbox\';
list = [{'IR2C'};{'NR3Cn'};{'OR3Cn'};{'TRall'};{'WTall'};{'WTcomb'};{'WR2C'}];

%% load data and pre-define matrix for data saving
data = load([fname 'Pacbio_analysis\all_final4\ud_dyad_pos.mat']);
wtupds = data.data.('WTall');
diffu = cell(169,1);
diffd = cell(169,1);
diffu1 = zeros(169,8);
diffd1 = zeros(169,8);
diffu2 = zeros(169,8);
diffd2 = zeros(169,8);
%% 
for i = 1:4
    name = list{i,1};
    upds = data.data.(name);
    %% calculate the difference between nucleosome position/occupancy of WT and mut samples
    for j = 1:169
        upd = upds{j,1};
        wtupd = wtupds{j,1};
        upd(upd(:,1)==0,:)=[];
        upd(upd(:,3)==-1208,:)=[];
        upd(upd(:,2)==-1208,:)=[];
        wtupd(wtupd(:,1)==0,:)=[];
        wtupd(wtupd(:,3)==-1208,:)=[];
        wtupd(wtupd(:,2)==-1208,:)=[];
        k=2;
        [diffu{j,1},diffu1(j,:),diffu2(j,:)] = diffud(upd(:,2),wtupd(:,2),k);
        [diffd{j,1},diffd1(j,:),diffd2(j,:)] = diffud(upd(:,3),wtupd(:,3),k);  
        
    end          
    all.(name).diffu = diffu;
    all.(name).diffu1 = diffu1;%upstream first nucleosome
    all.(name).diffu2 = diffu2;%upstream second nucleosome
    all.(name).diffd = diffd;
    all.(name).diffd1 = diffd1;%downstream first nucleosome
    all.(name).diffd2 = diffd2;%downstream second nucleosome
end
    save([fname 'Pacbio_analysis\all_final4\ud_dyad_shift.mat'],'all');
end

function [diff,diff1,diff2]= diffud(upd,wtupd,k) 
%% nucleosome was clustered into two clusters based on their position
%% mut sample
[idx,cl] = kmeans(upd,k);
p = zeros(k,1);
for x = 1:k
    p(x,1) = sum(idx==x)/length(idx);
end
st = sortrows([(1:k)',p,cl],3);
%% wt sample
[wtidx,wtcl] = kmeans(wtupd,k);
wtp = zeros(k,1);
for y = 1:k
    wtp(y,1) = sum(wtidx==y)/length(wtidx);
end
wtst = sortrows([(1:k)',wtp,wtcl],3);
%% difference
diff = [wtst,st,st(:,2:3)-wtst(:,2:3)];
diff1 = diff(1,:);
diff2 = diff(2,:);
end
