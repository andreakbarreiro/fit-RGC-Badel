% From J's original code
% Original is in ../Julijana_code_orig_2010

% File w/ data
dfile = 'ON-parasol-dclamp.mat';
load(dfile);

%Vlist = lowc_same(k,:);
Vlist_array = highc_same; 
gelist = exc_high_same(1:end-1);
gilist = inh_highc_same(1:end-1);


% Parameters
dt    = .0001;
Ve = 0;
Vi = -.08;
Aexc = 30;
Ainh = 40;
    


[np,nq]=size(Vlist_array);

% Assume longer one is # time steps:
totalT  = (max(np,nq)-1)*dt;
nList   = min(np,nq);

tlist = 0:dt:totalT;


tslist_cell  = cell(nList,1);
tplist_cell  = cell(nList,1);
peaklist_cell = cell(nList,1);

for k=1:nList
    
    Vlist = Vlist_array(k,:);
    
    spikeparam=[];
    [tslist,tplist,peaklist]=get_spike_times_fn(Vlist,spikeparam);
       
    tslist_cell{k} = tslist;
    tplist_cell{k} = tplist;
    peaklist_cell{k}=peaklist;
   
end
  
click_through_sp_times_flag = 0;
if (click_through_sp_times_flag)
    nSub=nList;
    % Plot all spike times
    figure;
    for j1=1:nSub%8
        subplot(nSub,1,j1);plot(tlist,Vlist_array(j1,:));hold on;
        plot(tslist_cell{j1},0,'r*');
        plot(tslist_cell{j1}+.002,0,'b.');
        plot(tplist_cell{j1},peaklist_cell{j1},'g.');
        grid on;
    end
    maxT = 0.1;
    for k1=1:totalT/maxT
        for j1=1:nSub %8
            subplot(nSub,1,j1);
            xlim([(k1-1)*maxT, k1*maxT]);
        end
        pause;
    end
end % Click_through
