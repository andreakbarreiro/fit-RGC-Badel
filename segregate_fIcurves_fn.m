function [ Vwin_cell,f_of_Vwin_cell ] = segregate_fIcurves_fn( Vlist_array,tslist_cell,...
    win_array,gelist,gilist,estCap,segparam)

% segregate fI curves
%   
 if (isfield(segparam,'Aexc'))
     Aexc = segparam.Aexc;
 else
     error('segregate_fIcurves_fn: You must specify segparam.Aexc');
 end
 if (isfield(segparam,'Ainh'))
     Ainh = segparam.Ainh;
 else
     error('segregate_fIcurves_fn: You must specify segparam.Ainh');
 end
 if (isfield(segparam,'Ve'))
     Ve = segparam.Ve;
 else
     error('segregate_fIcurves_fn: You must specify segparam.Ve');
 end
 if (isfield(segparam,'Vi'))
     Vi = segparam.Vi;
 else
     error('segregate_fIcurves_fn: You must specify segparam.Vi');
 end
 if (isfield(segparam,'tau_abs'))
     tau_abs = segparam.tau_abs;
 else
     error('segregate_fIcurves_fn: You must specify segparam.tau_abs');
 end

 % Default parameters
 dt = 0.0001;
 
 %
 remove_f_outliers = 0;
 if (isfield(segparam,'dt')); dt = segparam.dt; end; 
 if (isfield(segparam,'remove_f_outliers')); remove_f_outliers = segparam.remove_f_outliers; end;
  
 if (remove_f_outliers)
     if (isfield(segparam,'min_f_outlier'))
         min_f_outlier = segparam.min_f_outlier;
     else
         error('segregate_fIcurves_fn: You set segparam.remove_f_outliers, but did not set a min value'); 
     end
     if (isfield(segparam,'max_f_outlier'))
         max_f_outlier = segparam.max_f_outlier;
     else
         error('segregate_fIcurves_fn: You set segparam.remove_f_outliers, but did not set a max value'); 
     end
 end
 
% Check viability of windows
if (any(diff(win_array)<=0))
    error('segregate_fIcurves_fn: Windows must be increasing');
end



totalT  = (length(Vlist_array)-1)*dt;
tlist = 0:dt:totalT;

% Beginning, end of windows
win_begin   = [tau_abs win_array];
win_end     = [win_array Inf];

nWin = length(win_begin);

% Set up cells to store f values
Vwin_cell       = cell(nWin, 1);
f_of_Vwin_cell  = cell(nWin, 1);

% Are there multiple spike trains?
% i.e. Did I actually pass a cell of spike times?
tlistcell_flag = 0;  
if (iscell(tslist_cell)); tlistcell_flag = 1; end;

% How many spike time lists to process?
if (tlistcell_flag)
    nList = length(tslist_cell);
else
    nList=1;
end

for k=1:nList
%
    if (tlistcell_flag)
        tslist = tslist_cell{k};
    else
        tslist = tslist_cell;
    end
    
    % CONVERT VLIST TO 
    Vlist  = Vlist_array(k,:);
    
    
    % Make sure we know dt, Ve, Vi
    dVdtlist = diff(Vlist)/dt; %derivative of V
    %truncate all lists:  drop last point, so length matches dVdtlist
    Vlist=Vlist(1:end-1);
    
    %compute or load Iapplist
    Iapplist = -Aexc*gelist.*(Vlist-Ve) - Ainh*gilist.*(Vlist-Vi);
    
    % WE have estimated capacitance
    f_of_Vlist = dVdtlist - Iapplist/estCap;  

    win_labels=zeros(1,length(tlist));
    
     %march through the list of spike times.  Label timesteps w/in certain
    %ranges of every spike (moving forward) as to what window they fall into.

    for j=1:length(tslist)-1
        ts=tslist(j);
        tsnext=tslist(j+1);
        
        % Excludes a couple of time steps before spike time; must account
        % for our uncertainty in IDing it.
        tpad = 2*dt; 
        
        % Make sure we have tlist: relies on totalT and dt
        for j1=1:nWin
            win_labels( ts+win_begin(j1) <= tlist & tlist < min(ts+win_end(j1),tsnext-tpad) )=j1;
        end
        %pause;
        
    end
    % Record ...
    for j1=1:nWin
       Vtemp    = Vlist(win_labels==j1);
       ftemp    = f_of_Vlist(win_labels==j1);

       %%%%%
       %  IF ANY EXTRA PROCESSING HERE.....
       %%%%%
       if (remove_f_outliers)
          okind = find((ftemp <= max_f_outlier) & (ftemp >= min_f_outlier));
          %[length(okind) length(Vtemp) min(ftemp) max(ftemp)]
          
          %if (length(okind) < length(Vtemp))
          %    fprintf('Point removed (%d): orig. %d, new %d\n',j1,length(Vtemp),length(okind))
          %end
          Vtemp = Vtemp(okind);
          ftemp = ftemp(okind);
       end
       
       
        Vwin_cell{j1} = [Vwin_cell{j1} Vtemp];
        f_of_Vwin_cell{j1} = [f_of_Vwin_cell{j1} ftemp];
    end
    
end


end

