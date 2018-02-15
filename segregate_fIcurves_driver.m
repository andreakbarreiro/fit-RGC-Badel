%seek to fit:
%C dV/dt = I_ion + I_app
%rewrite
%dV / dt = f(V) + I_app /C
%

% ASSUME: 
%
% 1) We have spike times 
% i.e. you have already called
%
% get_spike_time_driver
%
% 2) We have an estimated capacitance
%
% i.e. you have already called
%
% estimate_capacitance_driver
%
% WE ALREADY HAVE: 
% 
% tslist_cell
% Vlist_array
% gelist
% gilist
% estCap
%

% Decide on windows
%win_array = [5 10 20 30 50 75 100 150 200]/1000;
win_array = [5 10 20 30 50 100 150]/1000;

% Set up parameters for fn
segparam.Aexc = Aexc;
segparam.Ainh = Ainh;
segparam.Ve = Ve;
segparam.Vi = Vi;
segparam.tau_abs = tau_ref;

[ Vwin_cell,f_of_Vwin_cell ] = segregate_fIcurves_fn( Vlist_array,tslist_cell,...
    win_array,gelist,gilist,estCap,segparam);



nWin = length(Vwin_cell);
% Beginning, end of windows (USE FOR PLOTS)
win_begin   = [segparam.tau_abs win_array];
win_end     = [win_array Inf];

% Fit FI curves...
pfit_cell = cell(nWin,1);
for j1=1:nWin
    pfit_cell_temp = fit_f_fminsearch(1000*Vwin_cell{j1},f_of_Vwin_cell{j1});
    pfit_cell{j1}=pfit_cell_temp;
    %pfit_cell_temp
end

% Another option: get averages. Then fit THIS
avg_f_cell = cell(nWin,1);
avg_V_cell = avg_f_cell;
for j1=1:nWin
    
    % Voltage in mV
    [Varray,farray]=get_fIcurve_AKB( Vwin_cell{j1}*1000, f_of_Vwin_cell{j1});

    avg_f_cell{j1} = farray;
    avg_V_cell{j1} = Varray;
    
end

plot_fIcurves_each_window=0;
if (plot_fIcurves_each_window)
    % Set up figures for f-I curve
    nC = ceil(nWin/2);
    nR = 2;

    figure;set(gcf,'Position',[100,200,1000,400]);
    for j1=1:nWin
        subplot(nR,nC,j1);
        plot(Vwin_cell{j1},f_of_Vwin_cell{j1},'.');
        grid on;set(gca,'FontSize',18);
        title(sprintf('%g-%g ms',1000*win_begin(j1),1000*win_end(j1)));

        % Round up/down to nearest 10 mV
        minV = floor(100*min(Vwin_cell{j1}))/100;
        maxV = ceil(100*max(Vwin_cell{j1}))/100;
        xlim([minV,maxV]);

        minF = max(min(f_of_Vwin_cell{j1}),-30);
        maxF = min(max(f_of_Vwin_cell{j1}),50);
        ylim([minF,maxF]);

        % Plot the fit f-I curve
        Varray = linspace(1000*minV,1000*maxV,100);

        % We fit in mV
        params = pfit_cell{j1};
        farray = params(1)*(params(2)-Varray+params(3)*exp((Varray-params(4))/params(3)));
        hold on;
        plot(Varray/1000,farray,'r-');
        plot(Varray/1000,zeros(size(Varray)),'k--');

        plot(avg_V_cell{j1}/1000,avg_f_cell{j1},'m.-');
    end
end
 