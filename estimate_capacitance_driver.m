% CALL AFTER get_spike_times_driver

% IF you want to override default parameters
capoutparam = [];

[estCap,moreInfo] = estimate_capacitance_Outer_fn(tslist_cell, Vlist_array, ...
    gelist, gilist, capoutparam);
estCap 
moreInfo

plot_estCap_all_flag = 1;
if (plot_estCap_all_flag)
    Vrest_keep          = moreInfo.Vrest_keep;
    estCap_all          = moreInfo.estCap_all;
    Var_array_all       = moreInfo.Var_array_all;
    
    [moreInfo.Vrest_keep' moreInfo.numV_array']
    
    Cap_array = [1:100 110:10:300]/1000;
    
    figure;set(gcf,'Position',[100,200,1000,400]);

    subplot(1,2,1);surf(Vrest_keep,Cap_array,log(Var_array_all));shading flat;
    view(0,90);colorbar;set(gca,'FontSize',16);
    xlabel('Vrest');ylabel('Potential C');

    subplot(1,2,2);
    plot(Vrest_keep,estCap_all,'*-');
    xlabel('Vrest');ylabel('Estimated C');
end