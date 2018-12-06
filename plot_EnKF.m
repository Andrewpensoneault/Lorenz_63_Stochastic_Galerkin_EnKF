function[] = plot_EnKF(folder,plot_which,ma,C,true,ma2,C2,true2,Dt,n_cycle,nbv,Q,R,Q_init)
cwd = pwd;
cd(folder)
if (plot_which == 0 || plot_which == 2)
    %% Plots for EnKF
    figure()
    hold on
    plot3(ma(1,:),ma(2,:),ma(3,:))
    plot3(true(1,:), true(2,:), true(3,:), 'o')
    title({['Plot of the SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
    legend('Assimilated solution','Truth')
    %savefig(['Plot_3d_EnKF_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.fig'])
    saveas(gcf,['Plot_3d_EnKF_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.png'])
    pause(.1)
    close all
    
    figure()
    hold on
    errorbar(0:Dt:Dt*n_cycle,ma(1,:),squeeze(sqrt(C(1,1,:))))
    plot(0:Dt:Dt*n_cycle,true(1,:),'o')
    title({['Plot of the first dimension of the SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
    legend('Assimilated solution','Truth')
    xlabel(['Assimilation cycles (\Delta t =' num2str(Dt) ')'])
    ylabel('Units')
    %savefig(['Plot_1_EnKF_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.fig'])
    saveas(gcf,['Plot_1_EnKF_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.png'])
    pause(.1)
    close all
    
    figure()
    hold on
    errorbar(0:Dt:Dt*n_cycle,ma(2,:),squeeze(sqrt(C(2,2,:))))
    plot(0:Dt:Dt*n_cycle,true(2,:),'o')
    title({['Plot of the second dimension of the SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
    legend('Assimilated solution','Truth') 
    xlabel(['Assimilation cycles (\Delta t =' num2str(Dt) ')'])
    ylabel('Units')
    %savefig(['Plot_2_EnKF_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.fig'])
    saveas(gcf,['Plot_2_EnKF_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.png'])
    pause(.1)
    close all
    
    figure()
    hold on
    errorbar(0:Dt:Dt*n_cycle,ma(3,:),squeeze(sqrt(C(3,3,:))))
    plot(0:Dt:Dt*n_cycle,true(3,:),'o')
    title({['Plot of the third dimension of the SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
    legend('Assimilated solution','Truth')
    %savefig(['Plot_3_EnKF_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.fig'])
    saveas(gcf,['Plot_3_EnKF_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.png'])
    pause(.1)
    close all
end

if (plot_which == 1 || plot_which == 2)
    %% Comparison plots PCE-SREnKF
    figure()
    hold on
    plot3(ma2(1,:),ma2(2,:),ma2(3,:))
    plot3(true2(1,:), true2(2,:), true2(3,:), 'o')
    title({['Plot of the PCE-SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
    legend('Assimilated solution','Truth')
    %savefig(['Plot_3d_PCE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.fig'])
    saveas(gcf,['Plot_3d_PCE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.png'])
    pause(.1)
    close all
    
    
    
    figure()
    hold on
    errorbar(0:Dt:Dt*n_cycle,ma2(1,:),squeeze(sqrt(C2(1,1,:))))
    plot(0:Dt:Dt*n_cycle,true2(1,:),'o')
    title({['Plot of the first dimension of the PCE-SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
    legend('Assimilated solution','Truth')
    xlabel(['Assimilation cycles (\Delta t =' num2str(Dt) ')'])
    ylabel('Units')
    %savefig(['Plot_1_PCE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.fig'])
    saveas(gcf,['Plot_1_PCE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.png'])
    pause(.1)
    close all
    
    figure()
    hold on
    errorbar(0:Dt:Dt*n_cycle,ma2(2,:),squeeze(sqrt(C2(2,2,:))))
    plot(0:Dt:Dt*n_cycle,true2(2,:),'o')
    title({['Plot of the second dimension of the PCE-SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
    legend('Assimilated solution','Truth')
    xlabel(['Assimilation cycles (\Delta t =' num2str(Dt) ')'])
    ylabel('Units')
    %savefig(['Plot_2_PCE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.fig'])
    saveas(gcf,['Plot_2_PCE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.png'])
    pause(.1)
    close all
    
    figure()
    hold on
    errorbar(0:Dt:Dt*n_cycle,ma2(3,:),squeeze(sqrt(C2(3,3,:))))
    plot(0:Dt:Dt*n_cycle,true2(3,:),'o')
    title({['Plot of the third dimension of the PCE-SEnKF mean of the Lorenz 63 model with '] ['nbv=' num2str(nbv) ' and Q=' num2str(Q), ' and R=' num2str(R)]})
    legend('Assimilated solution','Truth')
    xlabel(['Assimilation cycles (\Delta t =' num2str(Dt) ')'])
    ylabel('Units')
    %savefig(['Plot_3_PCE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.fig'])
    saveas(gcf,['Plot_3_PCE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.png'])
    pause(.1)
    close all
end

if plot_which == 2 || plot_which == 3
%% Comparison Plots
    hold on
    plot(0:Dt:Dt*n_cycle,sqrt(sum((true-ma).^2)/3))
    plot(0:Dt:Dt*n_cycle,sqrt(sum((true2-ma2).^2)/3))
    title({['RMSE of truth and assimilated solution for at ensemble sizes ' num2str(nbv)] ['Q=' num2str(Q), ' and R=' num2str(R)]})
    legend(['SREnKF'],['PCE-EnKF'])  
    xlabel(['Assimilation cycles (\Delta t =' num2str(Dt) ')'])
    ylabel('RMSE')
    %savefig(['RMSE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.fig'])
    saveas(gcf,['RMSE_nbv_' num2str(nbv) '_Q_' num2str(Q) '_R_' num2str(R) '_Qinit_' num2str(Q_init) '.png'])
    pause(.1)
    close all
end
cd(cwd)
end
