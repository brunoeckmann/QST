function plotBehaviour(P_Ideal,P_Noise, simIndex,filename,titleString)
% PLOTBEHAVIOUR creates a 2D Plot of the P(ab|xy) vector (matrix)

    fig=figure('units','normalized','outerposition',[0 0 0.8 0.8]);
    
    plot(reshape(P_Ideal,36,1),'r-o')    
    hold on

    plot(reshape(P_Noise{simIndex,1},36,1),'c--d')
    plot(reshape(P_Noise{simIndex,2},36,1),'b--d')
    

    h=legend('$P(ab|xy)$ $\rho_{Ideal}$', '$P(ab|xy)$ $\rho_{Noise}$', '$P(ab|xy)$ $\rho_{Proj}$')
    
    set(h,'Interpreter','latex')
    title(titleString,'interpreter','latex')
    ylabel('P(ab|xy)')
    
%     set(gca,'TickLabelInterpreter', 'latex');
%     set(gca,'FontSize',18)
%     set(gca, 'YTickLabel', {'$\left| 00\right\rangle$' '' '' '$\left| 11\right\rangle$'})
%     set(gca, 'XTickLabel', {'$\left\langle 00\right|$' '' '' '$\left\langle 11\right|$'})
%     
%     zlim([-0.25,1])
%     zticks([-0.25 0 0.25 0.5 1])
%     % esgibt noch \middle|
%     view(230,30)
    
    fig.PaperOrientation='landscape';
    saveas(fig,filename);
    %close(fig)

end