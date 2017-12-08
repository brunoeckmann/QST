function plotRho(rho,filename,titleString)
% PLOTRHO creates a 3D Histogram Plot of a density matrix rho
% Plot is saved to filename

    fig=figure()
    bar3(rho);
    title([titleString],'interpreter','latex')
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontSize',18)
    set(gca, 'YTickLabel', {'$\left| 00\right\rangle$' '' '' '$\left| 11\right\rangle$'})
    set(gca, 'XTickLabel', {'$\left\langle 00\right|$' '' '' '$\left\langle 11\right|$'})
    
    zlim([-0.25,1])
    zticks([-0.25 0 0.25 0.5 1])
    % esgibt noch \middle|
    view(230,30)
    saveas(fig,filename);
    close(fig)
end