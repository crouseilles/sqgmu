function fct_plot(model,fft_b_adv_part,day)
% This function creates some plot online and save it

% Grid
x = model.grid.x;
y = model.grid.y;
My = model.grid.MX(2);
% Other parameters
taille_police = 12;id_part=1;
X0=[0 0];
T_adv_part = real(ifft2( fft_b_adv_part(:,:,1,id_part) ));

width = 3.3;height = 3.2;
figure1=figure(1);
set(figure1,'Units','inches','Position',[X0(1) X0(2) width height],'PaperPositionMode','auto');
imagesc(x,y,T_adv_part');

caxis([-1 1]*1e-3);
set(gca,'Units','normalized','FontUnits','points','FontWeight','normal',...
    'FontSize',taille_police,'FontName','Times')
ylabel('y(m)','FontUnits','points','interpreter','latex',...
    'FontSize',taille_police,'FontName','Times')
xlabel('x(m)','interpreter','latex','FontUnits','points',...
    'FontWeight','normal','FontSize',taille_police,'FontName','Times')
title({'One realization',['\hspace{0.5cm} $t=' num2str(day) '$ day ']},...
    'FontUnits','points','FontWeight','normal','interpreter','latex',...
    'FontSize',12,'FontName','Times')
axis xy; axis equal
colorbar
drawnow

%% Spectrum
X0=[3.3 1];
close(figure(4))
figure4=figure(4);

widthtemp = 12;heighttemp = 6;
set(figure4,'Units','inches','Position',[X0(1) X0(2) widthtemp heighttemp], ...
    'PaperPositionMode','auto');
fct_spectrum( model,fft_b_adv_part(:,:,:,id_part),'b');
set(gca,'XGrid','on','XTickMode','manual');
width = 4;height = 3;
set(figure4,'Units','inches','Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');
set(gca,'YGrid','on')

set(gca,'Units','normalized','FontUnits','points',...
    'FontWeight','normal','FontSize',taille_police,'FontName','Times')
ylabel('$|\hat{b}(\kappa)|^2$','FontUnits','points',...
    'interpreter','latex','FontSize',taille_police,'FontName','Times')
xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$','FontUnits','points',...
    'FontWeight','normal','FontSize',taille_police,...
    'interpreter','latex','FontName','Times')
title({'Spectrum of' '\hspace{0.5cm} one realization'},...
    'FontUnits','points','FontWeight','normal',...
    'interpreter','latex','FontSize',12,'FontName','Times')
drawnow

end

