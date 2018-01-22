function [E_cut_eff, nmbr_emttrs,Ephot_up] = plot_phs_mtch_prsr_p(lam,w0,b,Ions,Ephotmax,del_n,eta_crt)
        w0 = w0*1E-6;
        lam = lam*1E-6;
        E_cut_eff = Ephotmax(b);
        eta_up = Ions(1:b-1);
        Ephot_up = Ephotmax(1:b-1);
        Phs_mtch_prsr = 1013*lam^2./(2*pi^2*w0.^2*del_n*(1-eta_up/eta_crt));
        V_foc = pi^2 *w0.^4/(2*lam);
        nmbr_dnsty = 100*(Phs_mtch_prsr)/(1.38E-23*294); % 100 factor to convert to Pa
        nmbr_emttrs = V_foc .*nmbr_dnsty;
        
%         set(0,'DefaultAxesFontSize',22)
%         [AX,H1,H2] = plotyy(Ephot_up,nmbr_emttrs*1E-14,Ephot_up,Phs_mtch_prsr);
%         set(H1,'color','k')
%         set(H2,'color','k')
%         set(AX(2),'ycolor','k','linewidth',4,'xlim',[min(Ephot_up) 1.1*max(Ephot_up)],'ylim',[min(Phs_mtch_prsr)*0.8 max(Phs_mtch_prsr)*1.1],'xtick',linspace(0,1.1*max(Ephot_up),5),'xminortick','on','yminortick','on','box','off')
%         set(AX(1),'ycolor','k','linewidth',4,'xlim',[min(Ephot_up) 1.1*max(Ephot_up)],'ylim',[min(nmbr_emttrs)*0.8 1.1*max(nmbr_emttrs)]*1E-14,'xtick',round(linspace(0,1.1*max(Ephot_up),5)),'xminortick','on','yminortick','on','box','off')
%         set(AX(2), 'XTickLabel','','XAxisLocation','Top')
%         ylabel(AX(2),'Phasematching pressure (mbar)')
%         ylabel(AX(1),'# emitters (\times10^{14})')
%         xlabel('Cut-off photon energy (eV)')