function [H,NMSE,h_hat_new] = AMP_NNSPL(y,PHI,G,Nt,N_subcarrier,Pilot_subcarrier)
kappa_old = rand(N_subcarrier,G);
z_old = zeros(N_subcarrier,G);
Gamma_old = rand(N_subcarrier,Nt);
Zeta_old = zeros(N_subcarrier,Nt);
Kai_old = rand(N_subcarrier,Nt);
Lambda_old = rand(N_subcarrier,Nt);
Mu_old = zeros(N_subcarrier,Nt);
Tao_old = rand(N_subcarrier,Nt);
h_hat_old = zeros(N_subcarrier,Nt);
v_old =rand(N_subcarrier,Nt);
rho_old= rand(N_subcarrier,Nt);
Eta_old = rand;
sigma_old = rand;
threshold1 = 40;
NMSE = zeros(1,threshold1);
for index_loop = 1:threshold1
%     tic
%% Update kappa
    for index_G = 1:G
        PHI1 = reshape(PHI(index_G,:,:),[Nt,N_subcarrier]);
        kappa_new(:,index_G) = sum(abs(PHI1.').^2.*v_old,2);
    end
%% Update z
    for index_G = 1:G
        PHI1 = reshape(PHI(index_G,:,:),[Nt,N_subcarrier]);
        z_new(:,index_G) = sum(PHI1.'.*h_hat_old,2)-kappa_new(:,index_G)./(sigma_old+kappa_old(:,index_G))...
            .*(y(index_G,:).'-z_old(:,index_G));
    end
%% Update Gamma
    Gamma_new = zeros(N_subcarrier,Nt);
    for index_Nt = 1:Nt
        PHI2 = reshape(PHI(:,index_Nt,:),[G,N_subcarrier]);
        Gamma_new1(:,index_Nt) = sum(abs(PHI2.').^2./(sigma_old+kappa_new),2);
    end
    Gamma_new = 1./Gamma_new1;
%% Update Zeta
    parameter1 = zeros(N_subcarrier,Nt);
    for index_Nt = 1:Nt
        PHI2 = reshape(PHI(:,index_Nt,:),[G,N_subcarrier]);
        parameter1(:,index_Nt) = sum(conj(PHI2).*(y-z_new.')...
            ./(sigma_old+kappa_new.'),1);
    end
    Zeta_new = h_hat_old+Gamma_new.*parameter1;
%% Update Kai, Lambda, Mu, Tao, h_hat, v
    Kai_new = log(Gamma_new./(Gamma_new+Eta_old))+abs(Zeta_new).^2./Gamma_new-abs(Zeta_new).^2./(Gamma_new+Eta_old);
    Lambda_new = rho_old./(rho_old+(1-rho_old).*(exp(-Kai_new)));
    Mu_new = Eta_old.*Zeta_new./(Gamma_new+Eta_old);
    Tao_new = Eta_old.*Gamma_new./(Eta_old+Gamma_new);
    h_hat_new = Lambda_new.*Mu_new;
    v_new = Lambda_new.*(abs(Mu_new).^2+Tao_new)-abs(h_hat_new).^2;
%% Update rho (divide 4 is necessary for all elements)
    rho_new = zeros(N_subcarrier,Nt);
    rho_new(1,1) = (Lambda_new(1,2)+Lambda_new(2,1))/4;
    rho_new(N_subcarrier,Nt) = (Lambda_new(N_subcarrier,Nt-1)+Lambda_new(N_subcarrier-1,Nt))/4;
    rho_new(1,Nt) = (Lambda_new(1,Nt-1)+Lambda_new(2,Nt))/4;
    rho_new(N_subcarrier,1) = (Lambda_new(N_subcarrier-1,1)+Lambda_new(N_subcarrier,2))/4;
    for index_N = 2:N_subcarrier-1
        rho_new(index_N,1) = (Lambda_new(index_N-1,1)+Lambda_new(index_N+1,1)+Lambda_new(index_N,2))/4;
        rho_new(index_N,Nt) = (Lambda_new(index_N-1,Nt)+Lambda_new(index_N+1,Nt)+Lambda_new(index_N,Nt-1))/4;
    end
    for index_Nt = 2:Nt-1
        rho_new(1,index_Nt) = (Lambda_new(1,index_Nt-1)...
            +Lambda_new(1,index_Nt+1)+Lambda_new(2,index_Nt))/4;
        rho_new(N_subcarrier,index_Nt) = (Lambda_new(N_subcarrier,index_Nt-1)...
            +Lambda_new(N_subcarrier,index_Nt+1)+Lambda_new(N_subcarrier-1,index_Nt))/4;
    end
    for index_N = 2:N_subcarrier-1
        for index_Nt = 2:Nt-1
            rho_new(index_N,index_Nt) = (Lambda_new(index_N-1,index_Nt)+Lambda_new(index_N+1,index_Nt)...
                +Lambda_new(index_N,index_Nt-1)+Lambda_new(index_N,index_Nt+1))/4;
        end
    end
%% Update Eta
    Eta_new = sum(Lambda_new(Pilot_subcarrier,:).*(abs(Mu_new(Pilot_subcarrier,:)).^2+Tao_new(Pilot_subcarrier,:)),'all');
    parameter1 = sum(Lambda_new(Pilot_subcarrier,:),'all');
    Eta_new = Eta_new/parameter1;
%% Update m
    m_new = (z_new*sigma_old+transpose(y).*kappa_new)./(sigma_old+kappa_new);
%% Update V
    V_new = sigma_old*kappa_new./(sigma_old+kappa_new);
%% Update sigma
    sigma_new = sum(abs(y(:,Pilot_subcarrier).'-m_new(Pilot_subcarrier,:)).^2+V_new(Pilot_subcarrier,:),'all');
    sigma_new = sigma_new/(length(Pilot_subcarrier)*G);
%% Check convergence
    if index_loop ~=1
        NMSE(index_loop) = sum(abs(h_hat_old-h_hat_new).^2,'all')/sum(abs(h_hat_old).^2,'all');
    else
        NMSE(1) = 1000;
    end
    if NMSE(index_loop)<1e-4
        break;
    end
%% Update for next loop
    kappa_old = kappa_new;
    z_old = z_new;
    Gamma_old = Gamma_new;
    Zeta_old = Zeta_new;
    Kai_old = Kai_new;
    Lambda_old = Lambda_new;
    Mu_old = Mu_new;
    Tao_old = Tao_new;
    h_hat_old = h_hat_new;
    v_old = v_new;
    rho_old = rho_new;
    Eta_old = Eta_new;
    sigma_old = sigma_new;
    parameter_test = (sigma_old+v_old)/G*Nt;
%     toc
%     disp(['运行时间: ',num2str(toc)]);
end
% figure(2);
% plot(abs(h_hat_new(1,:)));
DFT_TX = dftmtx(Nt)/sqrt(Nt);
H = DFT_TX*transpose(h_hat_new);
end