%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation for 《Estimation of Sparse Massive MIMO-OFDM Channels
% With Approximately Common Support》
% Input:
% Nt: Number of transmitting antenna
% P: Number of pilots
% Bandwidth: Bandwidth of OFDM symbols
% G: Number of measurement
% SNR_dB: --
% L: Number of pathes
% path_delay/path_power: Note that there are K path_delays and
%                       for each element there are L/K paths with different
%                       angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
Nt = 128;
P = 64;
N_subcarrier = 2048;
Bandwidth = 1e9;
delta_B = floor(N_subcarrier/P);
Pilot_subcarrier = (1:P)*delta_B;
G = 20;
SNR_dB = 0:10:30;
sample = 1e2;
noise_power = 1e-4;
L = 60;
fc = 2e9;
fN = fc+Bandwidth/2;
DFT_TX = dftmtx(Nt)/sqrt(Nt);
path_delay = [0,0.3,8.9,12.9,17.1,20]*1e-6;
path_power = [-2.5,0,-12.8,-10,-25.2,-16];
Path_each_delay = L/length(path_delay);
MSE = zeros(1,length(SNR_dB));
%% Simulation
for index_snr = 1:length(SNR_dB)
    SNR = 10.^(SNR_dB(index_snr)/10);
    for index = 1:sample
        % Generate AWGN,symbol and channel
        w = (randn(G,N_subcarrier)+1j*randn(G,N_subcarrier))*sqrt(noise_power/2);
        X = (randn(G,Nt,N_subcarrier)+1j*randn(G,Nt,N_subcarrier))...
            *sqrt(SNR*noise_power/2);
        H = zeros(Nt,N_subcarrier);
        for index_p = 1:length(path_delay)
            phi_tx = unifrnd(50,70,1,Path_each_delay)/180*pi;
            for index_subcarrier = 1:N_subcarrier
                fn = fN+(index_subcarrier-N_subcarrier)*Bandwidth/N_subcarrier;
                gain = 10.^(path_power(index_p)/10)*exp(1j*2*pi*rand(Path_each_delay,1))*exp(-1j*2*pi*path_delay(index_p)*fn);
                Num_Nt = (0:Nt-1)';
                A_TX = exp(-1j*2*pi*Num_Nt*sin(phi_tx))/sqrt(Nt);
                H(:,index_subcarrier) = H(:,index_subcarrier)+A_TX*gain;
            end
        end
        % Transmission
        y = zeros(G,N_subcarrier);
        for index_subcarrier = 1:N_subcarrier
            H_ang(:,index_subcarrier) = DFT_TX'*H(:,index_subcarrier);
            y(:,index_subcarrier) = X(:,:,index_subcarrier)*H(:,index_subcarrier)+w(:,index_subcarrier);
            PHI(:,:,index_subcarrier) = X(:,:,index_subcarrier)*DFT_TX;
        end
        [H_hat,NMSE,h_hat_new]= AMP_NNSPL(y,PHI,G,Nt,N_subcarrier,Pilot_subcarrier);
        %         plot(real(h_hat_new(1,:)));
        MSE(index_snr) = MSE(index_snr) + 1/Nt*sum(abs(H_hat-H).^2,'all');
    end
end
MSE = 10*log10(MSE/sample);
plot(SNR_dB,MSE);