clc ;
close all ;
clear ;
rng(0) ;
%% Define Parameters
subcarriers = 64 ;
cp = 16 ; % cyclic prefix
data_symbols = 200 ;
pilot_symbols = 4 ;
OFDM_symbols = data_symbols + pilot_symbols ;
Tx = 4 ; 
Rx = 4 ;
M = 4 ; % QPSK
bits_per_symbol = log2(M) ;
SNR_db = 0:4:24 ;
num_frames = 2 ;
taps = 6 ;
tap_delay = 8 ;
symbols_per_frame = OFDM_symbols ;
QAM_symbols = subcarriers*data_symbols*Tx ;
bits_per_frame = QAM_symbols*bits_per_symbol ;
%% PILOT DESIGN:
if pilot_symbols < Tx 
    error('pilot symbols must be >=  no. of transmit antennas for orthogonal pilot-by-time-slot scheme..') ;
end
% pilot QAM randomized
pilot_seq = zeros(subcarriers,pilot_symbols,Tx) ;
for k = 1:Tx
    pilot_seq(:,:,k) = qammod(randi([0 M-1],subcarriers,pilot_symbols),M,'UnitAveragePower',true) ;
    pilot_seq(:,:,k) = pilot_seq(:,:,k)*exp(1*j*2*pi*(k-1)/Tx) ;
end
BER = zeros(length(SNR_db),1) ;
Capacity = zeros(length(SNR_db),1) ; % average sum capacity per subcarrier (bits/sec/Hz)
for idx = 1:length(SNR_db) % MONTE CARLO LOOP TECHNIQUE
    total_error = 0 ;
    total_bits = 0 ;
    cap_Acc = 0 ;
    linear_SNR = 10^(SNR_db(idx)/10) ;
    noise_var = 1/linear_SNR ;
    %% TRANSMITTER: MIMO ENCODER: 
    for frame = 1:num_frames
        bits = randi([0 1], bits_per_frame, 1);
        sym_vec = qammod(bi2de(reshape(bits,bits_per_symbol,[]).','left-msb'),M,'UnitAveragePower',true) ;
        tx_data = reshape(sym_vec,subcarriers,data_symbols,Tx) ;
        sym_tx = zeros(subcarriers,symbols_per_frame,Tx) ; %full frame symbol matrix
        for k = 1:Tx 
            pilot_slot_idx = k ;
            sym_tx(:, pilot_slot_idx,k) = pilot_seq(:,pilot_slot_idx,k) ;
        end
        sym_tx(:,pilot_symbols + 1:end,:) = tx_data ;
        % OFDM MODULATOR:
        tx = ifft(sym_tx,subcarriers,1) ;
        tx_cp = [tx(end-cp+1,:,:); tx] ;
        %% CHANNEL: TIME DOMAIN MULTIPATH 
        h_time = zeros(Rx,Tx,taps) ; % impulse resopnse taps
        for r = 1:Rx 
            for k = 1:Tx
                pdp = exp(-0:taps-1) ; % decay(exp)
                tap = (randn(1,taps) + 1*j*randn(1,taps))./sqrt(2) ;
                tap = tap.* sqrt(pdp/sum(pdp)) ; % total power normalized to 1
                h_time(r,k,:) = tap ;
            end
        end
        H_freq = zeros(Rx,Tx,subcarriers) ;
        for r = 1:Rx 
            for k = 1:Tx
                h_vec = squeeze(h_time(r,k,:)) ;
                H_freq(r,k,:) = fft(h_vec, subcarriers); % Frequency response
            end
        end
        % POWER NORMALIZATION:
        sym_tx = sym_tx/sqrt(Tx) ;
        tx_cp = tx_cp/sqrt(Tx) ;
        %% FREQ. SELECTIVE CHANNEL TRANSMISSION:
        rx_symb_freq = zeros(subcarriers,symbols_per_frame,Rx) ;
        for OFDM = 1:symbols_per_frame
            x_mat = squeeze(sym_tx(:,OFDM,:)).' ;
            for k = 1:subcarriers 
                H = squeeze(H_freq(:,:,k)) ;
                X = x_mat(:,k) ;
                noise = sqrt(noise_var/2)*(randn(Rx,1) + 1*j*randn(Rx,1)) ;
                Y = H * X + noise ;
                rx_symb_freq(k,OFDM,:) = Y ;
            end
        end
        %% RECIEVER:
        H_freq_est = zeros(Rx,Tx,subcarriers) ;
        for k = 1:Tx 
            pilot_slot = k ;
            Y_pilot = squeeze(rx_symb_freq(:,pilot_slot,:)) ;
            Y_pilot = Y_pilot.' ;
            pilot_symb = squeeze(sym_tx(:,pilot_slot,k)) ;
            for j = 1:subcarriers
                if abs(pilot_symb(j)) < 1e-12 
                    H_freq_est(:,k,j) = zeros(Rx,1) ;
                else
                    H_freq_est(:,k,j) = Y_pilot(:,j) ./ pilot_symb(j) ; % LS estimation
                end
            end
        end
        %% MIMO DECODER:
        % EQUALIZATION PER SUBCARRIER:
        symb_est = zeros(subcarriers,data_symbols,Tx) ;
        C_subcarrier = zeros(subcarriers,data_symbols) ; % per subcarrier capacity
        SNR_subcarrier = zeros(subcarriers,data_symbols) ; % per subcarrier SNR
        for OFDM_symb = 1:data_symbols
            sym_idx = pilot_symbols + OFDM_symb ;
            for j = 1:subcarriers
                Y = squeeze(rx_symb_freq(j,sym_idx,:)) ;
                H_est = squeeze(H_freq_est(:,:,j)) ;
                % MMSE Equalization:
                W = (H_est' * H_est + noise_var *eye(Tx)) \ (H_est') ;
                X_hat = W * Y ;
                symb_est(j,OFDM_symb,:) = X_hat ;
                H_true = squeeze(H_freq(:,:,j)) ;
                C_subcarrier(j,OFDM_symb) = real(log2(det(eye(Rx) + (linear_SNR/Tx) * (H_true*H_true')))) ;
                cap_Acc = cap_Acc + C_subcarrier(j,OFDM_symb) ;
                noise_vec = Y - H_true*X_hat ;
                SNR_subcarrier(j,OFDM_symb) = mean(abs(H_true*X_hat).^2)/mean(abs(noise_vec).^2) ;
            end
        end
        % Average capacity and SNR per subcarrier over all OFDM symbols
        C_subcarrier_avg = mean(C_subcarrier, 2) ;
        SNR_subcarrier_dB = 10*log10(mean(SNR_subcarrier, 2)) ;
        % DEMODULATION:
        symb_data = symb_est(:,pilot_symbols+1:end,:) ;
        rx_symb_vec = reshape(symb_est,[],1) ;
        rx_demod_idx = qamdemod(rx_symb_vec,M,'UnitAveragePower',true) ;
        rx_bits = de2bi(rx_demod_idx,bits_per_symbol,'left-msb')' ;
        rx_bits = rx_bits(:) ;

        tx_bits = bits(1:length(rx_bits)) ;
        errors = sum(tx_bits ~= rx_bits) ;
        total_error = total_error + errors ;
        total_bits = total_bits + length(rx_bits) ;
    end
    % BER CALCULATION:
    BER(idx) = total_error/total_bits ;
    Capacity(idx) = cap_Acc/(num_frames*subcarriers*data_symbols) ;
    fprintf('SNR %2d dB : BER = %.3e , Avg. Capacity = %.3f bits/sec/Hz \n',SNR_db(idx),BER(idx),Capacity(idx)) ;
end
% Constellation at a specific SNR:
snr_idx = 5 ;
OFDM_symb_plot = 1 ;
rx_symbols = squeeze(symb_est(:,OFDM_symb_plot,:)) ;
rx_vec = rx_symbols(:) ;
%% PLOTS:
figure(1) ;
semilogy(SNR_db,BER, '-o','LineWidth',1.6) ;
grid on ;
xlabel('SNR (dB)') ;
ylabel('Bit Error Rate') ;
title(sprintf('%d x %d MIMO - OFDM BER (MMSE : Min. Mean Squared Error)',Tx,Rx)) ;

figure(2) ;
plot(SNR_db,Capacity,'-s','LineWidth',1.6) ;
grid on ;
xlabel('SNR (dB)') ;
ylabel(' Average Sum Capacity (bits/sec/Hz)') ;
title(sprintf('Average Capacity per subcarrier (%d x %d)',Tx,Rx)) ;

figure(3) ;
plot(1:subcarriers,C_subcarrier_avg,'-o','LineWidth',1.5) ;
grid on ; 
xlabel('Subcarrier index') ;
ylabel('Capacity (bits/sec/Hz)') ;
title('Per Subcarrier Capacity (Average over OFDM Symbols)') ;

figure(4) ;
plot(1:subcarriers,SNR_subcarrier_dB,'-s','LineWidth',1.5) ;
grid on ; 
xlabel('Subcarrier index') ;
ylabel('SNR (dB)') ;
title('Per Subcarrier SNR (Average over OFDM Symbols)') ;

figure(5);
scatter(real(rx_vec), imag(rx_vec), '.') ;
xlabel('In-phase') ; 
ylabel('Quadrature') ;
title(sprintf('Constellation after MMSE Equalization at SNR = %d dB', SNR_db(snr_idx))) ;
grid on ;
