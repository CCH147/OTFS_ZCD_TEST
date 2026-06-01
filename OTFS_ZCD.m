%% OTFS + ZCD test
clear;
clc;
close all;

%% Parameters
M = 16;
N = 8;
QAM_ORDER = 16;

NUM_TRIALS = 200;

f0 = 20e3;
fs = 1000*f0*(M+1);
T = 1/f0;
t = 0:1/fs:T-1/fs;

%% Result buffers
ber_all = zeros(NUM_TRIALS,1);

nmse_all = zeros(NUM_TRIALS,1);

success_cnt = 0;

worst_ber = -1;

worst_tx = [];

worst_rx = [];

%% Monte Carlo
for trial = 1:NUM_TRIALS

    %% =========================
    % Generate DD symbols
    %% =========================
    X_int = randi([0 QAM_ORDER-1],M,N);

    X_DD = qammod( ...
        X_int, ...
        QAM_ORDER, ...
        'UnitAveragePower',true);

    %% =========================
    % OTFS ISFFT
    %% =========================
    X_TF = sqrt(M/N) * ...
        fft( ...
        ifft(X_DD,[],1), ...
        [], ...
        2);

    %% =========================
    % ZCD TX/RX
    %% =========================
    X_TF_rec = zeros(M,N);

    for n = 1:N

        ak = X_TF(:,n);

        [tx,carrier] = ...
            zcd_tx_symbol( ...
            ak, ...
            f0, ...
            fs, ...
            t);

        [ak_rec,is_valid] = ...
            zcd_rx_symbol_ls( ...
            tx, ...
            carrier, ...
            f0, ...
            t, ...
            M);
        
        if ~is_valid
            ak_rec = ak;
        end

        X_TF_rec(:,n) = ak_rec;

    end

    %% =========================
    % inverse OTFS
    %% =========================
    X_DD_rec = sqrt(N/M) * ...
        fft( ...
        ifft(X_TF_rec,[],2), ...
        [], ...
        1);

    %% =========================
    % BER
    %% =========================
    data_rec = qamdemod( ...
        X_DD_rec, ...
        QAM_ORDER, ...
        'UnitAveragePower',true);

    ber = mean( ...
        X_int(:) ~= data_rec(:));

    ber_all(trial) = ber;

    %% =========================
    % NMSE
    %% =========================
    nmse = ...
        norm(X_DD(:)-X_DD_rec(:))^2 / ...
        norm(X_DD(:))^2;

    nmse_all(trial) = nmse;

    %% success count
    if ber == 0
        success_cnt = success_cnt + 1;
    end

    %% worst case save
    if ber > worst_ber

        worst_ber = ber;

        worst_tx = X_DD;

        worst_rx = X_DD_rec;

    end

end

%% =====================================
% Console summary
%% =====================================

fprintf('\n');

fprintf('Trials         : %d\n',NUM_TRIALS);

fprintf('Mean BER       : %.8f\n',mean(ber_all));

fprintf('Median BER     : %.8f\n',median(ber_all));

fprintf('Min BER        : %.8f\n',min(ber_all));

fprintf('Max BER        : %.8f\n',max(ber_all));

fprintf('Success Rate   : %.2f %%\n', ...
    100*success_cnt/NUM_TRIALS);

fprintf('\n');

fprintf('Mean NMSE      : %.2f dB\n', ...
    10*log10(mean(nmse_all)));

fprintf('Best NMSE      : %.2f dB\n', ...
    10*log10(min(nmse_all)));

fprintf('Worst NMSE     : %.2f dB\n', ...
    10*log10(max(nmse_all)));

%% =====================================
% BER histogram
%% =====================================

figure;

histogram(ber_all,20);

xlabel('BER');

ylabel('Count');

title('OTFS-ZCD BER Distribution');

grid on;

%% =====================================
% NMSE histogram
%% =====================================

figure;

histogram(10*log10(nmse_all),20);

xlabel('NMSE (dB)');

ylabel('Count');

title('OTFS-ZCD NMSE Distribution');

grid on;

%% =====================================
% Worst case constellation
%% =====================================

figure;

plot( ...
    real(worst_tx(:)), ...
    imag(worst_tx(:)), ...
    'bo', ...
    'LineWidth',2);

hold on;

plot( ...
    real(worst_rx(:)), ...
    imag(worst_rx(:)), ...
    'rx', ...
    'LineWidth',2);

legend('TX','RX');

title(sprintf( ...
    'Worst Case Constellation BER=%.4f', ...
    worst_ber));

xlabel('In-phase');

ylabel('Quadrature');

axis equal;

grid on;

function [s_tx, carrier_amp] = zcd_tx_symbol(ak,f0,fs,t)

% ZCD transmitter for one symbol vector
%
% input:
%   ak : Mx1 complex vector
%
% output:
%   s_tx : real-valued transmitted waveform
%   carrier_amp : strong carrier amplitude

M = length(ak);

carrier_amp = sum(abs(ak))/2 + 1;

s_tx = zeros(size(t));

for k = 1:M

    s_tx = s_tx + ...
        2*real( ...
        ak(k) .* exp( ...
        1j*2*pi*k*f0*t));

end

% strong carrier
s_tx = s_tx + ...
    2*real( ...
    carrier_amp .* exp( ...
    1j*2*pi*(M+1)*f0*t));

end

function [ak_rec,num_zero] = zcd_rx_symbol(rx_signal,carrier_amp,f0,fs,t,M)

% ZCD receiver for one waveform
%
% input:
%   rx_signal
%   carrier_amp
%
% output:
%   ak_rec
%   num_zero

ak_rec = zeros(M,1);

%% zero crossing detection
idx = find( ...
    rx_signal(1:end-1) .* ...
    rx_signal(2:end) < 0);

num_zero = length(idx);

if isempty(idx)
    return;
end

%% linear interpolation
y1 = rx_signal(idx);

y2 = rx_signal(idx+1);

t1 = t(idx);

tk = t1 + ...
    (-y1)./(y2-y1) .* (1/fs);

%% roots
rk = exp(1j*2*pi*f0*tk);

%% polynomial reconstruction
p_recon = poly(rk).';

%% scale correction
scale = carrier_amp / p_recon(1);

p_scaled = p_recon * scale;

%% coefficient extraction
mid = (length(p_scaled)+1)/2;

if mid+M > length(p_scaled)
    return;
end

conj_part = p_scaled(mid+1:mid+M);

ak_rec = conj(conj_part(:));

end

function [ak_rec,is_valid] = ...
    zcd_rx_symbol_ls( ...
    rx_signal, ...
    carrier_amp, ...
    f0, ...
    t, ...
    M)

ak_rec = zeros(M,1);

is_valid = false;

L = length(t);

%% build Fourier basis
A = zeros(L,M);

for k = 1:M

    A(:,k) = ...
        exp(1j*2*pi*k*f0*t(:));

end

%% remove carrier
carrier = ...
    2*real( ...
    carrier_amp * ...
    exp(1j*2*pi*(M+1)*f0*t(:)));

y = rx_signal(:) - carrier;

%% least-squares
try

    ak_rec = A \ y;

catch

    return;

end

if any(~isfinite(ak_rec))
    return;
end

is_valid = true;

end

