clc; clear all; close all;

fc = 918e6;
fs = 2e6;
B = 125e3;
SF = 7
M = 2^SF
T_sym = M/B
PRE_LEN = 32; PYLD_LEN = 32;
N = 128;
x_s = @(n_, s_) exp(1i*2*pi*((n_.^2)/(2*N) + (s_/N - 1/2)*n_))

k = 0:N-1;
x_ch = x_s(k, 0);

samplesPerFrame = 3*128*16*(PRE_LEN + 2 + PYLD_LEN)

rxPluto = sdrrx('Pluto', ...
    'RadioID','usb:0', ...'ip:192.168.2.8', ...
    'GainSource','Manual', ...
    'Gain', 30, ...
    'CenterFrequency', fc, ...
    'OutputDataType', 'double', ...
    'EnableBurstMode', false, ...
    'SamplesPerFrame', samplesPerFrame, ...
    'BasebandSampleRate', fs, ...
    'UseCustomFilter', false);

IQdataRX_raw = rxPluto();
IQdataRX =  IQdataRX_raw(:);

figure;
spectrogram(IQdataRX, 2048, 2040, 2048, fs, 'centered', 'xaxis')

x_pre_rx = downsample(IQdataRX, fs/B)

y = zeros(N, floor((numel(x_pre_rx))/N));
y_dn = zeros(N, floor((numel(x_pre_rx))/N));

k = 1;
for i = 1:N:(numel(x_pre_rx) - N)
    win = x_pre_rx(i : i+N-1);
    y(:, k) = fftshift(fft(win.*conj(x_ch(:)), N));
    y_dn(:, k) = fftshift(fft(win.*(x_ch(:)), N));
    s(k) = find(abs(y(:, k)) == max(abs(y(:, k)))) - (N/2 + 0)
    s_dn(k) = find(abs(y_dn(:, k)) == max(abs(y_dn(:, k)))) - (N/2 + 0)
    k = k+1;
end

for j = 1:(numel(s)-PRE_LEN+1)
    s_win = s(j: (j+PRE_LEN-1 )).'
    s_ref = s_win(1)*ones(numel(s_win), 1);
    if ( sum(s_ref == s_win) == PRE_LEN || ...
         sum(s_ref == s_win) == (PRE_LEN - 1) || ...
         sum(s_ref == s_win) == (PRE_LEN + 1) )
        break
    end
end
i_pre = j

L_up = round(mean(s(i_pre + [0:PRE_LEN-1])))
L_dn = round(mean(s_dn(i_pre + PRE_LEN + [0:1])))


% ------------------------------------------------------------------------
% estimate coarse CFO and STO
cfo_int = round((L_up + L_dn)/2)
sto_int = mod((L_up - cfo_int), N)

% figure;
% surf(abs(y), 'EdgeColor','none')
% hold on
% surf(abs(y_dn), 'EdgeColor','none')
hd = figure(198);
imagesc(20*log10(abs(y(:, ind_pre + [0:(PRE_LEN+2+PYLD_LEN-1)]))), [-15 Inf]); hold on;
colormap bone
clrmp = colormap(hd);
colormap(hd, flipud(clrmp))
%
% -------------------------------------------------------------------------
% coarse CFO correction by integer number of DFT bins
% apply STO correction
t_n = (0:(numel(x_pre_rx)-1))/B;
x_pre_rx1 = x_pre_rx.*exp(-1i*2*pi*(cfo_int*B/N).*t_n.')
x_pre_rx2 = x_pre_rx.*exp(1i*2*pi*(cfo_int*B/N).*t_n.')

x_pre_rx1 = circshift(x_pre_rx1, sto_int)
x_pre_rx2 = circshift(x_pre_rx2, -sto_int)

%
% -------------------------------------------------------------------------
% recompute DFT after coarse CFO and coarse STO
y = zeros(N, floor((numel(x_pre_rx))/N));
y_dn = zeros(N, floor((numel(x_pre_rx))/N));
phi = zeros(N, 1);
k = 1;
for i = 1:N:(numel(x_pre_rx) - N)
    win1 = x_pre_rx1(i : i+N-1);
    win2 = x_pre_rx2(i : i+N-1);
    y(:, k) = fftshift(fft(win1.*conj(x_ch(:)), 128));
    y_dn(:, k) = fftshift(fft(win2.*(x_ch(:)), 128));
    s(k) = find(abs(y(:, k)) == max(abs(y(:, k)))) - (N/2 + 0);
    s_dn(k) = find(abs(y_dn(:, k)) == max(abs(y_dn(:, k)))) - (N/2 + 0);
    phi(k) = angle(y(64, k))
    k = k+1
end

% figure;
% surf(abs(y))
% hold on
% surf(abs(y_dn))
hd = figure(199);
imagesc(20*log10(abs(y(:, ind_pre + [0:(PRE_LEN+2+PYLD_LEN-1)]))), [-15 Inf]); hold on;
colormap bone
clrmp = colormap(hd);
colormap(hd, flipud(clrmp))
%
% -------------------------------------------------------------------------
% compute fractional CFO
for j = 1:(numel(s)-PRE_LEN+1)
    s_win = s(j: (j+PRE_LEN-1 )).'
    s_ref = s_win(1)*ones(numel(s_win), 1);
    if ( sum(s_ref == s_win) == PRE_LEN) % || ...
         %sum(s_ref == s_win) == (PRE_LEN - 1) || ...
         %sum(s_ref == s_win) == (PRE_LEN + 1) )
        break
    end
end
ind_pre = j
%cfo_frac = (mean(diff(unwrap(phi(ind_pre+PRE_LEN/2 + [0:1]))))/(2*pi))*B/N
cfo_frac = (mean(diff(unwrap(phi(ind_pre + [0:PRE_LEN-1]))))/(2*pi))*B/N
%
t_n = (0:(numel(x_pre_rx1)-1))/B;
x_pre_rx1_1 = x_pre_rx1.*exp(-1i*2*pi*(cfo_frac).*t_n.')
%
% -------------------------------------------------------------------------
% recompute DFT after fine CFO
y = zeros(N, floor((numel(x_pre_rx1_1))/N));
k = 1;
for i = 1:N:(numel(x_pre_rx1_1) - N)
    win1 = x_pre_rx1_1(i : i+N-1);
    y(:, k) = fftshift(fft(win1.*conj(x_ch(:)), 128));
    s(k) = find(abs(y(:, k)) == max(abs(y(:, k)))) - (N/2 + 0);
    phi2(k) = angle(y(64, k))
    k = k+1
end

figure; plot(unwrap(phi), 'DisplayName', '\phi w/o CFO_{frac}'); hold on; 
plot(phi2, 'DisplayName', '\phi CFO_{frac}'); legend show; grid on;
hd = figure(200);
imagesc(20*log10(abs(y(:, ind_pre + [0:(PRE_LEN+2+PYLD_LEN-1)]))), [-15 Inf]); hold on;
colormap bone
clrmp = colormap(hd);
colormap(hd, flipud(clrmp))
%
% -------------------------------------------------------------------------
% identify samples corresponding to preamble upchirp and downchirp
s_pre = s(ind_pre + [0:PRE_LEN-1])
%
s_pyld_rx = s(ind_pre + PRE_LEN + 2 + [0:PYLD_LEN-1])
bits_pyld_rx = mod(s_pyld_rx, N) + 1
bits_pyld_rx = (de2bi(bits_pyld_rx, 7)).'
bits_pyld_rx = bits_pyld_rx(:)

msg_rx = char(bin2dec(reshape(char(bits_pyld_rx(:)+'0'), 8, []).')).'
% msg_rx =
%     'Packet sent by linear chirps'
%