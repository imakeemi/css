clc; clear all; close all;

fc = 918e6;
fs = 2e6;
B = 125e3;
SF = 7
M = 2^SF
T_sym = M/(B)
PRE_LEN = 32
N = M*(fs/B);
x_s = @(n_, s_) exp(1i*2*pi*((n_.^2)/(2*N) + ((fs/B)*s_/N - 1/2)*n_))

k = ((0:1:N-1));
%
%plot(fs/(2*pi)*diff(unwrap(angle(x_s(k, 0)))/(fs/B)), 'DisplayName', 'S=0'); hold on; grid on;
%plot(fs/(2*pi)*diff(unwrap(angle(x_s(k, 32)))/(fs/B)), 'DisplayName', 'S=32'); legend show;
%plot(fs/(2*pi)*diff(unwrap(angle(x_s(k, 96)))/(fs/B)), 'DisplayName', 'S=96'); xlim([1 N])
%xlabel('n')
%ylim([-B/2 B/2]); yticks([-B/2 -B/4 0 B/4 B/2 ]); ylabel('Instantaneous Frequency (Hz)')
%
phi = unwrap(angle(x_s(k, 0)))/(fs/B);
x_ch_0 = exp(j*phi);
%
% -------------------------------------------------------------------------
% build standard PRE_LEN upchirps + 2 downchirps preamble

x_pre = repmat(x_ch_0(:), PRE_LEN, 1);
x_pre = [x_pre(:); conj(x_ch_0(:)); conj(x_ch_0(:))];
%
% -------------------------------------------------------------------------
% build 224 bits message
msg = 'Packet sent by linear chirps'
bits_d = reshape(dec2bin(msg, 8).'-'0',1,[]);

bits_d = reshape(bits_d, 7, []).';
s_d = bi2de(bits_d) - 1;

x_pyld = [];
for i = 1:numel(s_d)
    x_i = exp(j*unwrap(angle(x_s(k, s_d(i))))/(fs/B));
    x_pyld = [x_pyld; x_i(:)];
end
x_pk = [x_pre; x_pyld(:)];

IQdataTX = [N*zeros(1*128*16*(PRE_LEN + 2 + 28), 1); x_pk(:); N*zeros(1*128*16*(PRE_LEN + 2 + 28), 1)];
%IQdataTX = [x_pk(:);];

%
% -------------------------------------------------------------------------
% transmit signal
if max(abs(IQdataTX)) > 0
    IQdataTX = IQdataTX / max(abs(IQdataTX));
end

txPluto = sdrtx('Pluto', ...
    'RadioID','ip:169.254.27.106', ...
    'Gain', -15, ...
    'CenterFrequency', fc, ...
    'BasebandSampleRate', fs);

txPluto.ShowAdvancedProperties = true;

fprintf('Transmitting repeated LoRa packet...\n');
transmitRepeat(txPluto, IQdataTX);
