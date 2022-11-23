%%% 正交解调 %%%
% cosSig：输入的cos信号
% t：输入时间样点
% fc：信号的载波频率
% fs：信号的采样频率
% fs_out：输出信号的采样频率（需要大于基带带宽）

% 输出pluralSig：输出合成的复信号

function [pluralSig] = quadratureDemodulation(cosSig, t, fc, fs, fs_out)
  N = length(cosSig);
  half_win = round((fs_out/2)/(fs/N));
  L = 2*half_win;
  % Inphase channel
  I_sig = cosSig.*(2*cos(2*pi*fc*t));
  I_S = fft(I_sig);
  I_S_clip = [I_S(1:half_win) I_S(end-half_win+1:end)];
  I_sig_shallow = real(ifft(I_S_clip)*L/N);


  % Quadrature channel
  Q_sig = cosSig.*((-2)*sin(2*pi*fc*t));
  Q_S = fft(Q_sig);
  Q_S_clip = [Q_S(1:half_win) Q_S(end-half_win+1:end)];
  Q_sig_shallow = real(ifft(Q_S_clip)*L/N);

  % 合并成复信号
  pluralSig = I_sig_shallow + 1j*Q_sig_shallow;

end