%%% 扩展信号分析
% sig：输入需要分析的信号
% sample_freq：信号采样率
% mul_times：扩展信号的插值倍数
% win_length：样本窗长度，默认为64，且需输入偶数

function baseAnalyse(sig, sample_freq, mul_times, win_length)
  arguments 
    sig (1,:)
    sample_freq (1,1)
    mul_times (1,1)
    win_length (1,1) = 64
  end
  % 检查输入参数
  if (length(sig) < win_length)
    error('信号长度不能小于样本窗长度');
  end

  % 对信号进行窗截取
  [mainPeak, mLoc] = max(abs(sig));
  half = win_length/2;
  sig_win = sig(mLoc-half:mLoc+half-1);

  % 窗截取后扩展
  S_win = fft(sig_win);
  S_expend = [S_win(1:half), zeros(1, (mul_times-1)*win_length), S_win(half+1:end)];
  N = length(S_expend);
  sig_expend = ifft(S_expend)*N/win_length;

  %sig_win_uniform = 20*log10(abs(sig_win)/max(abs(sig_win)));
  sig_expend_uniform = 20*log10(abs(sig_expend)/max(abs(sig_expend)));
%   plot(sig_win_uniform);
%   hold on;
%   plot(linspace(1,win_length,N), sig_expend_uniform);
  plot(sig_expend_uniform, 'Color', [1 0.5 0], LineWidth=1);
  axis tight;ylim([-38, 0]);grid on;

end