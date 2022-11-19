%%% 
% sig：输入需要分析的信号
% sample_freq：信号采样率
% mul_times：扩展信号的插值倍数
% win_length：样本窗长度，默认为64

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
  half = round(win_length/2);
  sig_win = sig(mLoc-half:mLoc+half-1);

  % 窗截取后扩展
  S_win = fft(sig_win);
  S_expend = [S_win(1:half), zeros(1, mul_times-1), S_win(half+1:end)];
  sig_expend = fft(S_expend);

  plot(sig_win);
  hold on;
  plot(sig_expend);



  
end