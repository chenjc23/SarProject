% CSA对radarsat数据成像处理
% by Jc
% 2022/12/03

clear, clc, close all;

%% 数据导入
load CDdata1.mat data;
load CD_run_params.mat ...
  Fr Kr PRF R0 Tr c f0 Nrg_cells Nrg_lines_blk;         % 成像处理相关的参数

Nr = Nrg_cells;                 % 距离门数目
Na = Nrg_lines_blk;             % 距离线数目

%% 原始数据补零
s_ra = zeros(Na*2, Nr*2);       % 方位维和距离维频域近似两倍的插值
s_ra(1:Na, 1:Nr) = data;
s_ra = double(s_ra);            % 转换成浮点数处理

[Nfft_a, Nfft_r] = size(s_ra);  % 二维fft点数

figure();
imagesc(abs(s_ra));
title('原始回波数据');

%% 其余参数设置（部分需要知道卫星工作参数）
lambda = c/f0;                  % 信号波长
Kr = -Kr;                       % 发射脉冲负扫频
Vr = 7062;                      % 直线几何约束下有效雷达速率（忽略空变，即为参考速度）
Fa = PRF;                       % 方位向采样率
fn_c = -6900;                   % 多普勒中心频率
%%%% 为减小补余的距离徙动，参考方位向频率应为测绘带中心处的多普勒中心频率 %%%%
fn_ref = fn_c;              


% 时间轴定义
tr = 2*R0/c + (-Nfft_r/2:Nfft_r/2-1)/Fr;          % 距离时间轴
ta = (-Nfft_a/2:Nfft_a/2-1)/Fa;                   % 方位时间轴

%  频率轴定义
fr = (-Nfft_r/2:Nfft_r/2-1)*Fr/Nfft_r;            % 距离频率轴
fa = fn_c + (-Nfft_a/2:Nfft_a/2-1)*Fa/Nfft_a;     % 方位频率轴
 


%% 距离多普勒域----->补余距离徙动校正

% 1. 在忽略速度的距离空变时，整体徙动因子与参考速度徙动因子相同 %
% 但方位频率不可忽视方位空变，故整体徙动因子与参考频率徙动因子不同 %
% 2. Kr与K_src均忽略距离空变，但K_src不可忽略方位空变 %
% 3. 变频时域的线性相位不改变幅值，不影响最终成像 %

%%%% 原始回波从方位中心频率变频至基带，防止方位模糊 %%%%
s_ra_2 = s_ra.*repmat(exp(-1j*2*pi*fn_c*ta).', 1, Nfft_r);    

%%%% 相关参数 %%%% 
R_ref = R0;                                           % 参考距离
D_rd_ref = sqrt(1-lambda^2*(fn_ref^2)/(4*Vr^2));      % rd域参考徙动因子
D_rd = sqrt(1-lambda^2*(fa.^2)/(4*Vr^2));             % rd域的整体徙动因子
K_src = 2*Vr^2*f0^3*D_rd.^3./(c*R0*fa.^2);            % 二次距离压缩对应的调频率
Km = Kr./(1-Kr./K_src);                               % rd域距离向调频率

%%%% 矩阵存储方式与meshgrid/imagesc有转置关系 %%%%
[tr_m, ta_m] = meshgrid(tr, ta);                  % 时间2D坐标网格 
[fr_m, fa_m] = meshgrid(fr, fa);                  % 频率2D坐标网格
D_rd_m = repmat(D_rd.', 1, Nfft_r);               % 矩阵形式的整体徙动因子
Km_m = repmat(Km.', 1, Nfft_r);                   % 矩阵形式的rd域距离向调频率


%%%% 平移距离向时间，平移项含距离方位耦合 %%%%
tau_shift_m = tr_m - 2*R_ref./(c*D_rd_m);         % 以参考点时间作距离时间轴平移（矩阵形式）

%%%% 构造变标方程 %%%%
s_sc = exp(1j*pi*Km_m.*(D_rd_ref./D_rd_m-1).*tau_shift_m.^2);

%%%% rd域信号与变标方程相乘 %%%%
S_rd = fft(s_ra_2, [], 1);            % 方位向fft至距离多普勒
S_rd_sc = S_rd.*s_sc;            % 补余距离徙动校正

figure();
imagesc(abs(S_rd_sc));

%% 二维频域----->距离脉压，SRC，一致RCMC

% 1. 距离脉压，SRC，一致RCMC可通过一个指数项进行相位补偿实现 %

%%%% 距离补偿信号：一次二次压缩同时完成 %%%%
s_r_com =  exp(1j*pi*D_rd_m./(Km_m*D_rd_ref).*fr_m.^2);

%%%% 一致距离徙动补偿信号 %%%%
s_RCMC_com = exp(1j*4*pi/c*(1./D_rd_m-1./D_rd_ref)*R_ref.*fr_m);

%%%% 统一进行补偿，补偿器需fftshift至零频在两端 %%%%
S_df = fft(S_rd_sc, [], 2);                                 % rd域距离向fft至二维频谱
S_df_2 = S_df .* fftshift(s_r_com.*s_RCMC_com, 2);          % 一次二次脉压+一致RCMC，完成所有距离处理

figure();
imagesc(abs(S_df_2));

%% 距离多普勒域----->方位脉压，附加相位补偿

% 1. 方位压缩时，R0的空变性对方位聚焦影响较大，不可忽略R0的空变性 %

R0_var = repmat(tr*c/2, Nfft_a, 1);

%%%% 方位匹配滤波器 %%%%
s_a_com = exp(1j*4*pi*R0_var/lambda .* D_rd_m);

%%%% 附加相位补偿器 %%%%
s_extra_com = exp(-1j*4*pi*Km_m/c^2 .* (1-D_rd_m/D_rd_ref) ...
              .* (R0_var./D_rd_m - R_ref./D_rd_m).^2);

%%%% 统一补偿 %%%%
S_rd_2 = ifft(S_df_2, [], 2);                     % 距离ifft至rd域
S_rd_3 = S_rd_2 .* s_a_com .* s_extra_com;        % 方位压缩，附加相位补偿

figure();
imagesc(abs(S_rd_3));

%% 二维图像域
s_ra_final = ifft(S_rd_3, [], 1);                 % 方位ifft至二维图像

figure();
imagesc(abs(s_ra_final));

%% 对亮度非线性变换，减小对比度
s_enhance = 20*log10(abs(s_ra_final)/max(max(abs(s_ra_final)))+eps);

figure();
imagesc(s_enhance, [-80, 0]);





