%% SAR回波信号仿真（RD算法）
close all;
clear all, clc;
% 场景参数
c = 3e8;
range = 520;          % 距离向跨度
azimuth = 120;        % 方位向跨度
% 雷达参数  
theta_rc = 0;        % 雷达波束中心斜视角
%theta_rc = 0.9/180*pi;
% R_n0 = 10000;           % 雷达与场景中心距离，轨迹中点为二维零点，相应时刻为0时刻
R0 = 10000;              % 坐标中心横向距离(高度为0时为最短斜距)
R_nc = R0/cos(theta_rc);      % 波束中心穿越时刻对应的斜距

Vr = 100;             % 直线几何模型下雷达方位向运动速度
f0 = 10e9;            % 载频
Tr = 30e-6;           % 脉冲宽度

% 目标参数
obj_a_num = 1;      % 方位向目标个数
obj_r_num = 1;      % 距离向目标个数

loc_a = 0;
loc_r = R0;
% loc_a = linspace(-50, 50, obj_a_num);
% loc_r = R0 + linspace(-250, 250, obj_r_num);
[loc_r, loc_a] = meshgrid(loc_r, loc_a);        % 2D网格表征目标位置
loc_nc = loc_r/cos(theta_rc);                   % 波束中心穿越时刻时目标与雷达的斜距
nc = (loc_a-loc_r*tan(theta_rc))/Vr;            % 所有目标的波束中心穿越时刻

% 分辨率参数
rho_a = 1;     % 方位向分辨率
rho_r = 5;     % 距离向分辨率

% 其他参数定义
alpha = 0.886;                    % 1或0.886
alpha_osr = 1.4;                  % 距离过采样率
alpha_osa = 1.4;                  % 方位过采样率

lambda = c/f0;                                  % 信号波长
La = 2*rho_a;                                   % 天线长度
theta_bw = alpha*lambda/La;                     % 波束宽度
Ls = theta_bw*R_nc;                             % 合成孔径长度
Ta = Ls/Vr;                                     % 目标照射时间
f_dop = 2*Vr*cos(theta_rc)/lambda*theta_bw;     % 多普勒带宽
fnc = 2*Vr*sin(theta_rc)/lambda;                % 多普勒中心频率

Kr = c*alpha/(2*Tr*rho_r);                      % 距离向调频斜率
Ka = 2*Vr^2*cos(theta_rc)^2/lambda/R_nc;        % 方位向调频斜率
Fr = alpha_osr*Kr*Tr;                           % 距离向采样频率
Fa = alpha_osa*f_dop;                           % 方位向采样频率
%Fc = 6.8e8;                       
%Fr = Fc;

%% 生成区域回波信号
% 二维时间样点
tr = 2*(R0-range/2)/c-Tr/2:1/Fr:2*(R0+range/2)/c+Tr/2;   % 距离向时间
%ta = -azimuth/2/Vr-Ta:1/Fa:azimuth/2/Vr+Ta;             % 方位向时间
ta = 0:1/Fa:azimuth/2/Vr+Ta;
ta = [fliplr(-ta(2:end)) ta];                            % 方位向时间对准时间零点(考虑到方位向采样率低的影响)
%ta = fliplr(ta);      % 时间轴调转匹配频域变换得顺序
% lr = tr*c/2;        % 距离向距离
% la = ta*Vr;         % 方位向距离
Nr = length(tr); Na = length(ta);     % 距离向方位向采样点数
[t_r, t_a] = meshgrid(tr, ta);        % 距离向和方位向时间2D网格

% 二维时域信号构建
s0_tn = zeros(Na, Nr);                % 回波数据
for i=1:obj_a_num
  for j=1:obj_r_num
    Ls_diff = alpha*lambda/La*loc_nc(i,j);        % 目标照射长度(该点的合成孔径)
    Ts_diff = Ls_diff/(Vr*cos(theta_rc));                              % 目标照射时间
    Rn = sqrt(loc_r(i,j)^2+(loc_a(i,j)-Vr*t_a).^2);   % 目标的瞬时斜距，仅与方位时间有关
    %%% 距离向照射时间由脉宽决定，各距离相同 %%%
    %%% 而方位向照射时间由目标距离处合成孔径长度有关，实际上不同距离对应的照射时间不同 %%%
    % 距离向脉冲包络
    wr = abs(t_r-2*Rn/c)<Tr/2;                        
    % 方位向矩形方向图
    %wa = abs(t_a-nc(i,j))<Ts_diff/2;                       % 方位向方向图
    % 方位向sinc方向图
    theta = atan(Vr*(t_a-nc(i,j))/R_nc);              % 斜距平面内目标与雷达视线的夹角
    wa = sinc(alpha*theta/theta_bw).^2;                % 双程sinc方向图
    A0 = 1;     % 统一反射系数为1
    
    % 构建对应目标的原始基带信号回波
    s0_ij = A0*wr.*wa.*exp(-1j*4*pi*f0*Rn/c).*exp(1j*pi*Kr*(t_r-2*Rn/c).^2);
    
    % 载波调制的回波信号
    %s0_ij = A0*wr.*wa.*exp(1j*2*pi*f0*(t_r-2*Rn/c)).*exp(1j*pi*Kr*(t_r-2*Rn/c).^2);
    
    s0_tn = s0_tn + s0_ij;
  end
end

%s0_tn = s0_tn.*exp(-1j*2*pi*f0*(t_r));        % 复数域正交解调

% 以距离为尺度绘制二维时域
lr_sc = [tr(1), tr(end)]*c/2-R0;
la_sc = [ta(1), ta(end)]*Vr;
% figure(1);
% %mesh(t_r*c/2-R0, t_a*Vr, abs(s0_tn));
% imagesc(lr_sc, la_sc, abs(s0_tn));
% axis([-400, 400, -75,75]);
% % xlabel('距离向(相对场景中心)/m'), ylabel('方位向/m'), title('回波信号的二维时域(距离单位)');
% % colormap(summer), colorbar, shading flat;

% 以采样点为尺度绘制二维时域
figure(1);
% subplot(1,2,1);
% imagesc(abs(s0_tn));
% subplot(1,2,2);
imagesc(angle(s0_tn));

xlabel('距离向(采样点)'), ylabel('方位向(采样点)'), title('回波信号的二维时域');

%% 未压缩前的距离多普勒
Srd_uncompressed = fftshift(fft(s0_tn,[],1),1);
imagesc(angle(Srd_uncompressed));

S2df_uncompressed = fftshift(fft(Srd_uncompressed,[],2),2);
%imagesc(angle(S2df_uncompressed));

%% 距离向脉冲压缩
% 复制脉冲
tr_replica = -Tr/2:1/Fr:Tr/2-1/Fr;             % 复制发射脉冲时间样点
chirp_replica = exp(1j*pi*Kr*tr_replica.^2);   % 复制发射脉冲
Nc = length(chirp_replica);                    % 距离向chirp信号采样点数
Nf = Nr+Nc-1;                                  % 信号fft点数

% 复制脉冲时域补零至信号卷积的长度
chirp_extend = [chirp_replica, zeros(1, round((Nf-Nc)/2))];
chirp_extend = [zeros(1, Nf-length(chirp_extend)) chirp_extend];
chirp_extend = fftshift(chirp_extend);        % 翻折以对齐零时刻

% 回波距离向补零，进行距离匹配滤波
s0_extend = [s0_tn, zeros(Na, Nc-1)];               % 距离向补零，考虑频域相乘对应时域卷积
Sr_extend = fft(s0_extend,[],2);                    % 补零后求方位时间距离频谱
Sr_chirp = fft(chirp_extend);                       % 补零后chrip频谱
Sr_mul = Sr_extend .* repmat(conj(Sr_chirp),Na,1);  % 频域脉冲压缩
sr_compress = ifft(Sr_mul, [], 2);                 % 逆变换得距离压缩后回波的二维时域
sr_compress = sr_compress(:, 1:Nr);               % 舍去弃置区

figure();
%surf(t_r*c/2-R0, t_a*Vr, abs(sr_expressed));
%imagesc(lr_sc, la_sc, abs(sr_compress));
imagesc(abs(sr_compress));
xlabel('距离向(相对场景中心)/m'), ylabel('方位向/m'), title('二维时域距离向脉冲压缩(距离单位)');
colormap(turbo)
%axis([-400, 400, -75,75]);

%% 方位向傅里叶变换（距离多普勒域）
sr_compress_shift = sr_compress.*repmat(exp(-1j*2*pi*fnc*ta.'), 1, Nr); % 对应频谱搬移至0频
Sd = fft(sr_compress_shift,[],1);           % 方位向fft
Sd = fftshift(Sd, 1);                       % 频谱中心搬移至中心频率

%%% 此时多普勒已经转至零频，仅以多普勒中心频率构建多普勒域的频点范围作显示 %%%
f_Sd = 0:Fa/Na:Fa/2;
f_Sd = fnc + [fliplr(-f_Sd(2:end)) f_Sd];       

figure();
imagesc(lr_sc, [f_Sd(1) f_Sd(end)], abs(Sd));
title('距离多普勒域'), xlabel('距离向(相对场景中心)/m'), ylabel('多普勒频率/Hz')
colormap(turbo)

%% 距离徙动校正（插值法）
Rr = tr*(c/2);      % 每个距离门的对应距离
del_Rfn = lambda^2*(f_Sd.').^2*Rr/(8*Vr^2);         % 需要校正的距离徙动，随方位频率和距离变化

range_interval = (1/Fr)*c/2;                        % 距离单元间隔
del_Rfn_unit = del_Rfn/range_interval;              % 需要校正的距离间隔数

%%%%% 插值法距离徙动校正 %%%%% 
%%% 插值中的偏移量不进行量化，即不对sinc核深采样，直接使用浮点值，不存在量化的几何畸变， %%%
%%% 即不存在目标幅度调制，校正更准确。 %%%
core_len = 8;                 % sinc插值核长度
Sd_RCMC = zeros(Na, Nr);      % 存储RCMC后的距离多普勒 
for i = 1:Na
  for j = 1:Nr
    % 寻找校正对应点位
    idea_point = j+del_Rfn_unit(i, j);              % 需要校正到(i,j)位置的点(非量化点)
    idea_quant = round(idea_point);                 % 取整便于寻找位于距离门上的采样点
    % 确定信号和插值核的采样点位
    sig_sample_points = idea_quant + (-core_len/2:core_len/2-1);   % 信号上的插值采样点
    core_sample_points = idea_point - sig_sample_points;           % 插值核上的插值样点(浮点的)
    % 确定插值核
    sinc_core = sinc(core_sample_points);           % 插值系数
    sinc_core = sinc_core/sum(sinc_core);           % 插值系数归一化
    
    % 循环移位思想处理超出范围的信号插值采样点(算法简化版)
    exceed_loc = find(sig_sample_points>Nr);
    below_loc = find(sig_sample_points<1);
    sig_sample_points(exceed_loc) = sig_sample_points(exceed_loc) - Nr;
    sig_sample_points(below_loc) = sig_sample_points(below_loc) + Nr;
    % 最终校正步骤，插值加权求和
    Sd_RCMC(i, j) = sum(sinc_core.*Sd(i, sig_sample_points));
  end
end

figure();
imagesc(abs(Sd_RCMC));
title('距离多普勒域(RCMC校正后)'), xlabel('距离向(采样点)'), ylabel('方位向多普勒(采样点)')
colormap(turbo)

%% 方位向脉冲压缩（二维压缩结果）
ta_replica = -Ta/2:1/Fa:Ta/2-1/Fa;             % 复制方位脉冲时间样点
Nac = length(ta_replica);                      % chirp频谱样点数
Naf = Na+Nac-1;                                % 方位向fft点数

%%% 构建频域匹配滤波器 %%%
%%% 实际当Ka很小时，随中心时刻斜距的变化较明显，所以不能忽视距离向上中心时刻斜距的变化 %%%
R_nc_diff = tr*c/2/cos(theta_rc);
Ka_diff = 2*Vr^2*cos(theta_rc)^2/lambda./R_nc_diff; 
fd_extend = fftshift(linspace(-Fa/2, Fa/2, Naf));            
Haz_fn = exp(-1j*pi*(fd_extend.').^2*(1./Ka_diff));          % 匹配滤波器

%%% 压缩和RCMC校正后逆变换到时域补零做fft深采样 %%%
sra_compress_RCMC = ifft(fftshift(Sd_RCMC, 1));
Sd_RCMC_extend = fft(sra_compress_RCMC, Naf, 1);                 % 时域补零后的距离多普勒
%%% 方位向脉冲压缩 %%%
Sd_full_compress_RCMC = Sd_RCMC_extend .* (Haz_fn);              % 全压缩后的距离多普勒
sra_full_compress_RCMC = ifft(Sd_full_compress_RCMC, [], 1);     % 逆变换得RCMC和全压缩后的距离方位信号
sra_full_compress_RCMC = sra_full_compress_RCMC(1:Na, :);        % 舍去弃置区

figure();
imagesc(lr_sc, la_sc, abs(sra_full_compress_RCMC));
%mesh(abs(sa_expressed))
axis([-400, 400, -75,75]);
xlabel('距离向(相对场景中心)/m'), ylabel('方位向/m'), title('二维脉冲压缩结果(距离单位)');
colormap(turbo);

%% 二维频谱（最终之路！）
S_full_compress_RCMC = fft(Sd_full_compress_RCMC, [], 2);   % 距离向fft
S_full_compress_RCMC = fftshift(S_full_compress_RCMC, 2);
S_full_compress_RCMC = fftshift(S_full_compress_RCMC, 1);

figure();
imagesc(abs(S_full_compress_RCMC));















