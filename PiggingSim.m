function PiggingSim()
% =========================================================================
%  清管器清蜡运动规律模拟 - 单文件版
%  参考：李苗. 原油管道清管器运动规律研究[D]. 中国石油大学(北京), 2018.
%
%  功能：
%    模拟含蜡原油管道清管过程中清管器的速度、位置、管道沿线压力、
%    入口压力、出口流量、球前浆体长度和蜡浓度随时间的变化。
%
%  使用方法：
%    将本文件保存为 PiggingSim.m，在 MATLAB 命令行输入 PiggingSim 回车。
%
%  模块结构（全部为局部函数）：
%    1. get_parameters      - 参数输入
%    2. initialize_state    - 状态初始化
%    3. compute_forces      - 蜡层阻力 + 管壁摩擦力
%    4. locate_boundaries   - 清管器/浆体前缘边界定位
%    5. solve_fluid         - 一维瞬态流体有限差分求解
%    6. solve_pig_motion    - 清管器运动方程
%    7. update_slurry       - 浆体增长与物性更新
%    8. check_stuck         - 卡堵判断
%    9. record_history      - 历史数据记录
%   10. plot_results        - 结果绘图
% =========================================================================

clc; close all;
fprintf('===============================================\n');
fprintf('  清管器清蜡运动规律模拟\n');
fprintf('  基于李苗博士论文 (中国石油大学, 2018)\n');
fprintf('===============================================\n\n');

% ---------- 1. 参数输入 ----------
P = get_parameters();

% ---------- 2. 初始化 ----------
S = initialize_state(P);

% ---------- 3. 时间推进 ----------
fprintf('开始模拟...\n');
fprintf('管道长度 %.1f km | 管径 %.3f m | 沉积物厚度 %.3f m\n', ...
    P.L/1000, P.D, P.h_wax);
fprintf('时间步长 %.3f s | 空间步长 %.1f m | 总时长 %.0f s\n\n', ...
    P.dt, P.dx, P.t_end);

tic;
while S.t < P.t_end && ~S.stopped
    S = compute_forces(S, P);
    S = locate_boundaries(S, P);
    S = solve_fluid(S, P);
    S = solve_pig_motion(S, P);
    S = update_slurry(S, P);
    S = check_stuck(S, P);
    S = record_history(S, P);

    S.t    = S.t + P.dt;
    S.step = S.step + 1;

    if mod(S.step, P.print_interval) == 0
        fprintf(['t=%7.1fs | x_pig=%7.1fm | v_pig=%.3fm/s | ' ...
                 'P_in=%.3fMPa | L_s=%6.1fm | w=%.3f\n'], ...
            S.t, S.x_pig, S.v_pig, S.P(1)/1e6, S.L_s, S.w_slurry);
    end
end

if S.stopped
    fprintf('\n*** 模拟终止：%s ***\n', S.stuck_reason);
else
    fprintf('\n*** 模拟到达设定总时间 ***\n');
end
fprintf('总时间 %.1f s，共 %d 步\n', S.t, S.step);
toc;

% ---------- 4. 后处理 ----------
plot_results(S, P);
end


% =========================================================================
%  模块 1：参数输入
% =========================================================================
function P = get_parameters()
P = struct();

% ---- 管道参数 ----
P.L    = 30000;        % 管道长度 (m)
P.D    = 0.704;        % 管道内径 (m)
P.beta = 0;            % 管道倾角 (rad)，水平管取0
P.g    = 9.81;         % 重力加速度 (m/s^2)

% ---- 油品参数 ----
P.rho_oil = 850;                      % 密度 (kg/m^3)
P.mu_oil  = 30e-3;                    % 动力黏度 (Pa·s)
P.K_oil   = 1.5e9;                    % 体积模量 (Pa)
P.c_oil   = sqrt(P.K_oil/P.rho_oil);  % 声速 (m/s)

% ---- 蜡沉积物参数 ----
P.h_wax   = 0.05;      % 沉积物厚度 (m)
P.rho_wax = 900;       % 蜡沉积物密度 (kg/m^3)
P.tau_y   = 5;      % 屈服应力 (Pa)

% ---- 清管器参数 ----
P.m_pig    = 100;      % 质量 (kg)
P.oversize = 0.005;     % 过盈量 (2%)
P.D_v      = 0.1;     % 泄流孔直径 (m)，0 表示无泄流孔
P.K_f      = 2.5;      % 泄流孔局部阻力系数
P.E_clean  = 0.9;      % 清蜡效率

% ---- 泵与边界 ----
P.P0_pump = 3.5e6;     % 泵关死压力 (Pa)
P.a_pump  = 1.0e5;     % 泵特性系数 (Pa/(m^3/s)^2)
P.P_out   = 1.0e6;     % 出口压力 (Pa)
P.P_max   = 4.0e6;     % 最大允许压力 (Pa)
P.Q_min   = 0.5;       % 最小入口流速 (m/s)

% ---- 数值参数 ----
P.dx   = 500;                      % 空间步长 (m)
P.N    = round(P.L/P.dx) + 1;      % 节点数
P.dt   = 0.1;                      % 时间步长 (s)
P.t_end = 20000;                   % 总模拟时间 (s)
P.print_interval   = 2000;         % 屏幕输出间隔（步）
P.profile_interval = 5000;         % 保存压力剖面间隔（步）

% ---- 检查 CFL 条件 ----
CFL = P.c_oil * P.dt / P.dx;
if CFL > 1
    warning('CFL=%.3f > 1，计算可能不稳定，建议减小 dt 或增大 dx', CFL);
end
end


% =========================================================================
%  模块 2：状态初始化
% =========================================================================
function S = initialize_state(P)
S = struct();

% ---- 网格 ----
S.x = linspace(0, P.L, P.N)';
S.N = P.N;

% ---- 初始压力/速度场 ----
S.P = linspace(P.P0_pump, P.P_out, P.N)';
S.U = 1.4 * ones(P.N, 1);

% ---- 清管器 ----
S.x_pig     = 1.0;
S.v_pig     = 1.4;
S.i_pig     = 2;
S.theta_pig = 0;

% ---- 浆体 ----
S.L_s            = 0;
S.x_slurry_front = S.x_pig;
S.w_slurry       = 0;
S.eta_slurry     = P.mu_oil;
S.rho_slurry     = P.rho_oil;
S.i_sf           = 2;
S.theta_sf       = 0;

% ---- 阻力 ----
S.F_f  = 0;
S.F_w  = 0;
S.dP_w = 0;

% ---- 时间与状态 ----
S.t            = 0;
S.step         = 0;
S.stopped      = false;
S.stuck_reason = '';

% ---- 历史记录 ----
S.hist_t       = [];
S.hist_x_pig   = [];
S.hist_v_pig   = [];
S.hist_P_in    = [];
S.hist_Q_out   = [];
S.hist_L_s     = [];
S.hist_w       = [];
S.hist_P_profile   = {};
S.hist_P_profile_t = [];
end


% =========================================================================
%  模块 3：蜡层阻力 + 管壁摩擦力
% =========================================================================
function S = compute_forces(S, P)
% ---- 管壁摩擦力：F_f = 1140*ln(delta) + 2019，delta 单位 % ----
delta_pct = P.oversize * 100;
if delta_pct > 0
    S.F_f = 1140*log(delta_pct) + 2019;
else
    S.F_f = 0;
end
S.F_f = max(S.F_f, 0);

% ---- 蜡层破坏压差：dP_w = 0.74*(h/D)^(-0.17)*tau_y ----
h = P.h_wax;
D = P.D;
if h > 0
    S.dP_w = 0.74 * (h/D)^(-0.17) * P.tau_y;
else
    S.dP_w = 0;
end
A     = pi * D^2 / 4;
S.F_w = S.dP_w * A;
end


% =========================================================================
%  模块 4：清管器/浆体前缘边界定位
% =========================================================================
function S = locate_boundaries(S, P)
% ---- 清管器位置 ----
i = floor(S.x_pig / P.dx) + 1;
i = min(max(i, 1), P.N-1);
S.i_pig     = i;
S.theta_pig = (S.x_pig - S.x(i)) / P.dx;
S.theta_pig = min(max(S.theta_pig, 0), 1);

% ---- 浆体前缘位置 ----
S.x_slurry_front = S.x_pig + S.L_s;
j = floor(S.x_slurry_front / P.dx) + 1;
j = min(max(j, 1), P.N-1);
S.i_sf     = j;
S.theta_sf = (S.x_slurry_front - S.x(j)) / P.dx;
S.theta_sf = min(max(S.theta_sf, 0), 1);
end


% =========================================================================
%  模块 5：一维瞬态流体有限差分求解
% =========================================================================
function S = solve_fluid(S, P)
N     = S.N;
dx    = P.dx;
dt    = P.dt;

P_old = S.P;
U_old = S.U;
P_new = P_old;
U_new = U_old;

i_pig = S.i_pig;
i_sf  = S.i_sf;

% ---------- 内部节点 ----------
for i = 2:N-1
    if i < i_pig
        rho = P.rho_oil;    mu = P.mu_oil;    c = P.c_oil;
    elseif i > i_pig && i <= i_sf
        rho = S.rho_slurry; mu = max(S.eta_slurry, 1e-6);
        c   = sqrt(1.5e9 / max(rho, 1));
    else
        rho = P.rho_oil;    mu = P.mu_oil;    c = P.c_oil;
    end

    Ui = U_old(i);
    Re = rho * abs(Ui) * P.D / mu;
    if Re < 1
        f = 64;
    elseif Re < 3000
        f = 64 / Re;
    else
        f = 0.3164 * Re^(-0.25);
    end

    if Ui >= 0
        dPdx = (P_old(i) - P_old(i-1)) / dx;
        dUdx = (U_old(i) - U_old(i-1)) / dx;
    else
        dPdx = (P_old(i+1) - P_old(i)) / dx;
        dUdx = (U_old(i+1) - U_old(i)) / dx;
    end

    dPdt = -Ui * dPdx - rho * c^2 * dUdx;
    dUdt = -Ui * dUdx - dPdx / rho ...
           - f * Ui * abs(Ui) / (2*P.D) - P.g*sin(P.beta);

    P_new(i) = P_old(i) + dPdt * dt;
    U_new(i) = U_old(i) + dUdt * dt;
end

% ---------- 边界条件 ----------
A = pi * P.D^2 / 4;

% 入口：泵特性 P = P0 - a*Q^2
Q_in     = max(U_new(2), 0) * A;
P_new(1) = P.P0_pump - P.a_pump * Q_in^2;
P_new(1) = max(P_new(1), P.P_out);
U_new(1) = max(U_new(2), 0);

% 出口：定压
P_new(N) = P.P_out;
U_new(N) = U_new(N-1);

% 清管器尾部
if i_pig >= 1 && i_pig <= N
    U_new(i_pig) = S.v_pig;
    if i_pig > 1
        P_new(i_pig) = P_old(i_pig-1) ...
            - P.rho_oil * P.c_oil * (S.v_pig - U_old(i_pig-1));
        P_new(i_pig) = max(P_new(i_pig), 0);
    end
end

% 清管器头部：有/无泄流孔
if i_pig + 1 <= N
    if P.D_v > 0
        U_v    = S.v_pig * (P.D / P.D_v)^2;
        U_nose = S.v_pig - (P.D_v/P.D)^2 * (U_v - S.v_pig);
        dP_byp = 0.5 * S.rho_slurry * P.K_f * (U_v - S.v_pig)^2;
        U_new(i_pig+1) = U_nose;
        P_new(i_pig+1) = max(P_new(i_pig) - dP_byp, 0);
    else
        U_new(i_pig+1) = S.v_pig;
        P_new(i_pig+1) = max(P_new(i_pig) - S.dP_w, 0);
    end
end

% 浆体前缘：P、U 连续
if i_sf >= 2 && i_sf <= N-1
    U_new(i_sf) = U_new(i_sf-1);
    P_new(i_sf) = P_new(i_sf-1);
end

% 数值保护
P_new(~isfinite(P_new)) = P.P_out;
U_new(~isfinite(U_new)) = 0;
U_new = max(U_new, 0);
P_new = max(P_new, 0);

S.P = P_new;
S.U = U_new;
end


% =========================================================================
%  模块 6：清管器运动方程
% =========================================================================
function S = solve_pig_motion(S, P)
A  = pi * P.D^2 / 4;
i  = S.i_pig;
th = S.theta_pig;

% 插值求尾部压力
if i < P.N
    P_tail = (1-th) * S.P(i) + th * S.P(i+1);
else
    P_tail = S.P(i);
end
P_nose = P_tail - S.dP_w;

% 加速度
a_pig = ((P_tail - P_nose)*A - S.F_f - S.F_w ...
         - P.m_pig * P.g * sin(P.beta)) / P.m_pig;

% 更新速度、位置
S.v_pig = max(S.v_pig + a_pig * P.dt, 0);
S.x_pig = S.x_pig + S.v_pig * P.dt;

% 到达管末
if S.x_pig >= P.L
    S.x_pig = P.L;
    S.stopped = true;
    S.stuck_reason = '清管器到达管道末端';
end
end


% =========================================================================
%  模块 7：浆体增长与物性更新
% =========================================================================
function S = update_slurry(S, P)
% ---- 蜡浓度 ----
if P.D_v > 0
    U_v = S.v_pig * (P.D / P.D_v)^2;
    m   = 1.39;
    S.w_slurry = 1 / (1 + max(U_v / max(S.v_pig, 0.01), 0)^m);
    S.w_slurry = min(max(S.w_slurry, 0), 1);
else
    S.w_slurry = 1.0;
end

% ---- 浆体密度 ----
S.rho_slurry = P.rho_wax * S.w_slurry + P.rho_oil * (1 - S.w_slurry);

% ---- 浆体黏度（论文式，w 以 wt% 计）----
w_pct = S.w_slurry * 100;
S.eta_slurry = 59.1e-3 * exp(0.32 * w_pct);

% ---- 浆体增长 ----
h = P.h_wax;  D = P.D;  E = P.E_clean;
if S.v_pig > 1e-3 && S.w_slurry > 0.01
    dLsdt = (P.rho_wax / S.rho_slurry) ...
            * (4*h*(D-h) / D^2) * (E / S.w_slurry) * S.v_pig;
    S.L_s = S.L_s + dLsdt * P.dt;
end

S.x_slurry_front = S.x_pig + S.L_s;
end


% =========================================================================
%  模块 8：卡堵判断
% =========================================================================
function S = check_stuck(S, P)
if S.P(1) > P.P_max
    S.stopped = true;
    S.stuck_reason = sprintf('入口压力 %.2f MPa 超过最大允许压力 %.2f MPa', ...
        S.P(1)/1e6, P.P_max/1e6);
    return;
end

if S.t > 10 && S.v_pig < 0.01
    S.stopped = true;
    S.stuck_reason = sprintf('清管器速度趋零 (v=%.4f m/s)', S.v_pig);
    return;
end

if S.t > 10 && S.U(1) < P.Q_min
    S.stopped = true;
    S.stuck_reason = sprintf('入口流速 %.4f m/s 低于最小值 %.4f m/s', ...
        S.U(1), P.Q_min);
    return;
end
end


% =========================================================================
%  模块 9：历史数据记录
% =========================================================================
function S = record_history(S, P)
S.hist_t(end+1)     = S.t;
S.hist_x_pig(end+1) = S.x_pig;
S.hist_v_pig(end+1) = S.v_pig;
S.hist_P_in(end+1)  = S.P(1);
S.hist_Q_out(end+1) = S.U(end);
S.hist_L_s(end+1)   = S.L_s;
S.hist_w(end+1)     = S.w_slurry;

if mod(S.step, P.profile_interval) == 0
    S.hist_P_profile{end+1}   = S.P;
    S.hist_P_profile_t(end+1) = S.t;
end
end


% =========================================================================
%  模块 10：结果绘图
% =========================================================================
function plot_results(S, P)
figure('Position', [80, 80, 1200, 760], 'Name', '清管器清蜡模拟结果');

subplot(2,3,1);
plot(S.hist_t, S.hist_v_pig, 'b-', 'LineWidth', 1.4);
xlabel('时间 t (s)'); ylabel('清管器速度 (m/s)');
title('清管器运行速度'); grid on;

subplot(2,3,2);
plot(S.hist_t, S.hist_x_pig, 'r-', 'LineWidth', 1.4);
xlabel('时间 t (s)'); ylabel('清管器位置 (m)');
title('清管器位置'); grid on;

subplot(2,3,3);
yyaxis left;
plot(S.hist_t, S.hist_P_in/1e6, 'b-', 'LineWidth', 1.4);
ylabel('入口压力 (MPa)');
yyaxis right;
plot(S.hist_t, S.hist_Q_out, 'r-', 'LineWidth', 1.4);
ylabel('出口流速 (m/s)');
xlabel('时间 t (s)');
title('入口压力与出口流速'); grid on;

subplot(2,3,4);
plot(S.hist_t, S.hist_L_s, 'g-', 'LineWidth', 1.4);
xlabel('时间 t (s)'); ylabel('浆体长度 (m)');
title('浆体长度'); grid on;

subplot(2,3,5);
plot(S.hist_t, S.hist_w, 'm-', 'LineWidth', 1.4);
xlabel('时间 t (s)'); ylabel('浆体蜡浓度 w');
title('浆体蜡浓度'); grid on;

subplot(2,3,6);
hold on;
n = min(5, numel(S.hist_P_profile));
if n > 0
    idx = round(linspace(1, numel(S.hist_P_profile), n));
    for k = idx
        plot(S.x/1000, S.hist_P_profile{k}/1e6, 'LineWidth', 1.3, ...
            'DisplayName', sprintf('t=%.0fs', S.hist_P_profile_t(k)));
    end
end
xlabel('里程 (km)'); ylabel('压力 (MPa)');
title('管道沿线压力分布');
legend('Location','best'); grid on;
end