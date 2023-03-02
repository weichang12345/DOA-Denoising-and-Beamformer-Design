clear;

data = load('ASP_Final_Data');
theta_s_noisy = data.theta_s_noisy;
theta_i_noisy = data.theta_i_noisy;
matX = data.matX;
L = numel(theta_s_noisy);  % Time length
t = [1:L];                 % Time vector
N = numel(matX(:,1));      % Number of isotropic antennas
d = 0:(N-1);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
%% Denoising by method of EMD

imf_s = hht(theta_s_noisy, t, 0.2);
zero_imf_s = 2000*zerocrossrate(imf_s.',"Method","comparison");  % The zero crossings of each IMFs, can be treated as frequency of each IMFs
peak_s = islocalmax(theta_s_noisy,'MinProminence',2);            % Prominence peaks of theta_s_noisy
dip_s = islocalmin(theta_s_noisy,'MinProminence',2);             % Prominence peaks of theta_s_noisy
extreme_s = (sum(peak_s) + sum(dip_s))/2;                        % The average of prominence extremes of theta_s_noisy, can be treated as frequency of theta_s_noisy
range_s = sum(zero_imf_s < extreme_s);                           % Determine the number of IMFs to reconstruct the signal. The IMFs are adopted until the first IMF's frequency is larger than      ` theta_s_noisy's frequency
theta_s_hat = sum(imf_s(end-range_s:end,:));                     % Reconstruct the signal, and the result will be theta_s_hat

imf_i = hht(theta_i_noisy, t, 0.2);
zero_imf_i = 2000*zerocrossrate(imf_i.',"Method","comparison");  % The zero crossings of each IMFs, can be treated as frequency of each IMFs
peak_i = islocalmax(theta_i_noisy,'MinProminence',2);            % Prominence peaks of theta_i_noisy
dip_i = islocalmin(theta_i_noisy,'MinProminence',2);             % Prominence peaks of theta_i_noisy
extreme_i = (sum(peak_i) + sum(dip_i))/2;                        % The average of prominence extremes of theta_i_noisy, can be treated as frequency of theta_i_noisy
range_i = sum(zero_imf_i < extreme_i);                           % Determine the number of IMFs to reconstruct the signal. The IMFs are adopted until the first IMF's frequency is larger than theta_i_noisy's frequency
theta_i_hat = sum(imf_i(end-range_i:end,:));                     % Reconstruct the signal, and the result will be theta_i_hat

figure;
plot(t,theta_s_noisy, t,theta_i_noisy);
ylim([-10 20]);
legend('Signal $\tilde{\theta}_{s}(t)$','Interference $\tilde{\theta}_{i}(t)$','interpreter','latex','fontsize',14)
xlabel('Time')
ylabel('Degree')
title('DOA Without Denoised')

figure;
plot(t,theta_s_hat,t,theta_i_hat);
ylim([-10 20]);
legend('Signal $\hat{\theta}_{s}(t)$','Interference $\hat{\theta}_{i}(t)$','interpreter','latex','fontsize',14)
xlabel('Time')
ylabel('Degree')
title('DOA After Denoised')

%% Denoising by method of MVDR Spectrum

theta = linspace(-pi/2, pi/2, 500);
theta_deg = rad2deg(theta);
a = exp(j*pi*sin(theta).*d.');

delta = 0.01;
mu = 0.53;
k = zeros(N,1);
P = eye(N)/delta;
for l = 1 : L
    k = P*matX(:,l)/(mu + matX(:,l)'*P*matX(:,l));
    P = P/mu - k*matX(:,l)'*P/mu;
    for i = 1 : numel(a(1,:))
        MVDR_spectrum(i,l) = 1/(a(:,i)'*P*a(:,i));
    end
end

figure;  % Display MVDR Spectrum in image form
image(t, theta, abs(MVDR_spectrum)/max(max(abs(MVDR_spectrum)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree (Radians)')
title('MVDR Spectrum')

[value_s, index_s] = max(MVDR_spectrum(1:ceil(numel(theta)/2),:));
[value_i, index_i] = max(MVDR_spectrum(ceil(numel(theta)/2) + 1:end,:));
index_i = index_i + ceil(numel(theta)/2) - 1;

figure;  % Plot MVDR Spectrum with respect of Cartesian coordinate, and compare it with EMD denoising method
plot(t,index_s*180/(length(theta)-1)-90,t,theta_s_hat,t,index_i*180/(length(theta)-1)-90,t,theta_i_hat)
legend('MVDR Spectrum $\hat{\theta}_{s}(t)$','EMD denoising $\hat{\theta}_{s}(t)$',...
    ' MVDR Spectrum $\hat{\theta}_{i}(t)$','EMD denoising $\hat{\theta}_{i}(t)$','interpreter','latex','fontsize',10)
xlabel('Time')
ylabel('Degree')
title('Comparison between EMD denoising and MVDR Spectrum')

%% Estimate source signal by uniformly weighted beamformer

[y_uniform w_uniform B_uniform] = uniformly_weighted_beamformer(matX, a, N, L);

figure;
image(t, theta_deg, abs(B_uniform)/max(max(abs(B_uniform)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree')
title('Beampattern of uniformly weighted beamformer')

figure;
subplot(2,1,1)
plot(t,real(y_uniform));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by uniformly weight beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(y_uniform));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by uniformly weight beamformer','interpreter','latex','fontsize',18)

%% Estimate source signal by array steering beamformer

[y_array_steering w_array_steering B_array_steering] = array_steering_beamformer(matX, N, theta_s_hat, a, d, L);

figure;
image(t, theta_deg, abs(B_array_steering)/max(max(abs(B_array_steering)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree')
title('Beampattern of array steering beamformer')

figure;
subplot(2,1,1)
plot(t,real(y_array_steering));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by array steering beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(y_array_steering));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by array steering beamformer','interpreter','latex','fontsize',18)

%% Estimate source signal by MVDR beamformer
a_s = exp(j*pi*sind(theta_s_hat).*d.');
a_i = exp(j*pi*sind(theta_i_hat).*d.');

[y_MVDR w_MVDR B_MVDR]= MVDR_beamformer(matX, P, theta_s_hat, d, a, 0.01, 0.999, N, L);

figure;
image(t, theta_deg, abs(B_MVDR)/max(max(abs(B_MVDR)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree')
title('Beampattern of MVDR beamformer')

figure;
subplot(2,1,1)
plot(t,real(y_MVDR));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by MVDR beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(y_MVDR));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by MVDR beamformer','interpreter','latex','fontsize',18)

%% Estimate source signal by LCMV beamformer

[y_LCMV w_LCMV B_LCMV] = LCMV_beamformer(matX, P, theta_s_hat, theta_i_hat, d, a, 0.01, 0.999, N, L);

figure;
image(t, theta_deg, abs(B_LCMV)/max(max(abs(B_LCMV)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree')
title('Beampattern of LCMV beamformer')

figure;
subplot(2,1,1)
plot(t,real(y_LCMV));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by LCMV beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(y_LCMV));
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by LCMV beamformer','interpreter','latex','fontsize',18)

%% Estimate source signal by proposed beamformer

[s_t_hat w B] = proposed_beamformer(matX, N, theta_s_hat, theta_i_hat, a, d, w_array_steering, L);

figure;
image(t, theta_deg, abs(B)/max(max(abs(B)))*1000);
set(gca, "Ydir", "normal");
xlabel('Time')
ylabel('Degree')
title('Beampattern of proposed beamformer')

figure;
subplot(2,1,1)
plot(t,real(s_t_hat))
axis tight
xlabel('Time','fontsize',16)
ylabel('$Re\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Real part of $\hat{s}(t)$ by proposed beamformer','interpreter','latex','fontsize',18)
subplot(2,1,2)
plot(t,imag(s_t_hat))
axis tight
xlabel('Time','fontsize',16)
ylabel('$Im\left\{\hat{s}(t)\right\}$','interpreter','latex','fontsize',16)
title('Imaginary part of $\hat{s}(t)$ by proposed beamformer','interpreter','latex','fontsize',18)

%%

save('r11942140_張聿維_ASPFinal_PerformanceEvaluation.mat');
print -depsc myfig.eps
