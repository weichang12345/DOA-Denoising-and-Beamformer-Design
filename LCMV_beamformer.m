function [y w B] = LCMV_beamformer(matX, P, theta_s_hat, theta_i_hat, d, a, sigma_o, mu, N, L)

a_s = exp(j*pi*sind(theta_s_hat).*d.');
a_i = exp(j*pi*sind(theta_i_hat).*d.');
P = eye(N)/sigma_o;
w = a_s/N;
y = zeros(1, L);

for l = 1 : L
    k = P*matX(:,l)/(mu + matX(:,l)'*P*matX(:,l));
    P = P/mu - k*matX(:,l)'*P/mu;
    C = [a_s(:,l), a_i(:,l)];
    g = [1; 10^-6];
    w = P*C*inv(C'*P*C)*g;
    y(l) = w'*matX(:,l);
    for i = 1 : numel(a(1,:))
        B(i,l) = 20*log10(abs(w'*a(:,i)));
    end
end

end