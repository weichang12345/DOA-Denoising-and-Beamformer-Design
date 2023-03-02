function [y w B] = array_steering_beamformer(matX, N, theta_s_hat, a, d, L)

a_s = exp(j*pi*sind(theta_s_hat).*d.');
w = a_s/N;
y = zeros(1, numel(matX(1,:)));
for l = 1 : numel(matX(1,:))
    y(l) = w(:,l)'*matX(:,l);
    for i = 1 : numel(a(1,:))
        B(i,l) = 20*log10(abs(w(:,l)'*a(:,i)));
    end
end

end