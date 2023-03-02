function [y w B] = uniformly_weighted_beamformer(matX, a, N, L)

w = ones(N, 1)/N;
y = w'*matX;
for l = 1 : L
    for i = 1 : numel(a(1,:))
        B(i,l) = 20*log10(abs(w'*a(:,i)));
    end
end

end