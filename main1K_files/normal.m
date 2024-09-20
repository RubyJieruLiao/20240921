function w = normal(W)
    [N, ~] = size(W);
    u = zeros(N * (N - 1) / 2, 1);
    k = 0;
    for i = 1:N
        for j = i + 1:N
            k = k + 1;
            u(k, 1) = W(i, j);
        end
    end
    u = sort(u);
    uu = u(floor(length(u) * 0.95));
    w = W - diag(diag(W));
    idx = w > uu;
    w(idx) = uu;
    w = w / (max(max(w)));
end