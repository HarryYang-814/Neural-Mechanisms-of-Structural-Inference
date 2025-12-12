
% -------------------------------------------------------------------------
function [pVal, crit_t] = perm_t_test(vec1, vec2, nPerm)
% Two‑sample permutation t‑test (unequal n)
    real_t = abs(t_stat(vec1,vec2));
    combo  = [vec1; vec2];
    n1 = numel(vec1);
    perm_t = zeros(nPerm,1);
    parfor k = 1:nPerm
        randIdx  = randperm(numel(combo));
        perm_t(k) = abs(t_stat(combo(randIdx(1:n1)), combo(randIdx(n1+1:end))));
    end
    pVal = mean(perm_t >= real_t);
    crit_t = prctile(perm_t, 95);
end

% -------------------------------------------------------------------------
function t = t_stat(a,b)
% Helper: independent‑samples t value (welch‑Satterthwaite)
    s1 = var(a); s2 = var(b);
    n1 = numel(a); n2 = numel(b);
    t  = (mean(a)-mean(b)) / sqrt(s1/n1 + s2/n2);
end
