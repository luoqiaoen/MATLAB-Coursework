function [ fit ] = getDist(c, o)
c_n = size(c,2);
fit = norm(c(:, o(1)) - c(:,o(c_n)));
for i = 2:c_n
    fit = fit + norm(c(:, o(i)) - c(:,o(i -1)));
end
fit = 1/fit;
end

