function B = retain_topk(A, k)
% RETAIN_TOPK 保留矩阵每列的前k个最大值，其余置零
%   输入：
%       A - n*n的相似度矩阵
%       k - 每列保留的最大值个数
%   输出：
%       B - 处理后的矩阵

% 输入验证
[n, m] = size(A);
assert(n == m, '输入矩阵必须是方阵');
assert(isscalar(k) && isnumeric(k), 'k必须是数值标量');

k = floor(k); % 确保整数
k = max(0, min(k, n)); % 限制k范围

% 处理特殊情况
if k == 0
    B = zeros(n);
    return;
elseif k >= n
    B = A;
    return;
end

% 核心算法实现
[~, sorted_idx] = sort(A, 1, 'descend'); % 列方向降序排序

% 构建掩膜矩阵
rows = sorted_idx(1:k, :); % 前k行索引
cols = repmat(1:n, k, 1);  % 列索引模板
linear_ind = sub2ind([n,n], rows(:), cols(:)); % 转换为线性索引

mask = false(size(A));
mask(linear_ind) = true;

B = A .* mask; % 应用掩膜
end