function Gy=dyb(X)
% 使用边缘填充
padded_X = padarray(X, [1 1 0], 'replicate', 'both');

% 计算 Y 方向的梯度（非周期性边界）
Gy = padded_X(3:end, 2:end-1, :) - padded_X(2:end-1, 2:end-1, :);
end