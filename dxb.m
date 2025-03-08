function Gx=dxb(X)
% 使用边缘填充
padded_X = padarray(X, [1 1 0], 'replicate', 'both');

% 计算 X 方向的梯度（非周期性边界）
Gx = padded_X(2:end-1, 3:end, :) - padded_X(2:end-1, 2:end-1, :);
end