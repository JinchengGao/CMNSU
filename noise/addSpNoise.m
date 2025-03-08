function [X_noisy] = addSpNoise(X, Sp_rate)
    % 获取输入数据的尺寸
    [bands, N] = size(X);
    
    % 初始化加噪后的数据
    X_noisy = zeros(bands, N);
    
    % 对每个波段添加椒盐噪声
    for b = 1:bands
        % 提取当前波段的数据
        signal = X(b, :);
        
        % 创建一个与当前波段大小相同的矩阵
        Sp = 0.5 * ones(1, N);
        
        % 添加椒盐噪声
        Sp = imnoise(Sp, 'salt & pepper', Sp_rate);
        
        % 将椒盐噪声应用到当前波段
        signal(Sp == 0) = 0;
        signal(Sp == 1) = 1;
        
        % 将加噪后的波段数据存入输出矩阵
        X_noisy(b, :) = signal;
    end
end
