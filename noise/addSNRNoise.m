function HSI_noisy = addSNRNoise(HSI_input, SNR_range)
    % 获取输入数据的尺寸
    [bands, N] = size(HSI_input);
    
    % 初始化加噪后的数据
    HSI_noisy = zeros(bands, N);
    
    % 对每个波段添加不同的SNR噪声
    for b = 1:bands
        % 在指定范围内随机生成一个SNR值
        SNR = SNR_range(1) + (SNR_range(2) - SNR_range(1)) * rand;
        
        % 计算当前波段的信号能量
        signal = HSI_input(b, :);
        signal_power = mean(signal.^2);
        
        % 计算噪声的标准差
        noise_power = signal_power / (10^(SNR / 10));
        noise_std = sqrt(noise_power);
        
        % 生成噪声并添加到当前波段
        noise = noise_std * randn(1, N);
        HSI_noisy(b, :) = signal + noise;
    end
end
