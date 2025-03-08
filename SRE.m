function sre_value = SRE(x, x_hat)
    % 确保输入信号的维度一致
    if ~isequal(size(x), size(x_hat))
        error('原始信号和重构信号的尺寸不一致');
    end

    % 计算信号能量
    signal_energy = sum(x(:).^2);

    % 计算重构误差能量
    error_energy = sum((x(:) - x_hat(:)).^2);

    % 计算 SRE
    sre_value = 10 * log10(signal_energy / error_energy);
end