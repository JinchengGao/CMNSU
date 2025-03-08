function HSI_striped = addStNoise(img, im_size)
    M = im_size(1);
    N = im_size(2);
    B = size(img, 1);  % 使用逗号分隔语法以避免错误
    img_noisy = permute(reshape(img, [B, M, N]), [2, 3, 1]);
    W1_GT = ones(size(img_noisy));  % 初始化W1_GT
    
    %% 设置随机数生成器种子
    k = 42;  % 这里设置一个初始值，可以根据需要调整
    
    %% 模拟对角条纹噪声
    rng(k);  k = k + 1;
    b_idx_rand = randperm(B);
    stripe_band_num = round(B * 0.3);  
    band2 = b_idx_rand(1:stripe_band_num);
    rng(k); k = k + 1;
    stripnum2 = 5 + ceil(10 * rand(1, length(band2))); 
   
    for i = 1:length(band2)
        for no_stripes = 1:stripnum2(i)
            rng(k);  k = k + 1;
            idx_start_r = randi(M);
            idx_start_c = 1; 
            t = 0;
            rng(k);  k = k + 1;
            signt = sign(rand(1) - 0.5);
            for ir = idx_start_r:M
                t = t + 1;
                if idx_start_c + t < N
                    rng(k);  k = k + 1;
                    t1 = signt * (rand(size(img_noisy(ir, idx_start_c + t, band2(i)))) * 0.2 + 0.5);
                    img_noisy(ir, idx_start_c + t, band2(i)) = 1;
                    W1_GT(ir, idx_start_c + t, band2(i)) = 0;
                end
            end
        end
        %% 模拟更宽的条纹
        t = 0;
        rng(k);  k = k + 1;
        idx_start_r = randi(M);
        rng(k);  k = k + 1;
        signt = sign(rand(1) - 0.5);
        for ir = idx_start_r:M
            t = t + 1;
            if idx_start_c + t < N
                rng(k);  k = k + 1;
                t1 = signt * (rand(size(img_noisy(ir, (idx_start_c + t):min(idx_start_c + t + 4, N), band2(i)))) * 0.2 + 0.5);
                img_noisy(ir, (idx_start_c + t):min(idx_start_c + t + 4, N), band2(i)) = 1;
                W1_GT(ir, (idx_start_c + t):min(idx_start_c + t + 4, N), band2(i)) = 0;
            end
        end
        
        for no_stripes = 1:fix(stripnum2(i) / 2)
            rng(k);  k = k + 1;
            idx_start_c = randi(N);
            idx_start_r = 1;
            t = 0;
            rng(k);  k = k + 1;
            signt = sign(rand(1) - 0.5);
            for ic = idx_start_c:N
                t = t + 1;
                if idx_start_r + t < M
                    rng(k);  k = k + 1;
                    t1 = signt * (rand(size(img_noisy(idx_start_r + t, ic, band2(i)))) * 0.2 + 0.5);
                    img_noisy(idx_start_r + t, ic, band2(i)) = 1;
                    W1_GT(idx_start_r + t, ic, band2(i)) = 0;
                end
            end
        end
    end
    
    %% 输出加噪后的图像
    HSI_striped = permute(img_noisy, [3, 1, 2]);  % 恢复为原来的维度顺序
    HSI_striped = reshape(HSI_striped, [B, M * N]);  % 重新reshape为二维图像
    
end
