function X = CMN_TV_GPU(E,y,varargin)
q = 2;
p_f = 2;
p_s = 1;
gamma = 5;
mu1 = 1e-1;
mu2 = 1e-1;
step = 2;
Iter = 200;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'AL_ITERS'
                Iter = varargin{i+1};
            case 'MU_1'
                mu1 = varargin{i+1};
            case 'MU_2'
                mu2 = varargin{i+1};
            case 'GAMMA'
                gamma = varargin{i+1};
            case 'PF'
                p_f = varargin{i+1};
            case 'PS'
                p_s = varargin{i+1};
            case 'Q'
                q = varargin{i+1};
            case 'STEP'
                step = varargin{i+1};
            case 'IM_SIZE'
                im_size = varargin{i+1};
                n_row = im_size(1);
                n_col = im_size(2);
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end
[~, R] = size(E); % E represents the endmember matrix, R is the number of endmembers
[B, N] = size(y); % y represents the input data, N is the number of samples
X = gpuArray(zeros(R, N));
V1 = gpuArray(zeros(R, N));
V2 = gpuArray(zeros(R, N));
V3 = gpuArray(zeros(B, N));
V4 = gpuArray(zeros(R, N));
V5x = gpuArray(zeros(n_row,n_col,R));
V5y = gpuArray(zeros(n_row,n_col,R));
Lambda1 = gpuArray(zeros(R,N));
Lambda2 = gpuArray(zeros(R,N));
Lambda3 = gpuArray(zeros(B,N));
Lambda4 = gpuArray(zeros(R,N));
Lambda5x = gpuArray(zeros(n_row,n_col,R));
Lambda5y = gpuArray(zeros(n_row,n_col,R));
tol = 1e-4;
tol1 = sqrt(N * R) * tol;
tol2 = sqrt(N * R) * tol;
i = 1;
res_p = inf * ones(1,5);
res_d = inf * ones(1,5);
res_p_total = inf;
res_d_total = inf;

FDh = gpuArray(zeros(n_row,n_col));
FDh(1,1) = -1;
FDh(1,end) = 1;
FDh = fft2(FDh);
FDhH = conj(FDh);

FDv = gpuArray(zeros(n_row,n_col));
FDv(1,1) = -1;
FDv(end,1) = 1;
FDv = fft2(FDv);
FDvH = conj(FDv);

IL = 1./( FDhH.* FDh + FDvH.* FDv + 3);

Dh = @(x) real(ifft2(fft2(x).*FDh));
DhH = @(x) real(ifft2(fft2(x).*FDhH));
Dv = @(x) real(ifft2(fft2(x).*FDv));
DvH = @(x) real(ifft2(fft2(x).*FDvH));

%% 开始迭代
while (i <= Iter) && ((abs(res_p_total) >= tol1) || (abs(res_d_total) >= tol2))
    %%  update X
    temp1 = permute(reshape(V1 + V2 + V4 + (Lambda1 + Lambda2 + Lambda4)/gamma,[R,n_row,n_col]),[2,3,1]);
    temp2 = DhH(V5x + Lambda5x./gamma);
    temp3 = DvH(V5y + Lambda5y./gamma);
    X_reshape = permute(reshape(X,[R,n_row,n_col]),[2,3,1]);
    for r=1:R
        numerator_ft = fftn(temp1(:,:,r) + temp2(:,:,r) + temp3(:,:,r));
        X_reshape(:,:,r) = real(ifftn(numerator_ft.*IL));
    end
    X = permute(reshape(X_reshape,[n_row*n_col,R]),[2,1]);
    %% update V1
    V1_old = V1;
    V1 = soft_threshold(X - Lambda1./gamma , mu1./gamma);

    %% update V2
    V2_old = V2;
    V2 = (E'*E + eye(R)) \ (X + E'*V3 + E'*Lambda3./gamma - Lambda2 ./ gamma + E' * y);

    %% update V3
    V3_old = V3;
    V3 = Soft1(p_f, p_s, q, V2, V3_old, E,y,Lambda3, gamma);

    %% update V4
    V4_old = V4;
    V4 = max(X - Lambda4./gamma ,0);

    %% update V5
    V5x_old = V5x;
    V5y_old  =V5y;
    V5x =  soft_threshold(Dh(X_reshape) - Lambda5x./gamma,mu2./gamma);
    V5y =  soft_threshold(Dv(X_reshape) - Lambda5y./gamma,mu2./gamma);
    %% update Lambda1
    Lambda1 = Lambda1 + gamma*(V1 - X);
    res_p(1) = norm(V1-X,'fro');
    res_d(1) = norm(V1-V1_old,'fro');
    %% update Lambda2
    Lambda2 = Lambda2 + gamma*(V2 - X);
    res_p(2) = norm(V2-X,'fro');
    res_d(2) = norm(V2-V2_old,'fro');
    %% update Lambda3
    Lambda3 = Lambda3 + gamma*(V3 - E*V2 + y);
    res_p(3) = norm(V3 - E*V2 + y,'fro');
    res_d(3) = norm(V3-V3_old,'fro');
    %% update Lambda4
    Lambda4 = Lambda4 + gamma*(V4 - X);
    res_p(4) = norm(V4-X,'fro');
    res_d(4) = norm(V4-V4_old,'fro');
    %% update Lambda5
    Lambda5x = Lambda5y + gamma*(V5x - Dh(X_reshape));
    Lambda5y = Lambda5y + gamma*(V5y - Dv(X_reshape));
    res_p(5) = norm(V5x-Dh(X_reshape) + V5y-Dv(X_reshape),'fro');
    res_d(5) = norm(V5x + V5y-V5x_old - V5y_old,'fro');
%%
    res_p_total = mean(res_p);
    res_d_total = mean(res_d);
    if mod(i,10) ==0 
        if res_p_total> 10*res_d_total
                gamma = step*gamma;
            elseif res_d_total > 10 *res_p_total 
                gamma = gamma./step;
        end
        % Update gamma to keep primal and dual residuals within a factor of 10
        % fprintf('Iter: %d, res_p = %f, res_d = %f\n', i, res_p_total, res_d_total);
    end
    i = i+1;
end
X = gather(X);
end

function V = Soft1(p_f, p_s, q, X, V0, E,Y,Lambda, gamma)
        V0 = V0+1e-10;
        ln_v = log(abs(V0));

        % 计算 numerator 和 denominator
        numerator = abs(V0).^p_f .* (p_f * ln_v - 1) - abs(V0).^p_s .* (p_s * ln_v - 1);
        denominator = (p_f - p_s) * q * abs(V0).^q .* (ln_v .^ 2);

        % 计算 phi
        Phi = numerator ./ denominator;
        % 计算 Y
        Yy = -Lambda./gamma - Y + E*X;
        % disp(max(max(Phi)))
        % disp(min(min(Phi)))
    if q==1
        % 进行软阈值操作
        V = soft_threshold(Yy,  Phi/gamma);
    elseif q==2
        V = Yy ./ (2  * Phi/gamma + 1);
    end
end
function output = soft_threshold(input, lambda)
    % 软阈值函数
    output = sign(input) .* max(abs(input) - lambda, 0);
end

