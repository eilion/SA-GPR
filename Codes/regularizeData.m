function [data] = regularizeData(data)

% Regularize inputs:
X = data.X0;
[N,K] = size(X);

mu = sum(X,1)/N;
SIG = zeros(K);
for k = 1:K
    SIG(k,k) = var(X(:,k));
end

nu = 6;
DIFF = 1;
r = 0;
while r < 100 && DIFF > 1e-4
    r = r + 1;
    
    DD = sum((X-mu).*(SIG\((X-mu)'))',2);
    ZZ = (nu+K)./(nu+DD);
    
    mu_new = sum(ZZ.*X,1)/sum(ZZ);
    SIG = ((X'-mu')*(ZZ.*(X-mu)))/N;
    
    DIFF = sqrt(sum((mu_new-mu).^2)/K);
    mu = mu_new;
end

QQ = SIG\eye(K);
QQ = QQ^(1/2);
data.X = (X-mu)*QQ;
data.mu_X = mu;
data.sig_X = QQ\eye(K);


% Regularize outputs:
Y = data.Y0;
[N,K] = size(Y);

mu = sum(Y,1)/N;
SIG = zeros(K);
for k = 1:K
    SIG(k,k) = var(Y(:,k));
end

nu = 6;
DIFF = 1;
r = 0;
while r < 100 && DIFF > 1e-4
    r = r + 1;
    
    DD = sum((Y-mu).*(SIG\((Y-mu)'))',2);
    ZZ = (nu+K)./(nu+DD);
    
    mu_new = sum(ZZ.*Y,1)/sum(ZZ);
    SIG = ((Y'-mu')*(ZZ.*(Y-mu)))/N;
    
    DIFF = sqrt(sum((mu_new-mu).^2)/K);
    mu = mu_new;
end

QQ = SIG\eye(K);
QQ = QQ^(1/2);
data.Y = (Y-mu)*QQ;
data.mu_Y = mu;
data.sig_Y = QQ\eye(K);


data.F = [ones(N,1),data.X];


end