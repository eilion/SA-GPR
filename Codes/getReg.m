function [MU,NU] = getReg(A,X,Y,Z,ETA,GAMMA,RR,setting)

lambda0 = 1e-2;

X = X(Z==0);
Y = Y(Z==0);
RR = RR(Z==0);

N = size(X,1);

N0 = 200;
if N > N0
    % RR = max(RR,0.01);
    
    X0 = linspace(min(X),max(X),N0)';
    X0 = X0 + (max(X)-min(X))/N0*(rand()-1/2);
    % X0 = min(X) + (max(X)-min(X))*rand(N0,1);
    X0 = sort(X0,'ascend');
    if strcmp(setting.kernel,'SE') == 1
        K_00 = ETA^2*exp(-GAMMA^2*(X0-X0').^2) + lambda0^2*eye(N0);
        K_X0 = ETA^2*exp(-GAMMA^2*(X-X0').^2);
        K_A0 = ETA^2*exp(-GAMMA^2*(A-X0').^2);
    elseif strcmp(setting.kernel,'OU') == 1
        K_00 = ETA^2*exp(-GAMMA^2*abs(X0-X0')) + lambda0^2*eye(N0);
        K_X0 = ETA^2*exp(-GAMMA^2*abs(X-X0'));
        K_A0 = ETA^2*exp(-GAMMA^2*abs(A-X0'));
    elseif strcmp(setting.kernel,'SE[SE]') == 1
        K_00 = ETA^2*(1+GAMMA(2)^2*(1-exp(-GAMMA(1)^2*(X0-X0').^2))).^(-1/2) + lambda0^2*eye(N0);
        K_X0 = ETA^2*(1+GAMMA(2)^2*(1-exp(-GAMMA(1)^2*(X-X0').^2))).^(-1/2);
        K_A0 = ETA^2*(1+GAMMA(2)^2*(1-exp(-GAMMA(1)^2*(A-X0').^2))).^(-1/2);
    end
    WW = K_00 + (K_X0./sqrt(RR))'*(K_X0./sqrt(RR));
    
    MU = K_A0*(WW\(K_X0'*(Y./RR)));
    NU = ETA^2 + lambda0^2 + sum(K_A0.*(WW\K_A0'-K_00\K_A0')',2);
else
    if strcmp(setting.kernel,'SE') == 1
        K_AX = ETA^2*exp(-GAMMA^2*(A-X').^2);
        K_XX = ETA^2*exp(-GAMMA^2*(X-X').^2) + lambda0^2*eye(N);
    elseif strcmp(setting.kernel,'OU') == 1
        K_AX = ETA^2*exp(-GAMMA^2*abs(A-X'));
        K_XX = ETA^2*exp(-GAMMA^2*abs(X-X')) + lambda0^2*eye(N);
    elseif strcmp(setting.kernel,'SE[SE]') == 1
        K_AX = ETA^2*(1+GAMMA(2)^2*(1-exp(-GAMMA(1)^2*(A-X').^2))).^(-1/2);
        K_XX = ETA^2*(1+GAMMA(2)^2*(1-exp(-GAMMA(1)^2*(X-X').^2))).^(-1/2) + lambda0^2*eye(N);
    end
    K_XX = K_XX + diag(RR);
    
    MU = K_AX*(K_XX\Y);
    NU = ETA^2 + lambda0^2 - sum(K_AX.*(K_XX\K_AX')',2);
end


end