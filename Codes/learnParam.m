function [data] = learnParam(data,Samples,Profile,setting)

beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
gamma0 = 1e-3;

max_iters = 30000;

L = length(data);

if strcmp(setting.isLearn_R,'local') == 1
    for ll = 1:L
        Z = Samples.Z(Samples.ID==ll,:);
        X = data(ll).X;
        data(ll).R = mean((Z(end,:)-Z(1,:))./(X(end)-X(1)));
    end
elseif strcmp(setting.isLearn_R,'global') == 1
    AA = 0;
    BB = 0;
    for ll = 1:L
        Z = Samples.Z(Samples.ID==ll,:);
        X = data(ll).X;
        
        T = (Z(2:end,:)-Z(1:end-1,:))./(X(2:end)-X(1:end-1));
        
        AA = AA + sum(sum(T))/size(T,2);
        BB = BB + size(T,1) - 1;
    end
    for ll = 1:L
        data(ll).R = AA/BB;
    end
end

if strcmp(setting.isLearn_gamma,'local') == 1
    for ll = 1:L
        ALPHA = data(ll).alpha;
        BETA = data(ll).beta;
        
        Mw = zeros(3,1);
        Vw = zeros(3,1);
        
        Z = Samples.Z(Samples.ID==ll,:);
        X = data(ll).X;
        
        [N,K] = size(Z);
        
        T = (Z(2:end,:)-Z(1:end-1,:))./(X(2:end)-X(1:end-1));
        
        AA = sum(sum(log(T)))/K;
        BB = sum(sum(T))/K;
        
        for r = 1:max_iters
            PDEV = zeros(3,1);
            
            PDEV(1) = AA + (setting.s+N-1)*log(BETA) - (setting.r+N-1)*psi(ALPHA) + log(setting.p) - (N-1)*log(data(ll).R);
            PDEV(2) = ALPHA/BETA*(setting.s+N-1) - (setting.q+BB/data(ll).R);
            % PDEV(3) = BETA*R^(-2)*BB - (ALPHA-1)*(N-1)/R;
            
            Mw = beta1*Mw - (1-beta1)*PDEV;
            Vw = beta2*Vw + (1-beta2)*PDEV.*PDEV;
            
            ALPHA = abs(ALPHA - (gamma0*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon));
            BETA = abs(BETA - (gamma0*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon));
            % R = abs(R - (gamma0*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(3)./(sqrt(Vw(3))+epsilon));
            
            ALPHA = min(ALPHA,10);
            % ALPHA = max(1,min(ALPHA,10));
            % BETA = min(BETA,10);
        end
        
        data(ll).alpha = ALPHA;
        data(ll).beta = BETA;
    end
elseif strcmp(setting.isLearn_gamma,'global') == 1
    ALPHA = data(1).alpha;
    BETA = data(1).beta;
    
    AA = cell(L,1);
    BB = cell(L,1);
    NN = zeros(L,1);
    
    for ll = 1:L
        Z = Samples.Z(Samples.ID==ll,:);
        X = data(ll).X;
        
        [N,K] = size(Z);
        
        T = (Z(2:end,:)-Z(1:end-1,:))./(X(2:end)-X(1:end-1));
        
        AA{ll} = sum(sum(log(T)))/K;
        BB{ll} = sum(sum(T))/K;
        
        NN(ll) = N - 1;
    end
    
    Mw = zeros(3,1);
    Vw = zeros(3,1);
    
    for r = 1:max_iters
        PDEV = zeros(3,1);
        
        for ll = 1:L
            PDEV(1) = AA{ll} + (setting.s+NN(ll))*log(BETA) - (setting.r+NN(ll))*psi(ALPHA) + log(setting.p) - NN(ll)*log(data(ll).R);
            PDEV(2) = ALPHA/BETA*(setting.s+NN(ll)) - (setting.q+BB{ll}/data(ll).R);
        end
        
        Mw = beta1*Mw - (1-beta1)*PDEV;
        Vw = beta2*Vw + (1-beta2)*PDEV.*PDEV;
        
        ALPHA = abs(ALPHA - (gamma0*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(1)./(sqrt(Vw(1))+epsilon));
        BETA = abs(BETA - (gamma0*sqrt(1-beta2^r)/(1-beta1^r)).*Mw(2)./(sqrt(Vw(2))+epsilon));
        
        ALPHA = max(1,min(ALPHA,10));
        % BETA = min(BETA,10);
    end
    
    for ll = 1:L
        data(ll).alpha = ALPHA;
        data(ll).beta = BETA;
    end
end

if strcmp(setting.isLearn_shift,'true') == 1
    for ll = 1:L
        Z = Samples.Z(Samples.ID==ll,:);
        Y = data(ll).Y;
        
        [~,K] = size(Z);
        
        AA = 0;
        BB = 0;
        for k = 1:K
            MU = interp1(Profile(:,1),Profile(:,2),Z(:,k));
            NU = interp1(Profile(:,1),Profile(:,3).^2,Z(:,k));
            
            AA = AA + sum(1./NU);
            BB = BB + sum((Y-MU)./NU);
        end
        
        data(ll).h = BB/AA;
    end
end


end