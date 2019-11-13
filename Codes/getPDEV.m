function [PDEV] = getPDEV(Y,X,R,param)

PDEV = zeros(3,1);

N = size(Y,1);

W = getCov(X,X,param.eta,param.xi) + param.lambda^2*eye(N) + param.lambda0^2*eye(N) + diag(R);
W_inv = W\eye(N);
alpha = W_inv*Y;
ALPHA = alpha*alpha' - W_inv;

% eta:
WW = 2*param.eta*exp(-param.xi^2*(X-X').^2);
PDEV(1) = 0.5*trace(ALPHA*WW);

% xi:
WW = param.eta^2*exp(-param.xi^2*(X-X').^2).*(-2*param.xi*(X-X').^2);
PDEV(2) = 0.5*trace(ALPHA*WW);

% lambda:
WW = 2*param.lambda*eye(N);
PDEV(3) = 0.5*trace(ALPHA*WW);


end