function [K] = getCov(X1,X2,eta,xi)

K = eta^2*exp(-xi^2*(X1-X2').^2);


end