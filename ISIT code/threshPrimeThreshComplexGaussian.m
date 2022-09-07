function [eta,etaPrimeAvg] = threshPrimeThreshComplexGaussian(y,N,M,sigma,lambda,pLS)
eta = zeros(N,M);
etaPrimeAvg = zeros(M,M);
for n=1:N
    a = pLS(n)/(pLS(n)+sigma^2);
    b = (1-lambda)/lambda*((pLS(n)+sigma^2)/sigma^2)^M;
    c = pLS(n)/sigma^2/(pLS(n)+sigma^2);
    coeff0 = a^2/sigma^2;
    t0 = b*exp(-c*abs(y(n,:)*y(n,:)'));
    t = 1 + t0;
    coeff1 = a/t;
    eta(n,:) = coeff1*y(n,:);    
    etaPrimeMtx = coeff1*eye(M) + coeff0*(y(n,:)'*y(n,:))*t0/t^2;
    etaPrimeAvg = etaPrimeAvg + etaPrimeMtx;
end
etaPrimeAvg = etaPrimeAvg/N;
end
