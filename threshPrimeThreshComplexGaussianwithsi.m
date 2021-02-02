function [eta,etaPrimeAvg] = threshPrimeThreshComplexGaussianwithsi(y,N,M,sigma_1,tau,lambda,pLS,xnoise_withsi,p01,p10)
eta = zeros(N,M);
etaPrimeAvg = zeros(M,M);
for n=1:N
        a = pLS(n)/(pLS(n)+tau^2);
        b = (((pLS(n)+sigma_1^2)/sigma_1^2)^M)*exp(-(1/sigma_1^2-1/(pLS(n)+sigma_1^2))*abs(xnoise_withsi(n,:)*xnoise_withsi(n,:)'));
        %c = (1/(sigma_1^2))*exp(-abs(xnoise_withsi(n,:)*xnoise_withsi(n,:)')/sigma_1^2);
        %(1-p01)*lambda*b+p10*(1-lambda)*c;
        d = (p10+(1-p10)*b)/((1-p01)+p01*b);
        e = ((pLS(n)+tau^2)/tau^2)^M;
        f = 1/tau^2-1/(pLS(n)+tau^2);
        coeff0 = a^2/tau^2;
        t0 = d*e*exp(-f*abs(y(n,:)*y(n,:)'));
        t = 1 + (1-lambda)/lambda*t0;
        coeff1 = a/t;
        eta(n,:) = coeff1*y(n,:);
        etaPrimeMtx = coeff1*eye(M) + (1-lambda)/lambda*coeff0*(y(n,:)'*y(n,:))*t0/t^2;
        etaPrimeAvg = etaPrimeAvg + etaPrimeMtx;
%     a = pLS(n)/(pLS(n)+tau^2);
%     b = 1/(pLS(n)+sigma_1^2)^M*exp(-abs(xnoise_withsi(n,:)*xnoise_withsi(n,:)')/(pLS(n)+sigma_1^2));
%     c = 1/(sigma_1^2)^M*exp(-abs(xnoise_withsi(n,:)*xnoise_withsi(n,:)')/sigma_1^2);
% %     (1-p01)*lambda*b+p10*(1-lambda)*c;
%     d = (p01*lambda*b+(1-p10)*(1-lambda)*c)/((1-p01)*lambda*b+p10*(1-lambda)*c);
%     e = ((pLS(n)+tau^2)/tau^2)^M;
%     f = 1/tau^2-1/(pLS(n)+tau^2);
%     coeff0 = a^2/tau^2;
%     t0 = d*e*exp(-f*abs(y(n,:)*y(n,:)'));
%     t = 1 + t0;
%     coeff1 = a/t;
%     eta(n,:) = coeff1*y(n,:);
%     etaPrimeMtx = coeff1 + coeff0*(y(n,:)'*y(n,:))*t0/t^2;
%     etaPrimeAvg = etaPrimeAvg + etaPrimeMtx;
end
etaPrimeAvg = etaPrimeAvg/N;
end