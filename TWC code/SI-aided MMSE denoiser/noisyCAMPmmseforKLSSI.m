function [xnoise_withsi,x_withsi,mse,tau_est] =noisyCAMPmmseforKLSSI(A,N,M,L,y,xsig,maxN_itera,lambda,path_loss,p01,p10,sigma_1,xnoise_withsi,sigma_w)
% Complex Approximate message passing Fixed LS
z = y;
x_withsi = zeros(N,M);
sum_gain = 0;
for m=1:M
    sum_gain = sum_gain+(norm(y(:,m)))^2;
end
tau=0;
tau = sqrt(sum_gain)*sqrt(1/(M*L));  
mse = zeros(maxN_itera,1);
tau_real = zeros(maxN_itera,1);
tau_est = zeros(maxN_itera,1);
tau_real(1) = tau;
tau_est(1) = tau;
for i = 1:maxN_itera 
%     display(num2str(i));
    input = A'*z + x_withsi;
    [x_withsi,avgxprime] = threshPrimeThreshComplexGaussianwithsi(input,N,M,sigma_1,tau,lambda,path_loss,xnoise_withsi,p01,p10);
     if i==1
         xold=x_withsi;
     else
         x_withsi=0.95*x_withsi+0.05*xold;
     end
    xold=x_withsi;
    for n=1:N
        mse(i) = mse(i) + (norm(x_withsi(n,:) - xsig(n,:)))^2;
    end
    mse(i) = mse(i)/norm(xsig)^2;
    tau_real(i+1) = sqrt(sigma_w^2 + N/L*mse(i));
%     avgxprime = zeros(M,M);
%     for n=1:N
%         avgxprime = avgxprime + xprime(:,:,n);
%     end
%     avgxprime = avgxprime/N;
    z = y - A*x_withsi + N/L*z*avgxprime; 
    sum_gain = 0;
    for m = 1:M
        sum_gain = sum_gain + (norm(z(:,m)))^2;
    end
    tau=0;
    tau = sqrt(sum_gain)*sqrt(1/(M*L)); 
    if i ~= maxN_itera
        tau_est(i+1) = tau;
    end
end
aaaa=avgxprime;
xnoise_withsi = A'*z + x_withsi;
%plot(1:50,10*log10(abs(mse)));
end


