function [x_noise,x_est,mse,tau_est] = mmvsoftthreshoding(A,N,y,x,M,L1,maxN_itera,thre1,thre2,D_act_est_slot)
%%initialization 
r=y;
sum_gain=0;
for m = 1:M
    sum_gain = sum_gain+(norm(y(:,m)))^2;
end
tau = sqrt(sum_gain)*sqrt(1/(M*L1));
mse=zeros(maxN_itera,1);
x_est=zeros(N,M);
tau_est=zeros(maxN_itera,1);
%%start iteration
for i = 1:maxN_itera
    input = A'*r+x_est;
    tau_est(i)=tau;  
    tmp=0;
    etaprime=zeros(M,M);
    avgxprime=zeros(M,M);
    for n =1:N
        if D_act_est_slot(n)==1
            thre = thre1*tau;
        end
        if D_act_est_slot(n)==0
            thre = thre2*tau;
        end
        tmp = (norm(input(n,:))-thre);
        %%%eq. (15) in the paper
        if tmp>=0
            x_est(n,:) = input(n,:)-thre*input(n,:)/norm(input(n,:));
        else
            x_est(n,:) = zeros(1,M);
        end
        etaprime=(1*eye(M)-thre/norm(input(n,:))*eye(M)+0.5*thre*input(n,:)'*input(n,:)/(norm(input(n,:)))^3)*(tmp>=0);
        avgxprime=avgxprime+etaprime;  
    end
    %%%damping
    if i==1
        xold=x_est;
    else
        x_est=0.9*x_est+0.1*xold;
    end
    xold=x_est;
    avgxprime=avgxprime/(N);
    for n=1:N
        mse(i) = mse(i) + (norm(x_est(n,:) - x(n,:)))^2;
    end
    mse(i) = mse(i)/(N*M);
    r = y - A*x_est + N/L1*r*avgxprime;
    sum_gain = 0;
    %%Estimation of tau
    for m = 1:M
        sum_gain = sum_gain + (norm(r(:,m)))^2;
    end
    tau = sqrt(sum_gain)*sqrt(1/(M*L1));  
end
x_noise = A'*r + x_est;
end
