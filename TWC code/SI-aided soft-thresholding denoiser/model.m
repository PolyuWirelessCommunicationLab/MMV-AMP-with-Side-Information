clear;
tic;
N=2000; %total number of users
MTC=10000; %number of monte-carlo 
K=200; %number of active users at each coherence block
L=500; %length of pilot sequence
alpha=0.55; %p11, inactive to inactive
beta=0.05; %p10, inactive to active
p01=1-alpha; %transition probability
p10=beta; %transition probability
Block=10; %number of consecutive coherence blocks
D_abs_store=zeros(N,Block,MTC);
D_act_store=zeros(N,Block,MTC);
maxN_itera=50; %%set maxN_itera as 40 with two antennas
flag=[];
%MSE=zeros(maxN_itera,Block);
A = randn(L,N)*(1/L)^0.5*sqrt(1/2) + sqrt(-1)*randn(L,N)*(1/L)^0.5*sqrt(1/2); %pilot sequence
parfor mtc=1:MTC
    K=200;
    M=1; %number of antennas at the BS
    lambda=K/N; 
    L=500;
    display(strcat('MTC=',num2str(mtc)));
    D_act=zeros(N,Block);        
    D_abs=zeros(N,Block);
    D_act_est=zeros(N,1);
    tau_ministore=zeros(maxN_itera,Block);
    distance=zeros(N,1);
    path_loss=zeros(N,1);
    h=zeros(N,M);
    mse=zeros(50,1);
    for blk =1:Block
        N_fa=0;
        N_md=0;
        D_channel=zeros(N,M);
    %%%generate channel
        %generate users locations
        x_user=(rand(N,1)-0.5)*500;
        y_user=(rand(N,1)-0.5)*500;
        for n=1:N
           distance(n) = sqrt(x_user(n)^2+y_user(n)^2);
           path_loss(n) = 10^((-128.1 - 36.7*log10(distance(n)*1e-3))/10);
           h(n,:) = sqrt(1/2)*(randn(1,M) + sqrt(-1)*randn(1,M))*sqrt(path_loss(n));
           D_channel(n,:)=h(n,:);
        end
        %%%generate correlated activity pattern
        if blk ==1
           supp = randperm(N);
           D_act(supp(1:K),blk) = true;
           x = D_channel.*D_act(:,1);
        else
           D_act(:,blk) = D_act(:,blk-1);
           Inacts = find(D_act(:,blk-1) == 0);
           Inacts_Acts = Inacts(randperm(length(Inacts)));
           new_Acts = Inacts_Acts(1:round(length(Inacts)*beta));
           D_act(new_Acts,blk) = true;
           Acts = find(D_act(:,blk-1) == 1);
           Acts_Inacts = Acts(randperm(length(Acts)));
           new_Inacts = Acts_Inacts(1:round(length(Acts)*(1-alpha)));
           D_act(new_Inacts,blk) = false;
           supp= find(D_act(:,blk) == 1);
           K = numel(supp);
           x = D_channel.*D_act(:,blk);%signal 
        end
        D_act_store(:,:,mtc)=D_act;
        power = 10^(2.3)*10^(-3); %transmit power
        noise_power = 10^(-16.9)*10^(-3); %noise power
        B = 1e7; %bandwidth
        %%%generate noise
        noise = noise_power*B;
        w = randn(L,M)*sqrt(1/2) + sqrt(-1)*randn(L,M)*sqrt(1/2);
        noise_r = noise/power/L;
        sigma_w = sqrt(noise_r);
        z=w*sigma_w;
        %%%received signal at the BS
        y=A*x+z;
        %%% start the SI-aided AMP algorithm to recover x based on y
        
        %%%First, calculate the weight \theta for SI-aided soft threholding
        %%%denoiser
        if blk==1 %in the first coherence block, no SI is available
            thre1=calweight(1-lambda,lambda,M);
            thre2=calweight(1-lambda,lambda,M);
        else
            %in the remaining coherence blocks, SI is available. We
            %calculate weight based on Theorem 2.
            new_alpha1=(1-beta)*(1-lambda)*p_fa_est+(1-alpha)*lambda;
            new_alpha2=(beta*(1-lambda)*p_fa_est+alpha*(lambda));
            beta1 =(alpha*lambda+beta*(1-lambda)*(1-p_fa_est))/(lambda+(1-lambda)*(1-p_fa_est));
            thre1 = calweight(new_alpha1,new_alpha2,M);%active at the previous coherence block
            thre2 = calweight(1-beta,beta,M);%inactive at the previous coherence block
        end
        %%SI
        if blk ==1
            D_act_est_slot=D_act(:,1);
        else
            D_act_est_slot=D_act_est;
        end
        %%%SI-aided AMP, eq. (6) and (7) in the paper
         [x_noise,x_est,mse,tau_est] = mmvsoftthreshoding(A,N,y,x,M,L,maxN_itera,thre1,thre2,D_act_est_slot);
         tau_store(:,blk,mtc)=tau_est;
         for n=1:N
             D_abs(n,blk) = norm(x_noise(n,:));
         end
         %%%%estimate the activity in this block, which will be used to
         %%%%assist next block's detection;
         D_abs_sort = sort(D_abs(:,blk),'descend');
         threshold_est = D_abs_sort(K+1);%set a detection threshold
         for n = 1:N
             if D_abs(n,blk)>threshold_est
                 D_act_est(n) = true;
             else
                 D_act_est(n) = false;
             end
         end
         p_fa_est=igamma(M,(threshold_est/min(tau_est))^2)/gamma(M);%%calculate (27)
         D_abs_store(:,:,mtc)=D_abs;
    end
end
%%%%%%%start detection
for blk=1:2:10
    P_md_h=zeros(1,15);
    P_fa_h=zeros(1,15);
    d_K_sort = sort(D_abs_store(:,blk,:),'descend');
for j =1:15
    N_fa=0;
    N_md=0;
    for mtc = 1:MTC             
     thre=(1.5-0.1*j)*d_K_sort(K+1,mtc);   %%%%%note that the threshold should be changed case by case to generate the tradeoff curve   
        for n =1:N
              if D_abs_store(n,blk,mtc)>thre&&D_act_store(n,blk,mtc)==0
                   N_fa=N_fa+1;
              elseif D_abs_store(n,blk,mtc)<=thre&&D_act_store(n,blk,mtc)==1
                   N_md=N_md+1;
              end
        end
    end
    P_fa(j)=N_fa/MTC/(N-K);
    P_md(j)=N_md/MTC/K;
end
loglog(P_fa,P_md);
    hold on;
end




toc;
