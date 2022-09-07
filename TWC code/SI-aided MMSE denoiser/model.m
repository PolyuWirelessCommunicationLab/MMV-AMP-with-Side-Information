clear;
%%%generate side information
MTC=10000;%%Monte Carlo
N=2000;%user
L1=250;%pilot length
K=200;%active user
lambda=K/N;
M=2; %number of antenna
Block=10;%number of blocks in each realization
D_abs_without_store = zeros(N,Block,MTC);
D_abs_withsi_store = zeros(N,Block,MTC);
D_act_store = zeros(N,Block,MTC);
D_tau_without_store=zeros(Block,MTC);
D_tau_withsi_store=zeros(Block,MTC);
D_tau_dcs_store = zeros(Block, MTC);
path_loss_store=zeros(N,MTC);
Diff_L_MSE=zeros(60,Block,10);
Total_MSE=zeros(60,Block,MTC);
parfor mtc = 1: MTC
        N=2000;%user
        L1=250;%pilot length
        K=200;%active user
        M=2;
        Block=10;%number of blocks in each realization
        alpha=0.55;%p11, inactive to inactive
        beta=0.05;%p10, inactive to active
        p01=1-alpha;%transition probability
        p10=beta;%transition probability
        maxN_itera =60; %iteration in each block
        D_dcs_abs = zeros(N,Block);
        D_dcs_signal = zeros(N,M,Block);
        D_tau_withoutsi = zeros(Block,1);
        D_tau_withsi = zeros(Block,1);
        D_tau_dcs = zeros(Block,1);
        lambda_in=zeros(N,1);
        lambda_out = zeros(N,1);
        pi_out = zeros(N,1);
        phi_out=zeros(N,1);
        pi_in =zeros(N,1);
%%%%generate channel information
        path_loss = zeros(N,1);
        display(strcat('MTC=',num2str(mtc)));
        D_act = false(N,Block);%initialize act matrix for each user
        D_channel = zeros(N,M,Block);%initialize channel matrix
        D_abs_without=zeros(N,Block);
        D_abs_withsi=zeros(N,Block);
        D_signal_without = zeros(N,M,Block);%initialize denoised signal matrix
        D_signal_withsi = zeros(N,M,Block);
        MSE_withoutsi=zeros(maxN_itera,Block);%initialize MSE matrix
        MSE_withsi=zeros(maxN_itera,Block);
        %%generate measurement matrix
        A = randn(L1,N)*(1/L1)^0.5*sqrt(1/2) + sqrt(-1)*randn(L1,N)*(1/L1)^0.5*sqrt(1/2);
        %%calculate pathloss
         x_user=(rand(N,1)-0.5)*500;
         y_user=(rand(N,1)-0.5)*500;  % location
        distance = zeros(N,1)*1000; 
        for n=1:N
            distance(n)=sqrt(x_user(n)^2+y_user(n)^2);
            path_loss(n)=10^((-128.1-37.6*log10(distance(n)*1e-3))/10);
        end
        path_loss_store(:,mtc) = path_loss;
        power = 10^(2.3)*10^(-3);
        noise_power = 10^(-16.9)*10^(-3);
        B = 1e7;
        noise = noise_power*B;
        h=zeros(N,M);
        for blk = 1:Block
        %%generate channel model
            for n=1:N
                h(n,:) = sqrt(1/2)*(randn(1,M) + sqrt(-1)*randn(1,M))*sqrt(path_loss(n));
            end
            D_channel(:,:,blk) = h;
    %%%generate signal 
    %%generate active users set
    %%in the first block, P(active)= lambda, in next blocks, the
    %%probability will change according to transition probability, eg.
    %%alpha and beta
            if blk ==1
                supp = randperm(N);
                D_act(supp(1:K),blk) = true;
                 x = h.*D_act(:,1);
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
                x = h.*D_act(:,blk);%signal 
            end
     %%generate noise
            w = randn(L1,M)*sqrt(1/2) + sqrt(-1)*randn(L1,M)*sqrt(1/2);
            noise_r = noise/power/L1;
            sigma_w = sqrt(noise_r);
            Z= w*sigma_w;
     %%%Received signal at BS
            y = A*x + Z;
     %%%amp and amp-si algorithms, in the first block, there is no side
     %%information
     tic;
            if blk == 1
                [xnoise_withsi,xhat_withsi,mse_withsi,tau_est_withsi] = noisyCAMPmmseforKLS(A,N,M,L1,y,x,maxN_itera,K/N,path_loss,sigma_w);
                D_tau_withsi(1) = tau_est_withsi(end);
                D_signal_withsi(:,:,blk) = xnoise_withsi;
                MSE_withsi(:,blk) = mse_withsi;
                sigma_1 = tau_est_withsi(end);
            else
                [xnoise_withsi,xhat_withsi,mse_withsi,tau_est_withsi] = noisyCAMPmmseforKLSSI(A,N,M,L1,y,x,maxN_itera,K/N,path_loss,p01,p10,sigma_1,xnoise_withsi,sigma_w);   
                D_tau_withsi(blk) = tau_est_withsi(end);           
                D_signal_withsi(:,:,blk) = xnoise_withsi;
                MSE_withsi(:,blk) = mse_withsi;
                sigma_1 =  tau_est_withsi(end);
            end
        end
        toc;
     for blk = 1:Block
        for n = 1:N
            D_abs_withsi(n,blk) = norm(D_signal_withsi(n,:,blk));            
        end
     end
       D_act_store(:,:,mtc) = D_act;
       D_tau_withsi_store(:,mtc) = D_tau_withsi;
       D_abs_withsi_store(:,:,mtc) = D_abs_withsi;
       Total_MSE(:,:,mtc)=MSE_withsi;
end   
%%%detection
for blk=1:2:10
    P_md_h=zeros(1,15);
    P_fa_h=zeros(1,15);
    d_K_sort = sort(D_abs_withsi_store(:,blk,:),'descend');
    for j = 1:15
    N_fa=0;
    N_md=0;
        for mtc = 1:MTC
            thre=(1.5-0.1*j)*d_K_sort(K+1,mtc); %%%%%note that the threshold should be changed case by case to generate the tradeoff curve   
            for n = 1:N
                if D_abs_withsi_store(n,blk,mtc) > thre && D_act_store(n,blk,mtc)==0
                    N_fa = N_fa +1;
                 end
                 if D_abs_withsi_store(n,blk,mtc) <= thre && D_act_store(n,blk,mtc)==1
                    N_md = N_md + 1;
                end
           end 
        end
        P_md(j)=N_md/K/MTC ;
        P_fa(j)=N_fa/(N-K)/MTC;
    end
    loglog(P_fa,P_md);
    hold on
end
