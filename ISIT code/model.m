clear;
%%%generate side information
MTC=1e5;%%Monte Carlo
N=4000;%user
L1=450;%pilot length
K=400;%active user
M=2;
alpha=0.91;%p11, inactive to inactive
beta=0.01;%p10, inactive to active
Block=10;%number of blocks in each realization
D_abs_without_store = zeros(N,Block,MTC);
D_abs_withsi_store = zeros(N,Block,MTC);
D_act_store = zeros(N,Block,MTC);
D_tau_without_store=zeros(Block,MTC);
D_tau_withsi_store=zeros(Block,MTC);
path_loss_store=zeros(N,MTC);
%%%execute parallel computing
parfor mtc =  1: MTC
        display(strcat('MTC=',num2str(mtc)));
        K=400;%active user
        p01=1-alpha;%transition probability
        p10=beta;%transition probability
        maxN_itera =50; %iteration in each block
        D_tau_withoutsi = zeros(Block,1);
        D_tau_withsi = zeros(Block,1);     
        %initialize some matrix
        path_loss = zeros(N,1);
        D_act = false(N,Block);%initialize act matrix for each user
        D_channel = zeros(N,M,Block);%initialize channel matrix
        D_abs_without=zeros(N,Block);
        D_abs_withsi=zeros(N,Block);
        D_signal_without = zeros(N,M,Block);%initialize denoised signal matrix
        D_signal_withsi = zeros(N,M,Block);        
        %%generate measurement matrix
        A = randn(L1,N)*(1/L1)^0.5*sqrt(1/2) + sqrt(-1)*randn(L1,N)*(1/L1)^0.5*sqrt(1/2);
        %%generate pathloss for each user
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
            for n=1:N
                h(n,:) = sqrt(1/2)*(randn(1,M) + sqrt(-1)*randn(1,M))*sqrt(path_loss(n));
            end
            D_channel(:,:,blk) = h;
        %%%generate signal 
        %%generate active users set
        %%in the first block when there is no si, P(active)= lambda, then in next blocks, the
        %%probability will change according to transition probability, i.e.alpha and beta
            if blk ==1
                supp = randperm(N);
                D_act(supp(1:K),blk) = true;
                 x = zeros(N,M);
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
                x=zeros(N,M);
                x = h.*D_act(:,blk);%signal 
            end
     %%generate noise
            w = randn(L1,M)*sqrt(1/2) + sqrt(-1)*randn(L1,M)*sqrt(1/2);
            noise_r = noise/power/L1;
            sigma_w = sqrt(noise_r);
            Z= w*sigma_w;
     %%%Received signal at BS
            y = A*x + Z;
     %%%amp and amp-si algorithms, in the first block, there is no SI
            if blk == 1
                [xnoise_withoutsi,xhat_withoutsi,mse_withoutsi,tau_est_withoutsi] = noisyCAMPmmseforKLS(A,N,M,L1,y,x,maxN_itera,K/N,path_loss,sigma_w);
                [xnoise_withsi,xhat_withsi,mse_withsi,tau_est_withsi] = noisyCAMPmmseforKLS(A,N,M,L1,y,x,maxN_itera,K/N,path_loss,sigma_w);
                D_tau_withoutsi(blk) = min(tau_est_withoutsi);    
                D_tau_withsi(1) = tau_est_withsi(end);
                D_signal_without(:,:,blk) = xnoise_withoutsi;
                D_signal_withsi(:,:,blk) = xnoise_withsi;
                sigma_1 = tau_est_withsi(end);
    %%%%%%user activity detection
            else
                [xnoise_withoutsi,xhat_withoutsi,mse_withoutsi,tau_est_withoutsi] = noisyCAMPmmseforKLS(A,N,M,L1,y,x,maxN_itera,K/N,path_loss,sigma_w);
                [xnoise_withsi,xhat_withsi,mse_withsi,tau_est_withsi] = noisyCAMPmmseforKLSwithsi(A,N,M,L1,y,x,maxN_itera,K/N,path_loss,p01,p10,sigma_1,xnoise_withsi,sigma_w);
                D_tau_withoutsi(blk) = min(tau_est_withoutsi);   
                D_tau_withsi(blk) = tau_est_withsi(end);
                D_signal_without(:,:,blk) = xnoise_withoutsi;              
                D_signal_withsi(:,:,blk) = xnoise_withsi;
                sigma_1 = tau_est_withsi(end);
            end
        end
    %%%%%calculate abs for each user to conduct user activity detection
     for blk = 1:Block
        for n = 1:N
            D_abs_without(n,blk) = norm(D_signal_without(n,:,blk));
            D_abs_withsi(n,blk) = norm(D_signal_withsi(n,:,blk));
        end
     end
       D_act_store(:,:,mtc) = D_act;
       D_tau_without_store(:,mtc) = D_tau_withoutsi;
       D_tau_withsi_store(:,mtc) = D_tau_withsi;
       D_abs_without_store(:,:,mtc) = D_abs_without;
       D_abs_withsi_store(:,:,mtc) = D_abs_withsi;       
end   
 %%%%%User activity detection
 %%%%%the case withoutsi
for blk = 1:1
    P_md_without=zeros(Block,1);
    P_fa_without=zeros(Block,1);
    P_md_h=zeros(1,1);
    P_fa_h=zeros(1,1);
    N_fa_without = zeros(MTC,1);
    N_md_without = zeros(MTC,1);
        for mtc = 1:MTC
            %%%%tune on the threshold to get the curve, this needs to be
            %%%%dealt with case by case
                Th_h=sort(D_abs_withsi_store(:,blk,mtc),'Descend');
                    threshold_h =0.78*Th_h(400);
                    for n = 1:N
                       if D_abs_withsi_store(n,blk,mtc) > threshold_h && D_act_store(n,blk,mtc)==0
                            N_fa_without(mtc) = N_fa_without(mtc) +1;
                        end
                        if D_abs_withsi_store(n,blk,mtc) < threshold_h && D_act_store(n,blk,mtc)==1
                            N_md_without(mtc) = N_md_without(mtc) + 1;
                        end
                    end 
        end
                    P_md_h=N_md_without/K;  
                    P_fa_h=N_fa_without/(N-K);
                    P_md_without = sum(P_md_h)/MTC;
                    P_fa_without = sum(P_fa_h)/MTC;
end
    %%%The case with si
    %     %%%calculate LLR for block blk
P_md_withsi=zeros(4,1);
P_fa_withsi=zeros(4,1);
for blk = 3:2:9
    for mtc = 1:MTC
                for n = 1:N
                    a2 =M*log((D_tau_withsi_store(blk,mtc)^2/(path_loss_store(n,mtc)+D_tau_withsi_store(blk,mtc)^2)))+(1/D_tau_withsi_store(blk,mtc)^2-1/(path_loss_store(n,mtc)+D_tau_withsi_store(blk,mtc)^2))*D_abs_withsi_store(n,blk,mtc)^2;                  
                    b2 = ((D_tau_withsi_store(blk-1,mtc)^2/(path_loss_store(n,mtc)+D_tau_withsi_store(blk-1,mtc)^2))^M)*exp((1/D_tau_withsi_store(blk-1,mtc)^2-1/(path_loss_store(n,mtc)+D_tau_withsi_store(blk-1,mtc)^2))*D_abs_withsi_store(n,blk-1,mtc)^2);
                    if isinf(b2) == 1
                        c2 = log(91);
                    else
                        c2 = log((0.009+0.091*b2)/(0.891+0.009*b2)*9);
                    end
                    LLR(n,mtc) = a2+c2;
                end
    end
 %%set threshold
%display(strcat('MTC=',num2str(MTC)));     
            N_fa_withsi = zeros(MTC,1);
            N_md_withsi = zeros(MTC,1);
                for mtc=1:MTC
                    Th_op =sort(LLR(:,mtc),'Descend');
                    threshold_op = 0;
                    for n = 1:N
                               if LLR(n,mtc) > threshold_op && D_act_store(n,blk,mtc)==0
                                    N_fa_withsi(mtc) = N_fa_withsi(mtc) +1;
                                end
                                if LLR(n,mtc) < threshold_op && D_act_store(n,blk,mtc)==1
                                    N_md_withsi(mtc) = N_md_withsi(mtc) + 1;
                                end
                    end
                end
                P_md_op=N_md_withsi/K;   
                P_fa_op=N_fa_withsi/(N-K);
                P_md_withsi((blk-1)/2) = sum(P_md_op)/MTC;
                P_fa_withsi((blk-1)/2 ) = sum(P_fa_op)/MTC;
end


      

