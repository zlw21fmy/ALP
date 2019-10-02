

%% 得出文中的图15  不同信噪比下含噪信号直接提取频移和滤波后信号提取的关系图
clc;
clear all;
close all;
%% 生成20个谱  
fn = 201;         %
delatvb = 0.04;   %
vb_single =linspace(10.69,10.71,201)';
temp=linspace(10.71,10.69,201)';
vb_single =[vb_single;temp(2:end)];
vb = repmat(vb_single,1,1);   % vb 

Fres = linspace(10.66,10.74,fn); 
% N=round(delatvb/2/(Fres(2)-Fres(1)));%半峰一半对应点数
N=40;
Fres = Fres';
num_1  = size(vb,1);     % 空间采样点数目  1200

num_2 = 20;                      % 重复次数
snr_given_matrix =5:1:20;  %给定信噪比
for j= 1:num_1
Gain_clean(:,j) = 1./(1+4*((Fres- vb(j))./delatvb).^2);   % 1*fn
end              %生成一个干净的谱   fn行  num_1列（100行  1200列）

Gain_clean_series = repmat(Gain_clean,1,num_2);  %20个干净的谱 横向串联
%% 
M = 10;
mu  = 0.008;
num_repeat = 1;
%% 滤波后信号迭代

for hh = 1:length(snr_given_matrix)
    snr_given = snr_given_matrix(hh);
    for  tt = 1:num_repeat
        Gain_noisy_series =  zeros(fn,num_1*num_2);
    
    
    for i = 1:fn
        Gain_noisy_series(i,:)= awgn(Gain_clean_series(i,:),snr_given,'measured');
    end               %串联后的谱按行加噪
    %%   定义lms滤波器参数
    period = 1;                                 % 参考输入信号 延迟是1
    x_matrix = [zeros(fn,period) Gain_noisy_series];  % 延迟体现为补零
    d_matrix = [ Gain_noisy_series zeros(fn,period)];
    for i = 1:fn
        x  = x_matrix(i,:) ;     % Input to the filter
        d  = d_matrix(i,:);      % Desired signal  延迟了一步
        ha = adaptfilt.lms(M,mu);
        [y,e] = filter(ha,x,d);
        y_out(i,:) =  y;         %滤波输出矩阵
        e_out(i,:)=  e;
    end
    
    num_4 = num_2;
    y_out_temp  = y_out(:,2:end);   %
    y_out_object  = y_out_temp (:,(1-1)*num_1+1:end);  % 从第三个谱开始
      y_out_object(:,1)=y_out(:,3);
    for p = 1:(num_4)
        
        y_test = y_out_object(:,(p-1)*num_1+1:(p)*num_1);  %选出滤波后的第p个谱
        [m1 ,locs1] = max(y_test, [], 1);  % 找出各列的最大值
        
        for jj = 1:num_1
            locs1_object = locs1(jj);
            if  locs1_object> N && locs1_object<=fn- N
                indice1 = locs1_object-(N):locs1_object+(N); %从最大值从左到右各找半个线宽
            elseif locs1_object<= N
                indice1 = 1:locs1_object+(locs1_object-1); %从最大值从左到右各找半个线宽
            elseif locs1_object> fn-N
                indice1 = locs1_object-(fn-locs1_object) :fn; %从最大值从左到右各找半个线宽
            end
            
            Fres1_object = Fres(indice1);
            p1 = polyfit(Fres1_object, y_test(indice1,jj),2);
            BFS_filtered_enssemble(jj,p) = -p1(2)./(2*p1(1));%利用抛物线的对称轴计算中心频移
        end
        
        for ii = 1:fn
            SNR_filtered_enssemble(ii,p) = snr(y_test(ii,:),abs(y_test(ii,:)-Gain_clean(ii,:)));
            
        end
    end
    
    BFS_filtered_error(tt,:) = mean(abs( BFS_filtered_enssemble(:,3:end)-repmat(vb,1,(num_4-2)))'); %从第三个到最后一个谱每个位置的误差的多次平均值
    SNR_filtered_mean(tt,:)  = mean(SNR_filtered_enssemble');
    disp([ '信噪比为' num2str(snr_given_matrix(hh)) '滤波后'] )
    end
BFS_filtered_error_mean(hh,:) = mean(BFS_filtered_error);    
SNR_filtered_mean_mean(hh,:) = mean(SNR_filtered_mean);

end

%%  

    %%  

  for hh = 1:length(snr_given_matrix)
    snr_given = snr_given_matrix(hh);
    for  tt = 1:num_repeat
        Gain_noisy_series =  zeros(fn,num_1*num_2);
    
    
    for i = 1:fn
        Gain_noisy_series(i,:)= awgn(Gain_clean_series(i,:),snr_given,'measured');
    end     
    num_4 = num_2;
    for p = 1:num_4
        gain_noisy_test = Gain_noisy_series(:,(p-1)*num_1+1:(p)*num_1);  %% 叠加完选出的第p个矩阵
        [m1 ,locs1] = max(gain_noisy_test, [], 1);  % 找出各列的最大值  即给定一个位置处的洛伦兹曲线的峰值
        
        for jj = 1:num_1
            locs1_object = locs1(jj);   % 找出最大值 对应的位置
            
            if  locs1_object> N && locs1_object<=fn- N
                indice1 = locs1_object-(N):locs1_object+(N); %从最大值从左到右各找半个线宽
            elseif locs1_object<= N
                indice1 = 1:locs1_object+(locs1_object-1); %从最大值从左到右各找半个线宽
                
            elseif locs1_object> fn-N
                indice1 = locs1_object-(fn-locs1_object) :fn; %从最大值从左到右各找半个线宽
                
            end
            
            Fres1_object = Fres(indice1);         %找到拟合范围对应的频率点值
            p1 = polyfit(Fres1_object,gain_noisy_test(indice1,jj),2);
            BFS(jj,p) = -p1(2)./(2*p1(1));%利用抛物线的对称轴计算中心频移
        end
        for ii = 1:fn
            SNR(ii,p) = snr(gain_noisy_test(ii,:),abs(gain_noisy_test(ii,:)-Gain_clean(ii,:)));
        end
    end
    BFS_noisy_error(tt,:) =  mean(abs(BFS-repmat(vb,1,num_4))');
    SNR_noisy (tt,:) = mean(SNR');
    disp([ '信噪比为' num2str(snr_given_matrix(hh))  '含噪'])
end
BFS_noisy_error_mean(hh,:) = mean(BFS_noisy_error);    
SNR_noisy_mean_mean(hh,:) = mean(SNR_noisy);

end
%%
% 
matrix_plot_2 =1000*[BFS_noisy_error_mean  BFS_filtered_error_mean ]; %最后挑选了其中一次计算的结果画图
%%








