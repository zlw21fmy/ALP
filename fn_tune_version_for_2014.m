%% 这个程序主要用在循环进行优化步长和阶数  0924 新加了画图比较两种信号输入方式的部分，但是没有计算程序
clc;
clear all;
close all;
%% 生成40个横向串联的谱
fn = 201;         %扫频点个数
delatvb = 0.04;   %线宽
vb_single =linspace(10.69,10.71,201)';
temp=linspace(10.71,10.69,201)';
vb_single =[vb_single;temp(2:end)];
vb = repmat(vb_single,1,1);   % vb

Fres = linspace(10.66,10.74,fn);
N=round(delatvb/2/(Fres(2)-Fres(1)));%半峰一半对应点数

Fres = Fres';

num_2 =40;                      % 串联次数
num_4 = num_2 ;
%%
snr_given =10;                   %给定信噪比


%%
num_1  = size(vb,1);     % 空间采样点数目 401


for j= 1:num_1
    Gain_clean(:,j) = 1./(1+4*((Fres- vb(j))./delatvb).^2);   % 1*fn
end              %生成一个干净的谱   fn行  num_1列

Gain_clean_series = repmat(Gain_clean,1,num_2);  %40个干净的谱 横向串联

%%
M_matrix = [5;6;7;8;9;10;15;20;25;30;35;40;45;50;55;60]  ;            %16*1    % LMS 滤波器阶数
mu_matrix =[0.0005;0.001;0.003;0.005;0.007;0.008;0.01;0.015;0.02;0.025;0.03]'; %11*1  % 步长
for jj = 1:length(M_matrix)
    M = M_matrix(jj);
    for kk = 1:length(mu_matrix)
        
        mu  = mu_matrix(kk);
        for  ii = 1:20
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
                %                 disp(['第' num2str(i) '个频率点去噪']);
            end
            
            %% 滤波后的1到40个
            y_out_temp  = y_out(:,2:end);   %去掉第一列以后  是滤波得出的20个谱  因为之前输入信号也补了一个0
            y_out_object  = y_out_temp (:,(1-1)*num_1+1:end);  % 取出滤波后的第1到20个谱
            y_out_object(:,1)=y_out(:,3);%我改了
            
            for p = 1:num_4    %我改了
                y_test = y_out_object(:,(p-1)*num_1+1:(p)*num_1);  %选出滤波后的第p个谱
                [m1, locs1] = max(y_test, [], 1);  % 找出各列的最大值
                for qq = 1:num_1
                    locs1_object = locs1(qq);
                    
                    if  locs1_object> N && locs1_object<=fn- N
                        indice1 = locs1_object-(N):locs1_object+(N); %从最大值从左到右各找半个线宽
                    elseif locs1_object<= N
                        indice1 = 1:locs1_object+(locs1_object-1); %从最大值从左到右各找半个线宽
                        
                    elseif locs1_object> fn-N
                        indice1 = locs1_object-(fn-locs1_object) :fn; %从最大值从左到右各找半个线宽
                        
                    end
                    Fres1_object = Fres(indice1);                  % 找到一个线宽对应的扫频区间
                    p1 = polyfit(Fres1_object, y_test(indice1,qq),2); % 二次多项式拟合
                    BFS_filtered(qq,p) = -p1(2)./(2*p1(1));%利用抛物线的对称轴计算中心频移，处理后
                    
                end
                for b = 1:fn
                    SNR_filtered(b,p) = snr(y_test(b,:),abs(y_test(b,:)-Gain_clean(b,:)));
                end
            end
            %%  按 第1 种计算方式  最小误差为0.47MHz
            BFS_filtered_error_1(ii,:)=(mean(abs(BFS_filtered-repmat(vb,1,num_4))));%各个频差点处取平均
            BFS_filtered_error_1_mean = mean(BFS_filtered_error_1(ii,:));
           
            disp([ '阶数' num2str(M_matrix(jj))  '步长' num2str(mu_matrix(kk)) ]);
            disp(['第' num2str(ii) '次计算'   '第一种衡量方式下  滤波后频移结果为' num2str(BFS_filtered_error_1_mean)  ]);
       
        end
        BFS_error_average_1(jj,kk) = mean(mean(BFS_filtered_error_1)); 
    end
end



