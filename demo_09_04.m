
%%  %用来验证级联数量的影响  也需要每一个num_2都重复20次

clc;
clear all;
close all;
%% 生成20个谱  
fn = 201;         %老师，这里我改了  原来您写的是141  
delatvb = 0.04;   %线宽
% vb_single = [(repmat(10.68,250,1));(linspace(10.68,10.7,100))';(repmat(10.7,500,1));(linspace(10.7,10.68,100))';(repmat(10.68,250,1))];
vb_single =linspace(10.69,10.71,201)';
temp=linspace(10.71,10.69,201)';
vb_single =[vb_single;temp(2:end)];
vb = repmat(vb_single,1,1);   % vb 

Fres = linspace(10.66,10.74,fn); 
% N=round(delatvb/2/(Fres(2)-Fres(1)));%半峰一半对应点数
N=40;
Fres = Fres';
snr_given =10;                   %给定信噪比
num_2_matrix = [4;5;6;7;8;9;10;15;20;25;30;35;40;45;50 ];                     % 重复次数

% num_2_matrix = [4;5;6 ];  
num_1  = size(vb,1);     % 空间采样点数目  1200
for j= 1:num_1
Gain_clean(:,j) = 1./(1+4*((Fres- vb(j))./delatvb).^2);   % 1*fn
end              %生成一个干净的谱   fn行  num_1列（100行  1200列）
M = 10;
mu  = 0.008;

num_repeat = 30;

   BFS_error_mean = zeros(size(num_2_matrix,1),max(num_2_matrix));
    SNR_filtered_mean = zeros(size(num_2_matrix,1),max(num_2_matrix));
    

for hh = 1:size(num_2_matrix,1)
    BFS_filtered_error = zeros(num_repeat,num_2_matrix(hh));
    SNR_filtered_record = zeros(num_repeat,num_2_matrix(hh));
 
    for tt = 1:num_repeat
        
        num_2 = num_2_matrix(hh);
        Gain_clean_series = repmat(Gain_clean,1,num_2);  %20个干净的谱 横向串联
        Gain_noisy_series = zeros(size(Gain_clean_series));
        for i = 1:fn
            Gain_noisy_series(i,:)= awgn(Gain_clean_series(i,:),snr_given,'measured');
        end               %串联后的谱按行加噪
        %%   定义lms滤波器参数
        period = 1;                                 % 参考输入信号 延迟是1
        x_matrix = [zeros(fn,period) Gain_noisy_series];  % 延迟体现为补零
        d_matrix = [ Gain_noisy_series zeros(fn,period)];
        y_out = zeros(size(x_matrix));
        e_out = zeros(size(x_matrix));
        for i = 1:fn
            x  = x_matrix(i,:) ;     % Input to the filter
            d  = d_matrix(i,:);      % Desired signal  延迟了一步
            ha = adaptfilt.lms(M,mu);
            [y,e] = filter(ha,x,d);
            y_out(i,:) =  y;         %滤波输出矩阵
            e_out(i,:)=  e;
        end
        %%  含噪的20个谱各自进行提取频移
        num_4 =num_2;
        
        %% 滤波后的3到20个
        y_out_temp  = y_out(:,2:end);   %去掉第一列以后  是滤波得出的20个谱  因为之前输入信号也补了一个0
        y_out_object  = y_out_temp (:,(1-1)*num_1+1:end);  % 取出滤波后的第1到20个谱
        
        y_out_object(:,1)=y_out(:,3);%我改了
       % BFS_error_mean = zeros(size(num_2_matrix,1),max(num_2_matrix));
       % SNR_filtered_mean = zeros(size(num_2_matrix,1),max(num_2_matrix));
%         BFS_filtered = zeros(num_1,num_4);
        
%         SNR_filtered = zeros(fn,num_4);
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
        %         temp=(mean(abs(BFS_filtered-repmat(vb,1,num_4))));%我改了，看从哪里开始稳定了
        
        
        BFS_filtered_error(tt,:) =  mean(abs(BFS_filtered-repmat(vb,1,num_4)));
        SNR_filtered_record(tt,:) = mean(SNR_filtered);
         disp(['num_4为' num2str(num_2_matrix(hh))  '  第' num2str(tt) '次计算'])
    end
    BFS_error_mean(hh,1:num_2_matrix(hh)) = mean(BFS_filtered_error);
    SNR_filtered_mean(hh,1:num_2_matrix(hh))= mean(SNR_filtered_record);
    
    
end

% return
%% 
figure;
BFS_error_mean_object = [BFS_error_mean(2,:);BFS_error_mean(7:end,:)];

BFS_error_mean_object_small = BFS_error_mean_object(1:2:end,:);
linestyle = {'-s','-*',':','--',':','o','*','-.','-p','-h'};
  c = [1 0 0
      0 1 0
      0 0 1
      0 1 1
      1 0 1]; 

for m = 1:size(BFS_error_mean_object_small ,1)

plot( nonzeros(BFS_error_mean_object_small (m,:)),linestyle{m},'color',c(m,:));
hold on;
end
legend('级联5个谱','级联15个谱','级联25个谱','级联35个谱','级联45个谱');
xlabel('迭代谱个数');
ylabel('布里渊频移误差');
set(gca,'FontSize',16); 


%% 
figure;
 SNR_filtered_mean_object = [ SNR_filtered_mean(2,:); SNR_filtered_mean(7:end,:)];

 SNR_filtered_mean_object_small =  SNR_filtered_mean_object(1:2:end,:);
linestyle = {'s','d',':','+','.','o','*','-.','-p','-h'};
  c = [1 0 1;
      1 0 0;
      0 1 0;
      0 0 1;
      0 1 1;
       ]; 

for m = 1:size( SNR_filtered_mean_object_small  ,1)

plot( nonzeros( SNR_filtered_mean_object_small  (m,:)),linestyle{m},'color',c(m,:));
hold on;
end
legend('级联5个谱','级联15个谱','级联25个谱','级联35个谱','级联45个谱');
xlabel('迭代谱个数');
ylabel('信噪比');
set(gca,'FontSize',16); 
%% 





















linestyle={'-','--',':'};
hold on;
for i=1:length(Rb)
    Ps=[qfuncinv(ber)]^2*q*Rb(i)./(4*r*(1-k));
    plot(k,Ps,linestyle{i},'linewidth',2);
    grid on;
    hold on;
    legend('Rb=1Gbps','Rb=5Gbps','Rb=10Gbps');
end

figure;
k=(0.1:0.01:0.5);
r=1;
ber=1e-9;
q=1.6e-19;
%Q=qfuncinv(ber);
Rb=[1e9,5e9,1e10];
linestyle={'-','--',':'};
hold on;
for i=1:length(Rb)
    Ps=[qfuncinv(ber)]^2*q*Rb(i)./(4*r*(1-k));
    plot(k,Ps,linestyle{i},'linewidth',2);
    grid on;
    hold on; 
    legend('Rb=1Gbps','Rb=5Gbps','Rb=10Gbps');
end

%% 

figure;
qpskConstellation = [-1+1i 1+1i; -1-1i 1-1i]/sqrt(2);
qpsk = reshape(qpskConstellation,1,[]); 
Num  = 40;
outter = 60;
 for nn = 1:outter
 qpsk = qpsk * (outter-1)/outter;
  c = rand(Num,3);       %???ú?ú????12????????RGB???ú??
     for idx = 1:Num
         theta = pi/2/Num*idx;
         rou = [cos(theta) sin(theta);sin(theta) -cos(theta)];
     realPart = real(qpsk);
     imagPart = imag(qpsk);
     reim = rou * [realPart;imagPart];
     realPart2 = real(qpsk*0.3);
     imagPart2 = imag(qpsk*0.3);
     reim2 = rou * [realPart2;imagPart2]; 
     plot(reim(1,:),reim(2,:),'o','color',c(idx,:));
     hold on;
     plot(reim2(1,:),reim2(2,:),'.','color',c(idx,:));
     hold on;
     pause(0.005);
     end
 end

%%

































