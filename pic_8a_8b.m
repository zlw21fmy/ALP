
%% 画出图 8（a） 8（b）
clc;
clear all;
close all;
%% 

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

num_2 =40;                      % 重复次数
snr_given = 10;  %给定信噪比
for j= 1:num_1
Gain_clean(:,j) = 1./(1+4*((Fres- vb(j))./delatvb).^2);   % 1*fn
end              %生成一个干净的谱   fn行  num_1列（100行  1200列）

Gain_clean_series_1 = repmat(Gain_clean,1,num_2);  %40个干净的谱 横向串联
Gain_clean_series_2 = repmat(Gain_clean,num_2,1);  %40个干净的谱 纵向串联


    for i = 1:fn
        Gain_noisy_series_1(i,:)= awgn(Gain_clean_series_1(i,:),snr_given,'measured');
    end               %串联后的谱按行加噪
    
    Gain_clean_series_2 = repmat(Gain_clean,num_2,1); 
    
    for j = 1:num_1
        Gain_noisy_series_2(:,j)= awgn(Gain_clean_series_2(:,j),snr_given,'measured');
    end  
    %% 
    figure;
    plot(Gain_noisy_series_1(1,1:900),'b');
    hold on;
    plot(Gain_clean_series_1 (1,1:900),'r','LineWidth',2);
 %%    
 figure;
 plot(Gain_noisy_series_2(1:500,1),'b');
 hold on;
 plot(Gain_clean_series_2 (1:500,1),'r','LineWidth',2);

























