
%% �����˲�ǰ��Աȵ�ͼ    1��ͼ12  �������ߵĲ���ԨƵ����ȡ��� 2�� ͼ13    �ռ��7��λ�õ�Ĳ���Ԩ��
%                           3��ͼ14  ȥ��ǰ����ԨƵ���������������ı仯

clc;
clear all;
close all;
%% ����40����  
fn = 201;         %
delatvb = 0.04;   %
vb_single =linspace(10.69,10.71,201)';
temp=linspace(10.71,10.69,201)';
vb_single =[vb_single;temp(2:end)];
vb = repmat(vb_single,1,1);   % vb 

Fres = linspace(10.66,10.74,fn); 
% N=round(delatvb/2/(Fres(2)-Fres(1)));%���һ���Ӧ����
N=40;
Fres = Fres';
num_1  = size(vb,1);     % �ռ��������Ŀ  1200

num_2 =40;                      % �ظ�����
snr_given = 10;  %���������
for j= 1:num_1
Gain_clean(:,j) = 1./(1+4*((Fres- vb(j))./delatvb).^2);   % 1*fn
end              %����һ���ɾ�����   fn��  num_1�У�100��  1200�У�

Gain_clean_series = repmat(Gain_clean,1,num_2);  %20���ɾ����� ������
%% 
M = 15;
mu  = 0.005;
 
    for i = 1:fn
        Gain_noisy_series(i,:)= awgn(Gain_clean_series(i,:),snr_given,'measured');
    end               %��������װ��м���
    
    Gain_clean_series_2 = repmat(Gain_clean,num_2,1); 
    
    for j = 1:num_1
        Gain_noisy_series_2(:,j)= awgn(Gain_clean_series_2(:,j),snr_given,'measured');
    end  
    
    
    %%   ����lms�˲�������
    period = 1;                                 % �ο������ź� �ӳ���1
    x_matrix = [zeros(fn,period) Gain_noisy_series];  % �ӳ�����Ϊ����
    d_matrix = [ Gain_noisy_series zeros(fn,period)];
    for i = 1:fn
        x  = x_matrix(i,:) ;     % Input to the filter
        d  = d_matrix(i,:);      % Desired signal  �ӳ���һ��
        ha = adaptfilt.lms(M,mu);
        [y,e] = filter(ha,x,d);
        y_out(i,:) =  y;         %�˲��������
        e_out(i,:)=  e;
    end
    
    num_4 = num_2;
    y_out_temp  = y_out(:,2:end);   %
    
    y_out_object  = y_out_temp (:,(1-1)*num_1+1:end);  % �ӵ�1���׿�ʼ
     y_out_object(:,1)=y_out(:,3);  
    for p = 1:(num_4)
        
        y_test = y_out_object(:,(p-1)*num_1+1:(p)*num_1);  %ѡ���˲���ĵ�p����
        [m1 ,locs1] = max(y_test, [], 1);  % �ҳ����е����ֵ
        
        for jj = 1:num_1
            locs1_object = locs1(jj);
            locs1_object_matrix(p,jj) = locs1_object;
            if  locs1_object> N && locs1_object<=fn- N
                indice1 = locs1_object-(N):locs1_object+(N); %�����ֵ�����Ҹ��Ұ���߿�
            elseif locs1_object<= N
                indice1 = 1:locs1_object+(locs1_object-1); %�����ֵ�����Ҹ��Ұ���߿�
            elseif locs1_object> fn-N
                indice1 = locs1_object-(fn-locs1_object) :fn; %�����ֵ�����Ҹ��Ұ���߿�
            end
            
            Fres1_object = Fres(indice1);
            p1 = polyfit(Fres1_object, y_test(indice1,jj),2);
            BFS_filtered(jj,p) = -p1(2)./(2*p1(1));%���������ߵĶԳ����������Ƶ��
        end
        
        for ii = 1:fn
            SNR_filtered(ii,p) = snr(y_test(ii,:),abs(y_test(ii,:)-Gain_clean(ii,:)));
            
        end
    end
   
    %% ���ں����ź�
    
    for p = 1:num_4
        
        gain_noisy_test = Gain_noisy_series(:,(p-1)*num_1+1:(p)*num_1);  %% ������ѡ���ĵ�p������
        [m1, locs1] = max(gain_noisy_test, [], 1);  % �ҳ����е����ֵ  ������һ��λ�ô������������ߵķ�ֵ
        
        for jj = 1:num_1
            locs1_object = locs1(jj);   % �ҳ����ֵ ��Ӧ��λ��
            
            if  locs1_object> N && locs1_object<=fn- N
                indice1 = locs1_object-(N):locs1_object+(N); %�����ֵ�����Ҹ��Ұ���߿�
            elseif locs1_object<= N
                indice1 = 1:locs1_object+(locs1_object-1); %�����ֵ�����Ҹ��Ұ���߿�
                
            elseif locs1_object> fn-N
                indice1 = locs1_object-(fn-locs1_object) :fn; %�����ֵ�����Ҹ��Ұ���߿�
                
            end
            
            Fres1_object = Fres(indice1);         %�ҵ���Ϸ�Χ��Ӧ��Ƶ�ʵ�ֵ
            p1 = polyfit(Fres1_object,gain_noisy_test(indice1,jj),2);
            BFS_noisy(jj,p) = -p1(2)./(2*p1(1));%���������ߵĶԳ����������Ƶ��
            
        end
        for ii = 1:fn
            SNR_noisy(ii,p) = snr(gain_noisy_test(ii,:),abs(gain_noisy_test(ii,:)-Gain_clean(ii,:)));   
        end
        
    end
    
    
 %%  ���ӻ�   ͼ11  �������ߵĲ���ԨƵ����ȡ���
   
   matrix_for_plot_1 = [ vb BFS_noisy(:,6)   ];     %������ȡ
    
  matrix_for_plot_2 = [ vb BFS_filtered(:,end)   ]; %�˲�����ȡ
%%  ͼ12    �ռ��7��λ�õ�Ĳ���Ԩ��
 
matrix_for_plot_3 = [ Gain_clean(:,7)  gain_noisy_test(:,7)  y_test(:,7)];

%%  ͼ13  ȥ��ǰ����ԨƵ���������������ı仯

    BFS_noisy_error = mean(abs(BFS_noisy-repmat(vb,1,num_4)));
    BFS_filtered_error = mean(abs( BFS_filtered-repmat(vb,1,num_4)));
    matrix_for_plot_4 = 1000*[ BFS_noisy_error;  BFS_filtered_error]' ;
 %%    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    