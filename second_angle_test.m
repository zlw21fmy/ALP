%% �������㣻�����ź����뷽ʽ�µ����������
clc;
clear all;
close all;

%% �������ֺ������
 
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

Gain_clean_series_1 = repmat(Gain_clean,1,num_2);  %40���ɾ����� ������
Gain_clean_series_2 = repmat(Gain_clean,num_2,1);  %40���ɾ����� ������


    for i = 1:fn
        Gain_noisy_series_1(i,:)= awgn(Gain_clean_series_1(i,:),snr_given,'measured');
    end               %��������װ��м���
    
    Gain_clean_series_2 = repmat(Gain_clean,num_2,1); 
    
    for j = 1:num_1
        Gain_noisy_series_2(:,j)= awgn(Gain_clean_series_2(:,j),snr_given,'measured');
    end  
%%    ���ڵ�һ���˲���ʽ   ����lms�˲�������

M = 15 ;                                    % LMS �˲�������
mu = 0.005;                                 % ����
period = 1;                                 % �ο������ź� �ӳ���1
    x_matrix = [zeros(fn,period) Gain_noisy_series_1];  % �ӳ�����Ϊ����
    d_matrix = [ Gain_noisy_series_1 zeros(fn,period)];
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
            BFS_filtered_1(jj,p) = -p1(2)./(2*p1(1));%���������ߵĶԳ����������Ƶ��
        end
        
        for ii = 1:fn
            SNR_filtered_1(ii,p) = snr(y_test(ii,:),abs(y_test(ii,:)-Gain_clean(ii,:)));
            
        end
    end
   BFS_filtered_1_error = mean(mean(abs(BFS_filtered_1-repmat(vb,1,num_2))))
   SNR_filtered_1_result = mean(mean( SNR_filtered_1))
    %% ���ں����ź�
    
    for p = 1:num_4
        
        gain_noisy_test = Gain_noisy_series_1(:,(p-1)*num_1+1:(p)*num_1);  %% ������ѡ���ĵ�p������
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
            BFS_noisy_1(jj,p) = -p1(2)./(2*p1(1));%���������ߵĶԳ����������Ƶ��
            
        end
        for ii = 1:fn
            SNR_noisy_1(ii,p) = snr(gain_noisy_test(ii,:),abs(gain_noisy_test(ii,:)-Gain_clean(ii,:)));   
        end
        
    end
  BFS_noisy_1_error = mean(mean(abs( BFS_noisy_1-repmat(vb,1,num_2))))
   SNR_noisy_1_result = mean(mean(SNR_noisy_1))

%%    ���ڵ�2���˲���ʽ   ����lms�˲�������

period = 1;                                 % �ο������ź� �ӳ���1
    x_matrix = [zeros(period,num_1); Gain_noisy_series_2];  % �ӳ�����Ϊ����
    d_matrix = [Gain_noisy_series_2; zeros(period,num_1)];
    for j = 1:num_1
        x  = x_matrix(:,j) ;     % Input to the filter
        d  = d_matrix(:,j);      % Desired signal  �ӳ���һ��
        ha = adaptfilt.lms(M,mu);
        [y,e] = filter(ha,x,d);
        y_out_2(:,j) =  y;         %�˲��������
        e_out_2(:,j)=  e;
    end
 
    
    
   y_out_object_2  = y_out_2(2:end,:);  % �ӵ�1���׿�ʼ  
    y_out_object_2(1,:)= y_out_object_2(2,:)
 for p = 1:(num_4)
        
        y_test_2 = y_out_object_2((p-1)*fn+1:(p)*fn,:);  %ѡ���˲���ĵ�p����
        [m1 ,locs1] = max(y_test_2, [], 1);  % �ҳ����е����ֵ
        
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
            p1 = polyfit(Fres1_object, y_test_2(indice1,jj),2);
            BFS_filtered_2(jj,p) = -p1(2)./(2*p1(1));%���������ߵĶԳ����������Ƶ��
        end
        
        for ii = 1:fn
            SNR_filtered_2(ii,p) = snr(y_test_2(ii,:),abs(y_test_2(ii,:)-Gain_clean(ii,:)));
            
        end
    end

BFS_filtered_2_error = mean(mean(abs( BFS_filtered_2-repmat(vb,1,num_2))))
SNR_filtered_2_result = mean(mean(SNR_filtered_2))

  %% �ڶ����˲���ʽ��   ���ں����ź�
%     
%     for p = 1:num_4
%         
%         gain_noisy_test_2 = Gain_noisy_series_2((p-1)*fn+1:(p)*fn,:);  %% ������ѡ���ĵ�p������
%         [m1, locs1] = max(gain_noisy_test_2, [], 1);  % �ҳ����е����ֵ  ������һ��λ�ô������������ߵķ�ֵ
%         
%         for jj = 1:num_1
%             locs1_object = locs1(jj);   % �ҳ����ֵ ��Ӧ��λ��
%             
%             if  locs1_object> N && locs1_object<=fn- N
%                 indice1 = locs1_object-(N):locs1_object+(N); %�����ֵ�����Ҹ��Ұ���߿�
%             elseif locs1_object<= N
%                 indice1 = 1:locs1_object+(locs1_object-1); %�����ֵ�����Ҹ��Ұ���߿�
%                 
%             elseif locs1_object> fn-N
%                 indice1 = locs1_object-(fn-locs1_object) :fn; %�����ֵ�����Ҹ��Ұ���߿�
%                 
%             end
%             
%             Fres1_object = Fres(indice1);         %�ҵ���Ϸ�Χ��Ӧ��Ƶ�ʵ�ֵ
%             p1 = polyfit(Fres1_object,gain_noisy_test_2(indice1,jj),2);
%             BFS_noisy_2(jj,p) = -p1(2)./(2*p1(1));%���������ߵĶԳ����������Ƶ��
%             
%         end
%         for ii = 1:fn
%             SNR_noisy_2(ii,p) = snr(gain_noisy_test_2(ii,:),abs(gain_noisy_test_2(ii,:)-Gain_clean(ii,:)));   
%         end
%         
%     end
    
    
% BFS_noisy_2_error = mean(mean(abs( BFS_noisy_2-repmat(vb,1,num_2))))
% SNR_noisy_2_result = mean(mean(SNR_noisy_2))
%% 
SNR_noisy_2_result= SNR_noisy_1_result;
BFS_noisy_2_error = BFS_noisy_1_error;
format long
matrix_for_out = [SNR_noisy_1_result SNR_filtered_1_result;
                  BFS_noisy_1_error  BFS_filtered_1_error;
                  SNR_noisy_2_result  SNR_filtered_2_result
                  BFS_noisy_2_error  BFS_filtered_2_error]
































