%%  ͼ14a  14b            �˲��������ֱ�Ӻ����źŵ��ӵ������  Ƶ������ͼ   
clc;
clear all;
close all;
%% ����20����  
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
snr_given_matrix =5:5:20;  %���������
for j= 1:num_1
Gain_clean(:,j) = 1./(1+4*((Fres- vb(j))./delatvb).^2);   % 1*fn
end              %����һ���ɾ�����   fn��  num_1�У�100��  1200�У�

Gain_clean_series = repmat(Gain_clean,1,num_2);  %20���ɾ����� ������
%% 
M = 15;
mu  = 0.005;
num_repeat = 5; 
%% �˲����źŵ���

for hh = 1:length(snr_given_matrix)
    snr_given = snr_given_matrix(hh);
for  tt = 1:num_repeat
    Gain_noisy_series =  zeros(fn,num_1*num_2);
    
    
    for i = 1:fn
        Gain_noisy_series(i,:)= awgn(Gain_clean_series(i,:),snr_given,'measured');
    end               %��������װ��м���
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
    y_out_object  = y_out_temp (:,(1-1)*num_1+1:end);  % �ӵ������׿�ʼ
      y_out_object(:,1)=y_out(:,3);
    for p = 1:(num_4)
        
        y_test = y_out_object(:,(p-1)*num_1+1:(p)*num_1);  %ѡ���˲���ĵ�p����
        [m1 ,locs1] = max(y_test, [], 1);  % �ҳ����е����ֵ
        
        for jj = 1:num_1
            locs1_object = locs1(jj);
            if  locs1_object> N && locs1_object<=fn- N
                indice1 = locs1_object-(N):locs1_object+(N); %�����ֵ�����Ҹ��Ұ���߿�
            elseif locs1_object<= N
                indice1 = 1:locs1_object+(locs1_object-1); %�����ֵ�����Ҹ��Ұ���߿�
            elseif locs1_object> fn-N
                indice1 = locs1_object-(fn-locs1_object) :fn; %�����ֵ�����Ҹ��Ұ���߿�
            end
            
            Fres1_object = Fres(indice1);
            p1 = polyfit(Fres1_object, y_test(indice1,jj),2);
            BFS_filtered_enssemble(jj,p) = -p1(2)./(2*p1(1));%���������ߵĶԳ����������Ƶ��
        end
        
        for ii = 1:fn
            SNR_filtered_enssemble(ii,p) = snr(y_test(ii,:),abs(y_test(ii,:)-Gain_clean(ii,:)));
            
        end
    end
    
    BFS_filtered_error(tt,:) = mean(abs( BFS_filtered_enssemble-repmat(vb,1,(num_4))));
    SNR_filtered_mean(tt,:) = mean(SNR_filtered_enssemble);
    disp([ '�����Ϊ' num2str(snr_given_matrix(hh)) '  �˲���' num2str(tt) '�μ���'])
end
BFS_filtered_error_mean(hh,:) = mean(BFS_filtered_error);    
SNR_filtered_mean_mean(hh,:) = mean(SNR_filtered_mean);


end

%% �����ź�ֱ�ӵ��� 

for hh = 1:length(snr_given_matrix)
    snr_given = snr_given_matrix(hh);
for  tt = 1:num_repeat
    Gain_noisy_series = zeros(fn,num_1*num_2);
    Gain_noisy_3d  = zeros(fn,num_1,num_2);
    for i = 1:fn
        Gain_noisy_series(i,:)= awgn(Gain_clean_series(i,:),snr_given,'measured');
    end               %��������װ��м���
    %%  ����������ε��� ���Խ�����ȡƵ��
    num_4 = num_2;
    for p = 1:num_4
        Gain_noisy_3d(:,:,p) = (Gain_noisy_series(:,(p-1)*num_1+1:(p)*num_1)); %Gain_noisy_series ��50���׺���ƴ�� Gain_noisy_3d �൱��50 ҳ��ÿһҳ��һ���� ��������Ϊ�˷������
        Gain_average(:,(p-1)*num_1+1:(p)*num_1) = sum(Gain_noisy_3d,3)/p;  %Gain_noisy_3d ��������ά�ȵķ������ ��ȡƽ��  Gain_average �ѵ���ƽ����Ľ����50���ף�����ƴ��
        
        gain_noisy_test = Gain_average(:,(p-1)*num_1+1:(p)*num_1);  %% ������ѡ���ĵ�p������
        [m1 ,locs1] = max(gain_noisy_test, [], 1);  % �ҳ����е����ֵ  ������һ��λ�ô������������ߵķ�ֵ
        
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
            BFS(jj,p) = -p1(2)./(2*p1(1));%���������ߵĶԳ����������Ƶ��
            
        end
        for ii = 1:fn
            SNR(ii,p) = snr(gain_noisy_test(ii,:),abs(gain_noisy_test(ii,:)-Gain_clean(ii,:)));
        end
    end
    BFS_noisy_error(tt,:) =  mean(abs(BFS-repmat(vb,1,num_4)));
    SNR_noisy(tt,:) = mean(SNR);
    disp([ '�����Ϊ' num2str(snr_given_matrix(hh))  ' �����źŵ�' num2str(tt) '�μ���'])
end
BFS_noisy_error_mean(hh,:) = mean(BFS_noisy_error);    
SNR_noisy_mean_mean(hh,:) = mean(SNR_noisy);

end



%% ��ͼ14a
matrix_plot_1 = 1000*[BFS_noisy_error_mean(2,:);BFS_filtered_error_mean(2,:)]';
%% ��ͼ14b
matrix_plot_2 = [SNR_noisy_mean_mean(2,:);SNR_filtered_mean_mean(2,:)]';





