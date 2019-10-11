%% ���������Ҫ����ѭ�������Ż������ͽ���  0924 �¼��˻�ͼ�Ƚ������ź����뷽ʽ�Ĳ��֣�����û�м������
clc;
clear all;
close all;
%% ����40������������
fn = 201;         %ɨƵ�����
delatvb = 0.04;   %�߿�
vb_single =linspace(10.69,10.71,201)';
temp=linspace(10.71,10.69,201)';
vb_single =[vb_single;temp(2:end)];
vb = repmat(vb_single,1,1);   % vb

Fres = linspace(10.66,10.74,fn);
N=round(delatvb/2/(Fres(2)-Fres(1)));%���һ���Ӧ����

Fres = Fres';

num_2 =40;                      % ��������
num_4 = num_2 ;
%%
snr_given =10;                   %���������


%%
num_1  = size(vb,1);     % �ռ��������Ŀ 401


for j= 1:num_1
    Gain_clean(:,j) = 1./(1+4*((Fres- vb(j))./delatvb).^2);   % 1*fn
end              %����һ���ɾ�����   fn��  num_1��

Gain_clean_series = repmat(Gain_clean,1,num_2);  %40���ɾ����� ������

%%
M_matrix = [5;6;7;8;9;10;15;20;25;30;35;40;45;50;55;60]  ;            %16*1    % LMS �˲�������
mu_matrix =[0.0005;0.001;0.003;0.005;0.007;0.008;0.01;0.015;0.02;0.025;0.03]'; %11*1  % ����
for jj = 1:length(M_matrix)
    M = M_matrix(jj);
    for kk = 1:length(mu_matrix)
        
        mu  = mu_matrix(kk);
        for  ii = 1:20
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
                %                 disp(['��' num2str(i) '��Ƶ�ʵ�ȥ��']);
            end
            
            %% �˲����1��40��
            y_out_temp  = y_out(:,2:end);   %ȥ����һ���Ժ�  ���˲��ó���20����  ��Ϊ֮ǰ�����ź�Ҳ����һ��0
            y_out_object  = y_out_temp (:,(1-1)*num_1+1:end);  % ȡ���˲���ĵ�1��20����
            y_out_object(:,1)=y_out(:,3);%�Ҹ���
            
            for p = 1:num_4    %�Ҹ���
                y_test = y_out_object(:,(p-1)*num_1+1:(p)*num_1);  %ѡ���˲���ĵ�p����
                [m1, locs1] = max(y_test, [], 1);  % �ҳ����е����ֵ
                for qq = 1:num_1
                    locs1_object = locs1(qq);
                    
                    if  locs1_object> N && locs1_object<=fn- N
                        indice1 = locs1_object-(N):locs1_object+(N); %�����ֵ�����Ҹ��Ұ���߿�
                    elseif locs1_object<= N
                        indice1 = 1:locs1_object+(locs1_object-1); %�����ֵ�����Ҹ��Ұ���߿�
                        
                    elseif locs1_object> fn-N
                        indice1 = locs1_object-(fn-locs1_object) :fn; %�����ֵ�����Ҹ��Ұ���߿�
                        
                    end
                    Fres1_object = Fres(indice1);                  % �ҵ�һ���߿��Ӧ��ɨƵ����
                    p1 = polyfit(Fres1_object, y_test(indice1,qq),2); % ���ζ���ʽ���
                    BFS_filtered(qq,p) = -p1(2)./(2*p1(1));%���������ߵĶԳ����������Ƶ�ƣ������
                    
                end
                for b = 1:fn
                    SNR_filtered(b,p) = snr(y_test(b,:),abs(y_test(b,:)-Gain_clean(b,:)));
                end
            end
            %%  �� ��1 �ּ��㷽ʽ  ��С���Ϊ0.47MHz
            BFS_filtered_error_1(ii,:)=(mean(abs(BFS_filtered-repmat(vb,1,num_4))));%����Ƶ��㴦ȡƽ��
            BFS_filtered_error_1_mean = mean(BFS_filtered_error_1(ii,:));
           
            disp([ '����' num2str(M_matrix(jj))  '����' num2str(mu_matrix(kk)) ]);
            disp(['��' num2str(ii) '�μ���'   '��һ�ֺ�����ʽ��  �˲���Ƶ�ƽ��Ϊ' num2str(BFS_filtered_error_1_mean)  ]);
       
        end
        BFS_error_average_1(jj,kk) = mean(mean(BFS_filtered_error_1)); 
    end
end



