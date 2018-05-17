function [QRS_amp,QRS_ind] = DS_detect( ecg_i,gr )
%% function [QRS_amp,QRS_ind]=DS_detect(ecg_i,gr)
%% 输入
% ecg_i : 原信号，一维向量
% gr : 绘图与否，0：不绘图，1：绘图
%% 输出
% QRS_amp:QRS波振幅.
% QRS_ind:QRS波索引.
% 绘制图像.
%% 作者：刘文涵（WHliu@whu.edu.cn）
%% 版本：1.1
%% 04-24-2018
%% 代码:
if nargin < 2
    gr = 1; 
    if nargin<1
           error('The algorithm need a input:ecg_i.');
    end
end
if ~isvector(ecg_i)
  error('ecg_i must be a row or column vector.');
end
fs=360;
if size(ecg_i,2)<round(1.5*fs)+1
    error('The algorithm need a longer input.');
end
tic,
s=ecg_i;
N=size(s,2);
ECG=s;
FIR_c1=[0.0041,0.0053,0.0068,0.0080,0.0081,0.0058,-0.0000,-0.0097,-0.0226,...   
   -0.0370,-0.0498,-0.0577,-0.0576,-0.0477,-0.0278,0,0.0318,0.0625,0.0867,...    
    0.1000,0.1000,0.0867,0.0625,0.0318,0,-0.0278,-0.0477,-0.0576,-0.0577,...   
    -0.0498,-0.0370,-0.0226,-0.0097,-0.0000,0.0058,0.0081,0.0080,0.0068,...
    0.0053,0.0041]; % 使用fdatool设计并导出的滤波器系数,带通FIR,15~25Hz,详情使用fdatool打开DS1.fda查看
FIR_c2=[0.0070,0.0094,0.0162,0.0269,0.0405,0.0555,0.0703,0.0833,0.0928,...    
    0.0979,0.0979,0.0928,0.0833,0.0703,0.0555,0.0405,0.0269,0.0162,0.0094,...    
    0.0070]; % 使用fdatool设计并导出的滤波器系数,低通FIR,截止频率5Hz,详情使用fdatool打开DS2.fda查看

l1=size(FIR_c1,2);
ECG_l=[ones(1,l1)*ECG(1) ECG ones(1,l1)*ECG(N)]; % 数据点延拓，防止滤波边缘效应；
ECG=filter(FIR_c1,1,ECG_l); % 使用filter滤波；
ECG=ECG((l1+1):(N+l1)); % 前面延拓了数据点，这里截取有用的部分；

%% 双斜率处理
a=round(0.015*fs);  % 两侧目标区间0.015~0.060s;
b=round(0.060*fs);
Ns=N-2*b;           % 确保在不超过信号长度；
S_l=zeros(1,b-a+1);
S_r=zeros(1,b-a+1);
S_dmax=zeros(1,Ns);
for i=1:Ns          % 对每个点双斜率处理；
    for k=a:b
        S_l(k-a+1)=(ECG(i+b)-ECG(i+b-k))./k;
        S_r(k-a+1)=(ECG(i+b)-ECG(i+b+k))./k;
    end
  S_lmax=max(S_l);
  S_lmin=min(S_l);
  S_rmax=max(S_r);
  S_rmin=min(S_r);
  C1=S_rmax-S_lmin;
  C2=S_lmax-S_rmin;
  S_dmax(i)=max([C1 C2]);
end

%% 再次进行低通滤波，思路与上述带通滤波一致
l2=size(FIR_c2,2);
S_dmaxl=[ones(1,l2)*S_dmax(1) S_dmax ones(1,l2)*S_dmax(Ns)];
S_dmaxt=filter(FIR_c2,1,S_dmaxl);
S_dmaxt=S_dmaxt((l2+1):(Ns+l2));

%% 滑动窗口积分
w=8;wd=7;
d_l=[zeros(1,w) S_dmaxt zeros(1,w)];  % 零延拓，确保所有的点都可以进行窗口积分
m=zeros(1,Ns);
   for n=(w+1):(Ns+w)                 % 滑动窗口；
      m(n-w)=sum(d_l(n-w:n+w));       % 积分；
   end
m_l=[ones(1,wd)*m(1) m ones(1,wd)*m(Ns)]; 

%% 双阈值检测与动态变化
QRS_buf1=[];   % 存储检测到的QRS波索引
AMP_buf1=[];   % 存储最近检测到的8个QRS波对应特征信号的波峰值
thr_init0=0.4;thr_lim0=0.23;
thr_init1=0.6;thr_lim1=0.3;  %% 阈值变化的初始值和下限设置
en=-1;        % 标记波峰检出情况，高于高阈值--1，高低阈值之间--0，未检出-- -1
thr0=thr_init0;
thr1=thr_init1;
thr1_buf=[]; % 阈值缓存，记录阈值变化情况；
thr0_buf=[];
for j=8:Ns
       t=1;
       cri=1;
       while t<=wd&&cri>0   % 检测候选波峰；
           cri=((m_l(j)-m_l(j-t))>0)&&(m_l(j)-m_l(j+t)>0);
           t=t+1;
       end
       if t==wd+1
           N1=size(QRS_buf1,2);               %N1:已经检测到的QRS波个数
           if m_l(j)>thr1                     % 高于高阈值时的处理
               if N1<2                        % N1小于2时直接存储；
                 QRS_buf1=[QRS_buf1 (j-wd)];  % j-wd 减去了滑动窗口积分带来的延迟；
                 AMP_buf1=[AMP_buf1 m_l(j)];
                 en=1;
               else
                 dist=j-wd-QRS_buf1(N1);
                 if dist>0.24*fs               % 检测波峰距离；
                     QRS_buf1=[QRS_buf1 (j-wd)]; 
                     AMP_buf1=[AMP_buf1 m_l(j)];
                     en=1;
                 else
                     if m_l(j)>AMP_buf1(end)   % 不应期处理
                         QRS_buf1(end)=j-wd;
                         AMP_buf1(end)=m_l(j);
                         en=1;
                     end     
                 end
               end
     
          else                                 % 特征峰值低于高阈值
               
              if N1<2&&m_l(j)>thr0             % 特征峰值在两阈值之间
                  QRS_buf1=[QRS_buf1 (j-wd)];
                  AMP_buf1=[AMP_buf1 m_l(j)];
                  en=0;
              else
                if m_l(j)>thr0                 % 特征峰值在两阈值之间
                  dist_m=mean(diff(QRS_buf1));
                  dist=j-wd-QRS_buf1(N1);
                  if dist>0.24*fs && dist>0.5*dist_m  % 不应期检测，并且，波峰要距离足够远（> 平均距离的一半）
                     QRS_buf1=[QRS_buf1 (j-wd)];
                     AMP_buf1=[AMP_buf1 m_l(j)];
                     en=0;
                  else
                      if m_l(j)>AMP_buf1(end)
                         QRS_buf1(end)=j-wd;
                         AMP_buf1(end)=m_l(j);
                         en=0;
                      end 
                  end
                else
                    en=-1;
                end
              end
           end
           N2=size(AMP_buf1,2);
           if N2>8
               AMP_buf1=AMP_buf1(2:9); % 确保只存储最近的8个特征波峰；
           end
		   % 下面的if与博文中的公式对应
           if en==1
              thr1=0.7*mean(AMP_buf1);
              thr0=0.25*mean(AMP_buf1);
           else
               if en==0
                   thr1=thr1-(abs(m_l(j)-mean(AMP_buf1)))/2;
                   thr0=0.4*m_l(j);
               end
           end
       end
       if thr1<=thr_lim1   % 确保阈值高于下限
           thr1=thr_lim1;
       end
       
       if thr0<=thr_lim0
           thr0=thr_lim0;
       end
       
      thr1_buf=[thr1_buf thr1]; 
      thr0_buf=[thr0_buf thr0];
end
delay=round(l1/2)-2*w+2;
QRS_ind=QRS_buf1-delay;   % 减去延迟，得到最终结果；
QRS_amp=s(QRS_ind);
toc
if gr==1    %绘图
   subplot(2,1,1);plot(m);axis([1 size(m,2) -0.3 1.6*max(m)]);
   hold on;title('Feature signal and thresholds');grid on;
   plot(QRS_buf1,m(QRS_buf1),'ro');
   plot(thr1_buf,'r');
   plot(thr0_buf,'k');
   legend('Feature Signal','QRS Locations','Threshold1','Threshold0');
   subplot(2,1,2);plot(s);%axis([1 size(s,2) min(s) 1.5*max(s)]);
   xlabel('n');ylabel('Voltage / mV');
   hold on;title('The result on the raw ECG');grid on;
   plot(QRS_ind,QRS_amp,'ro');
   legend('Raw ECG','QRS Locations');
end
end

