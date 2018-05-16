%% 截取全数据库目标心拍
clc; clear all;
Name_whole=[100,101,102,103,104,105,106,107,108,109,111,112,113,114,115,...
    116,117,118,119,121,122,123,124,200,201,202,203,205,207,208,209,...
    210,212,213,214,215,217,219,220,221,222,223,228,230,231,232,233,234];
INFO=[];
Nb=[];Lb=[];Rb=[];Vb=[];Ab=[];Sb=[];
Pl=100;Pr=150;
for na=1:48
    
%------ SPECIFY DATA ------------------------------------------------------
%------ 指定数据文件 -------------------------------------------------------

Name=num2str(Name_whole(na));
PATH= 'F:\MATLAB\MIT-BIH'; % 指定数据的储存路径  
HEADERFILE=strcat(Name, '.hea');% .hea 格式，头文件，可用记事本打开        
ATRFILE=strcat(Name, '.atr'); % .atr 格式，属性文件，数据格式为二进制数    
DATAFILE=strcat(Name, '.dat');% .dat 格式，ECG 数据
SAMPLES2READ=650000;          % 指定需要读入的样本数
                            % 若.dat文件中存储有两个通道的信号:
                            % 则读入 2*SAMPLES2READ 个数据 

%------ LOAD HEADER DATA --------------------------------------------------
%------ 读入头文件数据 -----------------------------------------------------
%
% 示例：用记事本打开的117.hea 文件的数据
%
%      117 2 360 650000
%      117.dat 212 200 11 1024 839 31170 0 MLII
%      117.dat 212 200 11 1024 930 28083 0 V2
%      # 69 M 950 654 x2
%      # None
%
%-------------------------------------------------------------------------
fprintf(1,'\\n$> WORKING ON %s ...\n', HEADERFILE); % 在Matlab命令行窗口提示当前工作状态
% 
% 【注】函数 fprintf 的功能将格式化的数据写入到指定文件中。
% 表达式：count = fprintf(fid,format,A,...)
% 在字符串'format'的控制下，将矩阵A的实数数据进行格式化，并写入到文件对象fid中。该函数返回所写入数据的字节数 count。
% fid 是通过函数 fopen 获得的整型文件标识符。fid=1，表示标准输出（即输出到屏幕显示）；fid=2，表示标准偏差。
%
signalh= fullfile(PATH, HEADERFILE);    % 通过函数 fullfile 获得头文件的完整路径
fid1=fopen(signalh,'r');    % 打开头文件，其标识符为 fid1 ，属性为'r'--“只读”
z= fgetl(fid1);             % 读取头文件的第一行数据，字符串格式
A= sscanf(z, '%*s %d %d %d',[1,3]); % 按照格式 '%*s %d %d %d' 转换数据并存入矩阵 A 中
nosig= A(1);    % 信号通道数目
sfreq=A(2);     % 数据采样频率
clear A;        % 清空矩阵 A ，准备获取下一行数据
for k=1:nosig           % 读取每个通道信号的数据信息
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
    dformat(k)= A(1);           % 信号格式; 这里只允许为 212 格式
    gain(k)= A(2);              % 每 mV 包含的整数个数
    bitres(k)= A(3);            % 采样精度（位分辨率）
    zerovalue(k)= A(4);         % ECG 信号零点相应的整数值
    firstvalue(k)= A(5);        % 信号的第一个整数值 (用于偏差测试)
end;
fclose(fid1);
clear A;

%------ LOAD BINARY DATA --------------------------------------------------
%------ 读取 ECG 信号二值数据 ----------------------------------------------
%
if dformat~= [212,212], error('this script does not apply binary formats different to 212.'); end;
signald= fullfile(PATH, DATAFILE);            % 读入 212 格式的 ECG 信号数据
fid2=fopen(signald,'r');
A= fread(fid2, [3, SAMPLES2READ], 'uint8')';  % matrix with 3 rows, each 8 bits long, = 2*12bit
fclose(fid2);
% 通过一系列的移位（bitshift）、位与（bitand）运算，将信号由二值数据转换为十进制数
M2H= bitshift(A(:,2), -4);        %字节向右移四位，即取字节的高四位
M1H= bitand(A(:,2), 15);          %取字节的低四位
PRL=bitshift(bitand(A(:,2),8),9);     % sign-bit   取出字节低四位中最高位，向右移九位
PRR=bitshift(bitand(A(:,2),128),5);   % sign-bit   取出字节高四位中最高位，向右移五位
M( : , 1)= bitshift(M1H,8)+ A(:,1)-PRL;
M( : , 2)= bitshift(M2H,8)+ A(:,3)-PRR;
if M(1,:) ~= firstvalue, error('inconsistency in the first bit values'); end;
switch nosig
case 2
    M( : , 1)= (M( : , 1)- zerovalue(1))/gain(1);
    M( : , 2)= (M( : , 2)- zerovalue(2))/gain(2);
    TIME=(0:(SAMPLES2READ-1))/sfreq;
case 1
    M( : , 1)= (M( : , 1)- zerovalue(1));
    M( : , 2)= (M( : , 2)- zerovalue(1));
    M=M';
    M(1)=[];
    sM=size(M);
    sM=sM(2)+1;
    M(sM)=0;
    M=M';
    M=M/gain(1);
    TIME=(0:2*(SAMPLES2READ)-1)/sfreq;
otherwise  % this case did not appear up to now!
    % here M has to be sorted!!!
    disp('Sorting algorithm for more than 2 signals not programmed yet!');
end;
clear A M1H M2H PRR PRL;
fprintf(1,'\\n$> LOADING DATA FINISHED \n');

%------ LOAD ATTRIBUTES DATA ----------------------------------------------
atrd= fullfile(PATH, ATRFILE);      % attribute file with annotation data
fid3=fopen(atrd,'r');
A= fread(fid3, [2, inf], 'uint8')';
fclose(fid3);
ATRTIME=[];
ANNOT=[];
sa=size(A);
saa=sa(1);
i=1;
while i<=saa
    annoth=bitshift(A(i,2),-2);
    if annoth==59
        ANNOT=[ANNOT;bitshift(A(i+3,2),-2)];
        ATRTIME=[ATRTIME;A(i+2,1)+bitshift(A(i+2,2),8)+...
                bitshift(A(i+1,1),16)+bitshift(A(i+1,2),24)];
        i=i+3;
    elseif annoth==60
        % nothing to do!
    elseif annoth==61
        % nothing to do!
    elseif annoth==62
        % nothing to do!
    elseif annoth==63
        hilfe=bitshift(bitand(A(i,2),3),8)+A(i,1);
        hilfe=hilfe+mod(hilfe,2);
        i=i+hilfe/2;
    else
        ATRTIME=[ATRTIME;bitshift(bitand(A(i,2),3),8)+A(i,1)];
        ANNOT=[ANNOT;bitshift(A(i,2),-2)];
   end;
   i=i+1;
end;
ANNOT(length(ANNOT))=[];       % last line = EOF (=0)
ATRTIME(length(ATRTIME))=[];   % last line = EOF
clear A;
ATRTIME= (cumsum(ATRTIME))/sfreq;
ind= find(ATRTIME <= TIME(end));
ATRTIMED= ATRTIME(ind);
ANNOT=round(ANNOT);
ANNOTD= ANNOT(ind);

%------ DISPLAY DATA ------------------------------------------------------
%{
figure(1); clf, box on, hold on
plot(TIME, M(:,1),'r');
if nosig==2
    plot(TIME, M(:,2),'b');
end;
for k=1:length(ATRTIMED)
    text(ATRTIMED(k),0,num2str(ANNOTD(k)));
end;
xlim([TIME(1), TIME(end)]);
xlabel('Time / s'); ylabel('Voltage / mV');
string=['ECG signal ',DATAFILE];
title(string);
fprintf(1,'\\n$> DISPLAYING DATA FINISHED \n');
%}
% -------------------------------------------------------------------------
fprintf(1,'\\n$> ALL FINISHED \n');

s=M(:,1);
s=s';
%%
[QRS_amp,QRS_ind] = DS_detect(s,0);% 调用QRS检测算法；
%%
Nt=size(QRS_ind,2);
R_TIME=ATRTIMED(ANNOTD==1 | ANNOTD ==2 | ANNOTD==3 | ANNOTD==5 | ANNOTD==8 |ANNOTD==9 );

REF_ind=round(R_TIME'.*360);
Nr=size(REF_ind,2);
ann=ANNOTD(ANNOTD==1 | ANNOTD ==2 | ANNOTD==3 | ANNOTD==5 | ANNOTD==8 |ANNOTD==9);

if Nt>Nr
    typ=0;
else
    typ=1;
end

if typ==0
    for n=1:Nr
        ref=REF_ind(n);
        for m=1:Nt
            act_ind=QRS_ind(m);
            if abs(ref-act_ind)<=54
               if act_ind<Pl || (SAMPLES2READ-act_ind)<Pr
                    break;
               else
                    SEG=s((act_ind-Pl+1):(act_ind+Pr));
               end
               switch ann(n)
                    case 1
                        Nb=[Nb;SEG];
                    case 2
                        Lb=[Lb;SEG];
                    case 3
                        Rb=[Rb;SEG];
                    case 5
                        Vb=[Vb;SEG];
                    case 8
                        Ab=[Ab;SEG];
                    case 9
                        Sb=[Sb;SEG];
               end
               
               break; 
            end 
        end
       
    end
else
    for n=1:Nt
        act_ind=QRS_ind(n);
        for m=1:Nr
            if abs(act_ind-REF_ind(m))<=54
                if act_ind<Pl || (SAMPLES2READ-act_ind)<Pr
                    break;
               else
                    SEG=s((act_ind-Pl+1):(act_ind+Pr));
               end
               switch ann(m)
                    case 1
                        Nb=[Nb;SEG];
                    case 2
                        Lb=[Lb;SEG];
                    case 3
                        Rb=[Rb;SEG];
                    case 5
                        Vb=[Vb;SEG];
                    case 8
                        Ab=[Ab;SEG];
                    case 9
                        Sb=[Sb;SEG];
               end
               
                break;
            end
        end
        
    end
   
end


end