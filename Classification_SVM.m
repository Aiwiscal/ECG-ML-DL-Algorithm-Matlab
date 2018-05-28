clear;clc;
%% 载入数据；
fprintf('Loading data...\n');
tic;
load('N_dat.mat');
load('L_dat.mat');
load('R_dat.mat');
load('V_dat.mat');
fprintf('Finished!\n');
toc;
fprintf('=============================================================\n');
%% 控制使用数据量，每一类5000，并生成标签；
fprintf('Data preprocessing...\n');
tic;
Nb=Nb(1:5000,:);Label1=ones(1,5000);%Label1=repmat([1;0;0;0],1,5000);
Vb=Vb(1:5000,:);Label2=ones(1,5000)*2;%Label2=repmat([0;1;0;0],1,5000);
Rb=Rb(1:5000,:);Label3=ones(1,5000)*3;%Label3=repmat([0;0;1;0],1,5000);
Lb=Lb(1:5000,:);Label4=ones(1,5000)*4;%Label4=repmat([0;0;0;1],1,5000);

Data=[Nb;Vb;Rb;Lb];
Label=[Label1,Label2,Label3,Label4];
Label=Label';

clear Nb;clear Label1;
clear Rb;clear Label2;
clear Lb;clear Label3;
clear Vb;clear Label4;
Data=Data-repmat(mean(Data,2),1,250); %使信号的均值为0，去掉基线的影响；
fprintf('Finished!\n');
toc;
fprintf('=============================================================\n');
%% 利用小波变换提取系数特征，并切分训练和测试集；
fprintf('Feature extracting and normalizing...\n');
tic;
Feature=[];
for i=1:size(Data,1)
    [C,L]=wavedec(Data(i,:),5,'db6');  %% db6小波5级分解；
    Feature=[Feature;C(1:25)];
end

Nums=randperm(20000);      %随机打乱样本顺序，达到随机选择训练测试样本的目的；
train_x=Feature(Nums(1:10000),:);
test_x=Feature(Nums(10001:end),:);
train_y=Label(Nums(1:10000));
test_y=Label(Nums(10001:end));

[train_x,ps]=mapminmax(train_x',0,1); %利用mapminmax内建函数特征归一化到0，1之间；
test_x=mapminmax('apply',test_x',ps);
train_x=train_x';test_x=test_x';
fprintf('Finished!\n');
toc;
fprintf('=============================================================\n');
%% 训练SVM，并测试效果；
fprintf('SVM training and testing...\n');
tic;
model=libsvmtrain(train_y,train_x,'-c 2 -g 1'); %模型训练；
[ptest,~,~]=libsvmpredict(test_y,test_x,model); %模型预测；

Correct_Predict=zeros(1,4);                     %统计各类准确率；
Class_Num=zeros(1,4);
Conf_Mat=zeros(4);
for i=1:10000
    Class_Num(test_y(i))=Class_Num(test_y(i))+1;
    Conf_Mat(test_y(i),ptest(i))=Conf_Mat(test_y(i),ptest(i))+1;
    if ptest(i)==test_y(i)
        Correct_Predict(test_y(i))= Correct_Predict(test_y(i))+1;
    end
end
ACCs=Correct_Predict./Class_Num;
fprintf('Accuracy_N = %.2f%%\n',ACCs(1)*100);
fprintf('Accuracy_V = %.2f%%\n',ACCs(2)*100);
fprintf('Accuracy_R = %.2f%%\n',ACCs(3)*100);
fprintf('Accuracy_L = %.2f%%\n',ACCs(4)*100);
toc;
