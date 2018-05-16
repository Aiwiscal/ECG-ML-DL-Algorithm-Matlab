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
%% 控制使用数据量，每一类5000，并生成标签,one-hot编码；
fprintf('Data preprocessing...\n');
tic;
Nb=Nb(1:5000,:);Label1=repmat([1;0;0;0],1,5000);
Vb=Vb(1:5000,:);Label2=repmat([0;1;0;0],1,5000);
Rb=Rb(1:5000,:);Label3=repmat([0;0;1;0],1,5000);
Lb=Lb(1:5000,:);Label4=repmat([0;0;0;1],1,5000);

Data=[Nb;Vb;Rb;Lb];
Label=[Label1,Label2,Label3,Label4];

clear Nb;clear Label1;
clear Rb;clear Label2;
clear Lb;clear Label3;
clear Vb;clear Label4;
Data=Data-repmat(mean(Data,2),1,250); %使信号的均值为0，去掉基线的影响；
fprintf('Finished!\n');
toc;
fprintf('=============================================================\n');

%% 数据划分与模型训练测试；
fprintf('Model training and testing...\n');
Nums=randperm(20000);      %随机打乱样本顺序，达到随机选择训练测试样本的目的；
train_x=Data(Nums(1:10000),:);
test_x=Data(Nums(10001:end),:);
train_y=Label(:,Nums(1:10000));
test_y=Label(:,Nums(10001:end));
train_x=train_x';
test_x=test_x';

cnn.layers = {
    struct('type', 'i') %input layer
    struct('type', 'c', 'outputmaps', 4, 'kernelsize', 31,'actv','relu') %convolution layer
    struct('type', 's', 'scale', 5,'pool','mean') %sub sampling layer
    struct('type', 'c', 'outputmaps', 8, 'kernelsize', 6,'actv','relu') %convolution layer
    struct('type', 's', 'scale', 3,'pool','mean') %subsampling layer
};
cnn.output = 'softmax';  %确定cnn结构；
                         %确定超参数；
opts.alpha = 0.01;       %学习率；
opts.batchsize = 16;     %batch块大小；
opts.numepochs = 30;     %迭代epoch；

cnn = cnnsetup1d(cnn, train_x, train_y);      %建立1D CNN;
cnn = cnntrain1d(cnn, train_x, train_y,opts); %训练1D CNN;
[er,bad,out] = cnntest1d(cnn, test_x, test_y);%测试1D CNN;

[~,ptest]=max(out,[],1);
[~,test_yt]=max(test_y,[],1);

Correct_Predict=zeros(1,4);                     %统计各类准确率；
Class_Num=zeros(1,4);                           %并得到混淆矩阵；
Conf_Mat=zeros(4);
for i=1:10000
    Class_Num(test_yt(i))=Class_Num(test_yt(i))+1;
    Conf_Mat(test_yt(i),ptest(i))=Conf_Mat(test_yt(i),ptest(i))+1;
    if ptest(i)==test_yt(i)
        Correct_Predict(test_yt(i))= Correct_Predict(test_yt(i))+1;
    end
end

ACCs=Correct_Predict./Class_Num;
fprintf('Accuracy = %.2f%%\n',(1-er)*100);
fprintf('Accuracy_N = %.2f%%\n',ACCs(1)*100);
fprintf('Accuracy_V = %.2f%%\n',ACCs(2)*100);
fprintf('Accuracy_R = %.2f%%\n',ACCs(3)*100);
fprintf('Accuracy_L = %.2f%%\n',ACCs(4)*100);
