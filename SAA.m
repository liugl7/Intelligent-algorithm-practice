clear ;
close all;
format long ;

%冷却表参数
times = 100;                      % 大循环的执行次数
MarkovLength = 20;          % 马可夫链长度
DecayScale = 0.95;            % 衰减参数 
da = 0.95;                         %第二种退火公式的参数
Temperature0 = 100;               % 初始温度
Temperature = Temperature0;
Tolerance = 1e-5;              % 容差
AcceptPoints = 0.0;           % Metropolis过程中总接受点

% 初始化矩阵
A= [  0.2,  0.2,  0.2,  0.2,  0.22;
      0.2,  0.2,  0.2,  0.2,  0.22;
      0.2,  0.2,  0.2,  0.2,  0.22;
      0.2,  0.2,  0.2,  0.2,  0.22;
      0.23,  0.23,  0.23,  0.23,  0.46];
[M,N] = size(A);
Sol = zeros(1,M);    %定义解向量

%随机选点.为迭代设立M个初始值，每个初始值的范围是[1,N]之间的任意整数    
SF=randperm(N);     %将1~N 随机打乱
INDEX=randperm(M);  %将1~M 随机打乱
for k=1:M
    Sol(k)=SF(INDEX(k)); 
end   
Fitness = ObjectFunction(Sol);
Best = Sol;             % 上一个最优解
fitness_best = Fitness;  % 在某一个温度下的最优解对应的“匹配度的值” 
history_fit_best = zeros(1,times);
t=0;

%每迭代一次退火一次(降温), 直到满足迭代条件为止 
while(t<times)
    t = t + 1;  %次数
    %Temperature = Temperature *DecayScale;  %这是原来的写法，也就是博客的第三种方法 ，相当于a设置为0.95
    Temperature = Temperature0/log(1+da*t);    %这是博客中的第二种写法，其中参数a设置为0.95
    
    AcceptPoints = 0.0;
    %在当前温度T下迭代loop(即MARKOV链长度)次
    for i= 1:1:MarkovLength
      %(1)  在当前解附近随机选取下一个解
  		% 变异原解，得到新解
        for r=1:1:M
            temp1 = rand;
            if (temp1<0.333) 
                ds = -1;
            elseif(temp1<0.666&&temp1>=0.333) 
                ds = 0;
            else
                ds = 1;
            end            
            Sol(r) = Sol(r) + ds;                               
        end

        % 解决越界的元素
        R1=setdiff(randperm(N),Sol);           %求补集
        i1=1;
  	    for k=1:M
  	        if(Sol(k)==0||Sol(k)==N+1)   %如果取的解有0或者N+1，那就用补集的元素替换
  	          Sol(k)=R1(i1);
  	          i1=i1+1;
  	        end
  	    end

        % 解决重复的元素
        [B, I] = unique(Sol, 'first');
        R2 = setdiff(1:numel(Sol), I); %求Sol的“重复集”
        temp_R2 = setdiff(randperm(N),unique(Sol) );  %求Sol的“标准补”
        for r = 1:1:length(R2)
            Sol(R2(r)) = temp_R2(r);  % 将重复元素用“标准补”的元素替换掉               
        end
        
        %计算当前变异解(新解)的适应度
        Fitness = ObjectFunction(Sol);
        
        %(2) 更新全局的最优解
        if (fitness_best > Fitness)
             Best = Sol; 
             fitness_best = Fitness;
        end 
        
        % (3) Metropolis过程 
        if( fitness_best <  Fitness  )
            % 接受, 此处lastPoint即下一个迭代的点以新接受的点开始
            Sol = Best;   %如果(2)中Best没有被更新，那么还原"新解"不变 
            AcceptPoints = AcceptPoints + 1;  
        else
            change = -1 * ( Fitness - fitness_best ) / Temperature ;
            if(exp(change)>rand)
                Best = Sol;
                AcceptPoints = AcceptPoints + 1;  
            end                 
        end
        
    end
    
    history_fit_best(t) = fitness_best; 
    
    disp('-----------------------------');
    disp('温度：');
    disp(Temperature);   %输出当前“温度”
    disp('最优解的函数值为：')
    disp(fitness_best);  %输出当前最优解对应的函数值
end

disp('-----------------------------');
disp('最优解为:');
disp(Best);
disp('最优解对应的函数值为:');
disp(ObjectFunction(Best));

plot(1:t,history_fit_best(1:t),'-.r') ; 
xlabel('迭代次数');                        %  坐标标注
ylabel('全局历史最佳适应值');
title('全局历史最佳适应值趋势图：');
legend('模拟退火算法');

function z = ObjectFunction(Sol)
% 初始化矩阵
   A = [  0.2,  0.2,  0.2,  0.2,  0.22;
          0.2,  0.2,  0.2,  0.2,  0.22;
          0.2,  0.2,  0.2,  0.2,  0.22;
          0.2,  0.2,  0.2,  0.2,  0.22;
          0.23,  0.23,  0.23,  0.23,  0.46];  
  [row,~] = size(A); 
  z = 0;
  for i=1:1:row
      z = z+ A(i,Sol(i));       
  end
  z = 1/z;
  
%{
A=  [0.8 , 0.1 , 0.1 , 0.2, 0.4, 0.2, 0.6, 0.7, 0.8, 0.1;
       0.3 , 0.9 , 0.8 , 0.2, 0.2, 0.3, 0.2, 0.4, 0.6, 0.7;  
       0.7, 0.6, 0.6 , 0.7, 0.2, 0.1, 0.5, 0.6, 0.8, 0.9;   
       0.4 , 0.6 , 0.1 , 0.2, 0.1, 0.4, 0.7, 0.9, 0.1, 0.2];
A=  [0.8 , 0.1 , 0.1 , 0.2, 0.4, 0.2;
      0.3 , 0.9 , 0.8 , 0.2, 0.2, 0.3;  
      0.7, 0.6, 0.6 , 0.7, 0.2, 0.1;   
      0.4 , 0.6 , 0.1 , 0.2, 0.1, 0.4];
%}
end
 