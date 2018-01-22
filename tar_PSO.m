%PSO的参数
clc;
close all;
clear;
c1 = 1.4944;     %c1  c2都是学习因子
c2 = 1.4944;
W_max = 0.9;
W_min = 0.4;
M = 100;         %PSO粒子进化代数  
num = 20;        %粒子数，好的粒子有10个，中等粒子有20个，差的粒子有10个
up_a = 20;       %a的上界，a的取值范围是[1,up_a]之间的随机整数
up_f_sum = 40;   %f1+f2的上界，f1\f2的取值范围是[1,up_f_sum-10]之间的随机整数

v = zeros(num,3);
x = zeros(num,3);
lbest = zeros(num,3);
fit = zeros(num,1);

lbest_x = zeros(M,3);
gbest_x = zeros(M,3);
lbest_x_fit = zeros(M,1);
gbest_x_fit = zeros(M,1);
gbest_poi = [0,0,0];

%生成初始种群
for i = 1:num
   temp_v = new_v();
   v(i,1) = temp_v(1);    v(i,2) = temp_v(2);   v(i,3) = temp_v(3);  
   temp_x = new_x(up_a,up_f_sum);
   x(i,1) = temp_x(1);    x(i,2) = temp_x(2);   x(i,3) = temp_x(3);
   lbest(i,1) = temp_x(1);    lbest(i,2) = temp_x(2);   lbest(i,3) = temp_x(3); %将当前粒子的最优值初始化为粒子本身
   
   fit(i) = fitness(temp_x);   %计算粒子对应的适应度  
end
[fmin,I] = min(fit);           % 找到fitness中的最小值
lbest_x(1,1) = x(I,1); lbest_x(1,2) = x(I,2);   lbest_x(1,3) = x(I,3);         % 初始化“每代最优值对应的个体”
lbest_x_fit(1)=  10^30;        % 初始化“每代最优值”   

gbest_poi(1) = x(I,1); gbest_poi(2) = x(I,2);gbest_poi(3) = x(I,3);        
gbest_fit =  10^30;         
gbest_x(1,1) = gbest_poi(1); gbest_x(1,2) = gbest_poi(2);gbest_x(1,3) = gbest_poi(3);    % 初始化“全局最优值对应的粒子个体”
gbest_x_fit(1) = gbest_fit;    % 初始化“全局最优值”
                       
%主循环，由于fitness计算的时候是使用了a,f1,f2三个参数，所以进化使用三个参数同时进化
for t = 1:M    
    
    %排序
    for i = 1:num-1
        for j =i+1:num
            if(fit(j)<fit(i))   %如果粒子i的适应度小于当前值，那就交换
                 temp1 = fit(i);  fit(i) = fit(j);  fit(j) = temp1; 
                 temp2_1 = x(i,1); x(i,1) = x(j,1); x(j,1) = temp2_1;
                 temp2_2 = x(i,2); x(i,2) = x(j,2); x(j,2) = temp2_2;
                 temp2_3 = x(i,3); x(i,3) = x(j,3); x(j,3) = temp2_3;
                    
                 temp3_1 = v(i,1); v(i,1) = v(j,1); v(j,1) = temp3_1;
                 temp3_2 = v(i,2); v(i,2) = v(j,2); v(j,2) = temp3_2;
                 temp3_3 = v(i,3); v(i,3) = v(j,3); v(j,3) = temp3_3;
                            
                 temp4_1 = lbest(i,1); lbest(i,1) = lbest(j,1); lbest(j,1) = temp4_1;
                 temp4_2 = lbest(i,2); lbest(i,2) = lbest(j,2); lbest(j,2) = temp4_2;
                 temp4_3 = lbest(i,3); lbest(i,3) = lbest(j,3); lbest(j,3) = temp4_3;                                         
            end            
        end        
    end
    
    %更新每个粒子的“历史最优值”
    for i = 1:num
        temp_lbest(1) = lbest(i,1); 
        temp_lbest(2) = lbest(i,2);
        temp_lbest(3) = lbest(i,3);
        
        if(fitness( temp_lbest) < fit(i)  )   %如果当前记录的最优个体没有“当前个体”好，那就更新
            lbest(i,:) = x(i,:); 
        end
    end
    
    %记录本代的最优值,由于前面已经排过序，所以这里直接取第一个就好。
    lbest_x(t,1) = x(1,1);  lbest_x(t,2) = x(1,2); lbest_x(t,3) = x(1,3);   % 记录“本代最优值对应的粒子个体”
    lbest_x_fit(t) =   fit(1);         % 记录“本代最优值”
    
    %更新全局最优值,由于前面已经排过序，所以这里直接取第一个就好。
    if(lbest_x_fit(t)<gbest_fit )
        gbest_fit = lbest_x_fit(t);   %更新全局最优值
        gbest_poi(1) = lbest_x(t,1); 
        gbest_poi(2) = lbest_x(t,2);
        gbest_poi(3) = lbest_x(t,3);  %更新全局最优值对应的粒子
    end 
    gbest_x_fit(t) = gbest_fit;       % 初始化“全局最优值”
    gbest_x(t,1) = gbest_poi(1); 
    gbest_x(t,2) = gbest_poi(2);
    gbest_x(t,3) = gbest_poi(3);      % 初始化“全局最优值对应的粒子个体”
    
    
    %好、中、差的粒子的数量；
    num_good = num * 0.25;
    num_mid = num * 0.5;
    num_bad = num * 0.25;
        
    
    %根据排序结果进化
    %好的粒子    
    for i1 = 1:num_good
        r1_1 = rand;                          %r1_1是[0,1]之间的随机数
        W_better = 0.5 + sqrt(0.3)*randn ;    %以0.5为均值，0.3为方差的随机正态分布    
        v(i1,1) = round( W_better*v(i1,1) + c1*r1_1*( lbest(i1,1)  - x(i1,1) )   ); %四舍五入取整 
        v(i1,2) = round( W_better*v(i1,2) + c1*r1_1*( lbest(i1,2)  - x(i1,2) )   ); %四舍五入取整 
        v(i1,3) = round( W_better*v(i1,3) + c1*r1_1*( lbest(i1,3)  - x(i1,3) )   ); %四舍五入取整 
        x(i1,1) = x(i1,1) + v(i1,1);
        x(i1,2) = x(i1,2) + v(i1,2);
        x(i1,3) = x(i1,3) + v(i1,3);
      %  legal(x(i1,1),x(i1,2),x(i1,3));  %检查是否越界，如果越界，将越界部分设为边界值
        
        if(x(i1,1)<1)  
            x(i1,1) = 1;  end
        if(x(i1,1)>up_a)  
            x(i1,1) = up_a;  end
        if(rand<=0.5)
            if(x(i1,2)<10) 
                x(i1,2) = 10; end
            if(x(i1,2)>up_f_sum-10) 
                x(i1,2) = up_f_sum-10; end 
            if(x(i1,3)<10) 
                x(i1,3) = 10; end 
            if(x(i1,3)>up_f_sum-x(i1,2)) 
                x(i1,3) = up_f_sum-x(i1,2);end %f2由f1来约束
        else            
            if(x(i1,3)<10) 
                x(i1,3) = 10; end
            if(x(i1,3)>up_f_sum-10) 
                x(i1,3) = up_f_sum-10; end 
            if(x(i1,2)<10) 
                x(i1,2) = 10; end 
            if(x(i1,2)>up_f_sum-x(i1,3)) 
                x(i1,2) = up_f_sum-x(i1,3);end %f1由f2来约束             
        end

        %如果不满足最后一个约束，就让f1\f2依次减1，直到满足约束
        [T1,T2] = computer_T1_T2(x(i1,1));    
        count = round(rand*100);
        while(T1*x(i1,2)/3600 + T2*x(i1,3)/3600 > 41)          
            if(mod(count,2)==1 )
                if(x(i1,2)>1)  
                    x(i1,2) = x(i1,2) - 1; end 
                count = count + 1;
            else
                if(x(i1,3)>1) 
                    x(i1,3) = x(i1,3) - 1; end 
                count = count + 1;
            end        
        end    
        
        
    end
    
    %中等粒子    
    %y_meand是第d维所有粒子的历史最优的平均值
    temp_sum_a = 0;
    temp_sum_f1 = 0;
    temp_sum_f2 = 0;
    for ii = 1:num
        temp_sum_a =  temp_sum_a + lbest(ii,1);
        temp_sum_f1 = temp_sum_f1 + lbest(ii,2);
        temp_sum_f2 = temp_sum_f2 + lbest(ii,3);
    end        
    y_meand.a = temp_sum_a / num;
    y_meand.f1 = temp_sum_f1 / num;
    y_meand.f2 = temp_sum_f2 / num;
    
    w = W_max - t/M*(W_max-W_min);        %惯性权重
    for i2 = num_good+1:num_good+num_mid
        r1_2 = rand;                          %r1  r2是[0,1]之间的随机数
        r2_2 = rand;
        v(i2,1) = round( w*v(i2,1) + c1*r1_2*( lbest(i2,1) - x(i2,1) ) + c2*r2_2*(y_meand.a -x(i2,1) )) ;
        v(i2,2) = round( w*v(i2,2) + c1*r1_2*( lbest(i2,2) - x(i2,2) ) + c2*r2_2*(y_meand.f1 -x(i2,2) )) ;
        v(i2,3) = round( w*v(i2,3) + c1*r1_2*( lbest(i2,3) - x(i2,3) ) + c2*r2_2*(y_meand.f2 -x(i2,3) )) ;
        x(i2,1) = x(i2,1) + v(i2,1);
        x(i2,2) = x(i2,2) + v(i2,2);
        x(i2,3) = x(i2,3) + v(i2,3);
    % legal(x(i2,1),x(i2,2),x(i2,3)); 
        if(x(i2,1)<1)  
            x(i2,1) = 1;  end
        if(x(i2,1)>up_a)  
            x(i2,1) = up_a;  end
        if(rand<=0.5)
            if(x(i2,2)<10) 
                x(i2,2) = 10; end
            if(x(i2,2)>up_f_sum-10) 
                x(i2,2) = up_f_sum-10; end 
            if(x(i2,3)<10) 
                x(i2,3) = 10; end 
            if(x(i2,3)>up_f_sum-x(i2,2)) 
                x(i2,3) = up_f_sum-x(i2,2);end %f2由f1来约束
        else            
            if(x(i2,3)<10) 
                x(i2,3) = 10; end
            if(x(i2,3)>up_f_sum-10) 
                x(i2,3) = up_f_sum-10; end 
            if(x(i2,2)<10) 
                x(i2,2) = 10; end 
            if(x(i2,2)>up_f_sum - x(i2,3)) 
                x(i2,2) = up_f_sum - x(i2,3);end %f2由f1来约束            
        end
        
        %如果不满足最后一个约束，就让f1\f2依次减1，直到满足约束
        [T1,T2] = computer_T1_T2(x(i2,1));    
        count = round(rand*100);
        while(T1*x(i2,2)/3600 + T2*x(i2,3)/3600 > 41)          
            if(mod(count,2)==1 )
                if(x(i2,2)>1)  
                    x(i2,2) = x(i2,2) - 1; end 
                count = count + 1;
            else
                if(x(i2,3)>1) 
                    x(i2,3) = x(i2,3) - 1; end 
                count = count + 1;
            end        
        end   
    
    
    end
    
    %差的粒子
    for i3 = num_good+num_mid+1:num
        r1_3 = rand;                          %r1_3是[0,1]之间的随机数
        poi = round(rand * (num_good-1)) + 1; %产生[1,num_good]之间的随机整数
        v(i3,1) = round (c1*r1_3* (x(poi,1) - x(i3,1) ));
        v(i3,2) = round (c1*r1_3* (x(poi,2) - x(i3,2) )) ;
        v(i3,3) = round (c1*r1_3* (x(poi,3) - x(i3,3) )) ;
        x(i3,1) = x(i3,1) + v(i3,1);
        x(i3,2) = x(i3,2) + v(i3,2);
        x(i3,3) = x(i3,3) + v(i3,3);
%        legal(x(i3,1),x(i3,2),x(i3,3));            
        if(x(i3,1)<1)  
            x(i3,1) = 1;  end
        if(x(i3,1)>up_a)  
            x(i3,1) = up_a;  end
        
        if(rand<=0.5)            
            if(x(i3,2)<10) 
                x(i3,2) = 10; end
            if(x(i3,2)>up_f_sum-10) 
                x(i3,2) = up_f_sum-10; end 
            if(x(i3,3)<10) 
                x(i3,3) = 10; end 
            if(x(i3,3)>up_f_sum-x(i3,2)) 
                x(i3,3) = up_f_sum-x(i3,2);end %f2由f1来约束
        else
            if(x(i3,3)<10) 
                x(i3,3) = 10; end
            if(x(i3,3)>up_f_sum-10) 
                x(i3,3) = up_f_sum-10; end 
            if(x(i3,2)<10) 
                x(i3,2) = 10; end 
            if(x(i3,2)>up_f_sum-x(i3,3)) 
                x(i3,2) = up_f_sum-x(i3,3);end %f2由f1来约束
        end        

        %如果不满足最后一个约束，就让f1\f2依次减1，直到满足约束
        [T1,T2] = computer_T1_T2(x(i3,1));    
        count = round(rand*100);
        while(T1*x(i3,2)/3600 + T2*x(i3,3)/3600 > 41)          
            if(mod(count,2)==1 )
                if(x(i3,2)>1)  
                    x(i3,2) = x(i3,2) - 1; end 
                count = count + 1;
            else
                if(x(i3,3)>1) 
                    x(i3,3) = x(i3,3) - 1; end 
                count = count + 1;
            end        
        end           
    end
    
    %更新进化后的粒子的适应度
    for i4 = 1:num
        temp_xi4(1) = x(i4,1); temp_xi4(2) = x(i4,2); temp_xi4(3) = x(i4,3);
        fit(i4) = fitness(temp_xi4 ); 
    end
    
end

draw(lbest_x_fit,gbest_x_fit,gbest_x,M,up_a,up_f_sum);   %画图
disp('the best result is:');               %输出最优值
disp(gbest_x_fit(M));
disp('The corresponding solution are:');
disp(gbest_x(M,1)); 
disp(gbest_x(M,2));
disp(gbest_x(M,3));
