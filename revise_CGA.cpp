#include <stdio.h>     //定义输入／输出函数
#include <stdlib.h>    //定义杂项函数及内存分配函数
#include <ctime>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//目标函数
double bar_func(double x){
    return  pow(x-0.5,2)+abs(x+3);
}

//得到编码长度
int get_m(int a,int b,int u){
    int m = 1;
    int temp = (b-a)*pow(10,u);
    while(1){
        if(temp<=pow(2,m)-1&&temp>pow(2,m-1)-1){
            return m;
            break;
        }
        else
            m = m+1;
    }
}

//适应度函数
double Fit(double x){
    int M = -3;
    return 1/(bar_func(x)+M);
}

//对种群进行初始化
void init(string B[],int N,int m){
    srand((unsigned)time(NULL) );
    for(int i=0;i<N;i++){
        B[i] = "";
        for(int j=0;j<m;j++){
            if(rand()%2==1) B[i]+="1";
            else B[i]+="0";
        }
    }
    return;
}

//double x = a+ (b-a)*x1/(2^m -1)
double get_x(string B[],int pos,int a,int b,int m){
    double res = 0,x1=0;
    for(int i=0;i<m;i++){
        x1= x1*2 + B[pos][i]-'0';
    }
    res = a + (b-a)*x1/(pow(2,m)-1);
    return res;
}

double get_x_best(string bestRes,int a,int b,int m){
    double res = 0,x1=0;
    for(int i=0;i<m;i++){
        x1= x1*2 + bestRes[i]-'0';
    }
    res = a + (b-a)*x1/(pow(2,m)-1);
    return res;
}

int main(){
    //0.参数
    int a = -1;       //自变量取值区间下界
    int b = 1;        //自变量取值区间上界
    int u = 5;        //精度要求
    int m ;           //编码长度
    double pc = 0.2;  //交叉的概率
    double pm = 0.001;//变异的概率
    int l = 5;        //交叉长度 m-l
    int t = 400;      //迭代次数
    int N = 40;       //种群规模
    string B[N];      //初始种群
    string bestRes;   //超级个体

    m = get_m(a,b,u); //得到编码长度
    init(B,N,m);      //对种群进行初始化
    bestRes = B[0];   //初始化超级个体
    srand((unsigned)time(NULL) );

    //每一代的操作都一样
    for(int o=0;o<t;o++){

        //1.选择进化的父母
        double sum_Fit =  0;
        double pi[N];      //每代被选中的概率
        double qi[N];      //计算累计概率
        string Next_B[N];  //用于进化的父母
        for(int i=0;i<N;i++) sum_Fit+= Fit(get_x(B,i,a,b,m) ) ;
        for(int i=0;i<N;i++) pi[i] = Fit(get_x(B,i,a,b,m))/sum_Fit;
        for(int i=0;i<N;i++){
            for(int j=0;j<=i;j++){
                qi[i] +=pi[j] ;
            }
        }
        for(int i=0;i<N;i++){
            double rand_temp0 = rand()/double(RAND_MAX);  //取得0-1之间的浮点数
            if(rand_temp0>=0&&rand_temp0<=qi[0]){
                    Next_B[i] = B[0];
                    continue;
            }
            for(int j=1;j<N;j++){
                if(rand_temp0>qi[j-1]&&rand_temp0<=qi[j] ){
                    Next_B[i] = B[j];
                    continue;
                }
            }
        }


        //2.交叉计算
        int num1 = 0;        //参与交叉的个体数量
        int num2 = 0;        //未参与交叉的个体的数量
        string temp1_B[N];   //参与交叉的个体
        string temp2_B[N];   //未参与交叉的个体
        string recombin_B[N];//交叉过的个体
        for(int i=0;i<N;i++){
            double rand_temp1 = rand()/double(RAND_MAX);  //取得0-1之间的浮点数
            if(rand_temp1<=pc ){
                temp1_B[num1++] = Next_B[i];
            }
            else
                temp2_B[num2++] = Next_B[i];
        }
        //如果上面选出的参与交叉的个体的数量为奇数，那么就去掉一个
        if(num1%2==1) {
            num1--;
            temp2_B[num2++] = temp1_B[num1];
        }
        //对于没有参与交叉的点直接看成本身的后代
        for(int i=0;i<num2;i++){
            recombin_B[i] = temp2_B[i];
        }
        //进行交叉操作
        for(int i=0;i<num1;i=i+2){
            string temp = temp1_B[i].substr(l,m-l);          //交叉起始点l在前面已经定义
            temp1_B[i] = temp1_B[i].substr(0,l)+temp1_B[i+1].substr(l,m-l);
            temp1_B[i+1] = temp1_B[i+1].substr(0,l)+temp;
        }
        //将交叉后的个体存入recombin_B
        for(int i=num2,j=0;i<N;i++,j++){
            recombin_B[i] = temp1_B[j];
        }


        //3.变异计算
        string varia_B[N];     //变异个体
        int varia_temp[N][m] ; //string不好操作，换成数组
        for(int i=0;i<N;i++){
            for(int j=0;j<m;j++){
                 varia_temp[i][j] = 0;
            }
        }
        for(int i=0;i<N;i++){
            for(int j=0;j<m;j++){
                 double rand_temp2 = rand()/double(RAND_MAX);  //取得0-1之间的浮点数
                 if(rand_temp2<=pm){
                    if((recombin_B[i][j]-'0')==1) varia_temp[i][j]==0;
                    else varia_temp[i][j]==1;
                 }
                 else varia_temp[i][j] =recombin_B[i][j]-'0';
            }
        }
        for(int i=0;i<N;i++){
            varia_B[i] = "";
            for(int j=0;j<m;j++){
                if(varia_temp[i][j]==1) varia_B[i]+="1";
                else varia_B[i]+="0";
            }
        }


        //4.选择下一代
        double sum_Fit_next =  0;
        double pi_next[N];      //每代被选中的概率
        double qi_next[N] ;     //计算累计概率
        for(int i=0;i<N;i++) sum_Fit_next+= Fit(get_x(varia_B,i,a,b,m) ) ;
        for(int i=0;i<N;i++) pi_next[i] = Fit(get_x(varia_B,i,a,b,m))/sum_Fit_next;
        for(int i=0;i<N;i++){
            for(int j=0;j<=i;j++){
                qi_next[i] +=pi_next[j] ;
            }
        }
        for(int i=0;i<N;i++){
            double rand_temp0 = rand()/double(RAND_MAX);  //取得0-1之间的浮点数
            if(rand_temp0>=0&&rand_temp0<=qi_next[0]){
                    B[i] = varia_B[0];
                    continue;
            }
            for(int j=1;j<N;j++){
                if(rand_temp0>qi_next[j-1]&&rand_temp0<=qi_next[j] ){
                    B[i] = varia_B[j];
                    continue;
                }
            }
        }

        //保留最好的个体作为超级个体
        for(int i=0;i<N;i++){
            double temp_best = Fit( get_x_best( bestRes,a,b,m));
            if(Fit( get_x(B,i,a,b,m) )> temp_best)  {
                bestRes = B[i];
                temp_best = Fit( get_x_best( bestRes,a,b,m));
            }
        }


    }

    //输出结果
    double x_res = get_x_best( bestRes,a,b,m);
    cout<<x_res<<endl<<bar_func(x_res)<<endl;
    system("pause");
    return 0;
}
