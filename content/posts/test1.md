+++
title = 'test1'
date = 2024-03-17T23:01:20+08:00
draft = false
tags = ['测试']
categories = ['测试']
+++


# 第四章 数学



## 1.质数

质数是指在大于1的自然数中，除了1和它本身以外不再有其他因数的自然数



### 1.1 试除法判定质数

问题描述：给定一个数，判断它是否是质数



枚举从2到 x-1 的自然数，如果有可以整除 x 的数，则 x 不是质数

实际上只需要枚举从2到根号 x 的数即可



循环的终止条件有三种写法

1.`i<=sqrt(x)` 循环每执行一次，就调用一次sqrt函数，而且要注意类型转换，最好提前存好`int t=sqrt(x)`

2.`i*i<=x` 缺点是容易爆`int`

3.`i<=x/i` 这种比较好



***代码：O(n^0.5^)***

```c++
bool check(int x){
    if(x<2) return false;
    for(int i=2;i<=x/i;i++)
        if(x%i==0) return false;
    return true;
}
```



### 1.2 分解质因数

问题描述：给定一个正整数，将这个数分解质因数，并按照质因数从小到大的顺序输出每个质因数的底数和指数。



每找到一个因子，就将它除尽，直到不含这个因子，还是枚举到根号 x 

注意：每个数最多只含有一个大于根号 x 的因子，最后要特殊判断



***代码：***

```c++
void divide(int x){
    for(int i=2;i<=x/i;i++){
        int cnt=0;
        while(x%i==0)
            cnt++, x/=i;
        if(cnt>0)
            printf("%d %d\n",i,cnt);
    }
    if(x>1) printf("%d %d\n",x,1);
}
```



### 1.3 筛法求质数

问题描述：给定一个正整数 n，请你求出 1∼n 中质数的个数。



2是质数，则4，6，8，10.....不是质数

3是质数，则 6，9，12，15.....不是质数

据此可以优化，不需要都用试除法判断



***代码（朴素筛）：O(nloglogn)***

```c++
#include<iostream>
using namespace std;

const int N=1000010;
int cnt,p[N/10];
bool vis[N];

void get_primes(int n){
    for(int i=2;i<=n;i++){
        if(vis[i]) continue;
        p[++cnt]=i;
        for(int j=i+i;j<=n;j+=i)
            vis[j]=true;
    }
}

int main(){
    
    int n;
    scanf("%d",&n);
    
    get_primes(n);
    
    printf("%d",cnt);
    
    return 0;
}
//O(n/2+n/3+n/5+n/7+……)=O(nloglogn)
```



筛掉合数时，如果保证只用它最小的质因子筛掉一次

可以保证时间复杂度是严格线性的



***代码（线性筛）：O(n)***

```c++
void get_primes(int n){
    for(int i=2;i<=n;i++){
        if(!vis[i]) p[++cnt]=i;
        for(int j=1;p[j]<=n/i;j++){
            vis[i*p[j]]=true;
            if(i%p[j]==0) break;    
        }
    }
}
```

____



## 2.约数

若整数 d 能整除整数 n，则 d 是 n 的约数，n 是 d 的倍数



### 2.1 试除法求约数

问题描述：给定 n 个正整数 **a~i~**，对于每个整数 **a~i~**，请你按照从小到大的顺序输出它的所有约数。



约数成对出现，还是只枚举到根号 a 即可，最后排序输出，注意完全平方数的特判



***代码：***

```c++
#include<iostream>
#include<algorithm>
using namespace std;

int ans[2000],cnt;

void get(int x){
    cnt=0;
    for(int i=1;i<=x/i;i++){
        if(x%i==0){
            ans[++cnt]=i;
            if(x/i!=i)    //完全平方数的因子只出现一次 
                ans[++cnt]=x/i;
        }
    }
    sort(ans+1,ans+cnt+1);
}

int main(){
    
    int n,x;
    scanf("%d",&n);
    
    while(n--){
        scanf("%d",&x);
        get(x);
        for(int i=1;i<=cnt;i++) printf("%d ",ans[i]);
        printf("\n");
    }
    
    return 0;
}
```



###2.2 约数个数

问题描述：给定 n 个正整数 **a~i~**，请你输出这些数的乘积的约数个数，答案对 1e9+7 取模。

  

将这个乘积分解质因数，设分解后的结果为
$$
p_1^{a_1}*p_2^{a_2}*p_3^{a_3}*……*p_m^{a_m}
$$
由乘法原理，约数的个数为
$$
(a_1+1)(a_2+1)(a_3+1)……(a_m+1)
$$
用哈希表存质因子



***代码：***

```c++
#include<iostream>
#include<unordered_map>
using namespace std;

const int mod=1e9+7;

unordered_map<int,int> h;

void get(int x){
    for(int i=2;i<=x/i;i++){
        while(x%i==0)
            h[i]++, x/=i;
    }
    if(x>1) h[x]++;
}

int main()
{
    int n,x;
    scanf("%d",&n);
    
    for(int i=1;i<=n;i++){
        scanf("%d",&x);
        get(x);
    }
    
    int ans=1;
    for(pair<int,int> p : h)
        ans=1LL*ans*(p.second+1)%mod;
     
    printf("%d",ans);
    
    return 0;
}
```



### 2.3 约数之和

问题描述：给定 n 个正整数 ai，请你输出这些数的乘积的约数之和，答案对 1e9+7 取模。



约数的和为
$$
(1+p_1+p_1^2+……+p_1^{a_1})(1+p_2+p_2^2+……+p_2^{a_2})……(1+p_m+p_m^2+……+p_m^{a_m})
$$


***代码：***

```c++
#include<iostream>
#include<unordered_map>
using namespace std;

const int mod=1e9+7;

unordered_map<int,int> h;

void get(int x){
    for(int i=2;i<=x/i;i++)
        while(x%i==0)
            h[i]++, x/=i;
    if(x>1) h[x]++;
}

int main(){
    
    int n,x;
    scanf("%d",&n);
    
    while(n--){
        scanf("%d",&x);
        get(x);
    }
    
    int ans=1;
    for(pair<int,int> p : h) {
        int a=p.first,b=p.second;
        int t=1;
        while(b--) t=(1LL*t*a+1)%mod;    //求1+p1+p1^2+……p1^a1
        ans=1LL*ans*t%mod;
    }
    
    printf("%d",ans);
    
    return 0;
}
```



### 2.4 欧几里得算法

问题描述：**最大公约数**

给定一对正整数 a，b，求出它们的最大公约数



下面证明**`a，b`，与`b，a%b`有相同的公因数集**

$$
假设 d 是 a，b 的公因数，\\

因为  a\%b 可以表示成 a-kb，因为 d 整除 b，a-kb 也可以表示成 a-kd  \\

又因为 d 整除 a ，a-kd 就是 d 的倍数，即 d 是b，a\%b的公因数\\

同理假设 d 是 b，a\%b 的公因数 \\

由 a-kb 是 d 的倍数，b 是 d 的倍数，所以 a 是 d 的倍数\\

所以  d 也是 a，b 的公因数\\
$$


由此可以得到欧几里得算法，非常简洁



***代码：***

```c++
int gcd(int a,int b){
    return b?gcd(b,a%b):a;
}
```

_____



问题描述：**最小公倍数**

***代码：***

```c++
int lcm(int a,int b){
    return a/gcd(a,b)*b;
}
```

_____





## 3.欧拉函数

$$
1∼N 中与 N 互质的数的个数被称为欧拉函数，记为 ϕ(N)。
$$

$$
若在算数基本定理中，N=p_1^{a_1}p_2^{a_2}p_3^{a_3}……p_m^{a_m}，则：\\
ϕ(N) = N×\frac{p_1-1}{p_1}×\frac{p_2-1}{p_2}×……×\frac{p_m-1}{p_m}
$$



证明用的是容斥原理
$$
\begin{aligned}
ϕ(N) &= N-\frac{N}{p_1}-\frac{N}{p_2}-\frac{N}{p_3}……\\
&+\frac{N}{p_1p_2}+\frac{N}{p_1p_3}+\frac{N}{p_2p_3}+……\\
&-\frac{N}{p_1p_2p_3}-\frac{N}{p_1p_2p_4}-\frac{N}{p_1p_3p_4}-\frac{N}{p_2p_3p_4}……\\
&+\frac{N}{p_1p_2p_3p_4}……
\end{aligned}
$$
展开即为上式

 

###3.1 朴素求欧拉函数

问题描述：给定一个数 x ，求它的欧拉函数值



循环的主要方式仍然类似分解质因数

***代码：***

```c++
int phi(int x){
    int ans=x;
    for(int i=2;i<=x/i;i++){
        if(x%i==0) ans=ans/i*(i-1);    //先除再乘
        while(x%i==0) x/=i;
    }
    if(x>1) ans=ans/x*(x-1);
    return ans;
}
```



 ### 3.2 筛法求欧拉函数

在线性筛质数的过程中，可以顺带求出欧拉函数，适合处理多次询问



问题描述：给定一个正整数 n，求 1∼n 中每个数的欧拉函数之和。

***代码：***

```c++
#include<iostream>
using namespace std;

const int N=1000010;

int p[N],phi[N],cnt;
bool vis[N];

void get_eulers(int n){
    phi[1]=1;    //1的欧拉函数值就是1
    for(int i=2;i<=n;i++){
        if(!vis[i]) p[++cnt]=i, phi[i]=i-1;    //质数i的欧拉函数值等于i-1
        for(int j=1;p[j]<=n/i;j++){
            vis[i*p[j]]=true;
            if(i%p[j]==0) {
                phi[i*p[j]]=phi[i]*p[j];    //含有质因子p[j]
                break;
            }
            else phi[i*p[j]]=phi[i]*(p[j]-1);    //不含质因子p[j]
        }
    }
}

int main(){
    
    int n;    
    scanf("%d",&n);
    
    get_eulers(n);
    
    long long ans=0;
    for(int i=1;i<=n;i++) ans+=phi[i];
    
    printf("%lld",ans);
    
    return 0;
}
```



### 3.3 欧拉定理

$$
若正整数a,n互质，则\;
a^{ϕ(n)}\equiv1\;(mod\;n)
$$

证明：
$$
在1到n的正整数中，一定可以找出ϕ(n)个与n互质的数\{a_1,a_2……a_{ϕ(n)}\}\\
考虑另一组数\{aa_1\;mod\;n,aa_2\;mod\;n,……,aa_{ϕ(n)}\;mod\;n\}\\
由于aa_k和n，aa_k\;mod\;n和n有相同的公因数集合(1\leqslant k\leqslant ϕ(n))\\
aa_k和n互质，所以aa_k\;mod\;n和n也互质\\
同时从这组数中任取两个，一定没有相同的\\
因为如果相同，则存在(a_i-a_j)\equiv0\;(mod\;n)\;\;\;(1\leqslant i,j\leqslant ϕ(n))\\\
而0< a_i-a_j< n，(不妨a_i>a_j)\\
故\{aa_1\;mod\;n,aa_2\;mod\;n……aa_{ϕ(n)}\;mod\;n\}也是均小于n的，ϕ(n)个互不相同的和n互质的数\\
和\{a_1,a_2……a_{ϕ(n)}\}实际上是同一组数\\
故a^{ϕ(n)}a_1a_2……a_{ϕ(n)}\equiv a_1a_2……a_{ϕ(n)}\;(mod\;n)\\
又a_1a_2……a_{ϕ(n)}与互质\\
故a^{ϕ(n)}\equiv1\;(mod\;n)，得证
$$
欧拉定理揭示了**同余的周期性**



### 3.4 费马小定理

$$
若p是质数，正整数a不是p的整数倍，则\;
a^{p-1}\equiv1\;(mod\;n)
$$

$$
如果p是质数，则ϕ(p)=p-1\\
费马小定理实际上是欧拉定理的推论
$$

____



## 4.快速幂

问题描述：给定**a，b，p**，求 **a^b^ % p**



暴力循环时间复杂度是线性的

把指数 b 看成二进制数，预处理 a^1^，a^2^，a^4^……，其实就是不断平方 a

二进制的 b 某一位是 1 时，加上对应预处理的结果就可以了

事实上预处理无需存储，在预处理的过程中就可以求出结果了，时间复杂度是对数级别的



***代码：O(log n)***

```c++
int qpow(int a,int b,int p){
    int ans=1;
    while(b){
        if(b&1) ans=1LL*ans*a%p;
        b>>=1; a=1LL*a*a%p;
    }
    return ans;
}
```

______



## 5.逆元

 ### 5.1 扩展欧几里得算法

问题描述：**二元一次不定方程**
$$
给定一对正整数 a,b，求出一组整数 x,y，使其满足ax+by=gcd(a,b)。\\
可以证明一定存在这样一组x,y
$$
方法：
$$
如果b=0，那么x=1,y=0就是一组解\\
由gcd(a,b)=gcd(b,a\%b)，如果反复使用欧几里得算法处理a,b\\
最终一定可以使b=0，只需用x=1,y=0反推出原始答案即可\\
\\
问题转化为如果已知方程bx+(a\%b)y=gcd(b,a\%b)的解\\
如何求出ax+by=gcd(a,b)的解\\
\\ 
将方程变形为bx'+(a-a/b*b)y'=gcd(b,a\%b)\\
即ay'+b(x'-a/b*y')=gcd(b,a\%b)\\
对比ax+by=gcd(a,b)\\
有x=y',y=x'-a/b*y
$$
如果等号右边的数**不是gcd(a,b)**，是**任意整数**
$$
若方程为ax+by=n\\
此时有整数解的条件是gcd(a,b)|n\\
否则将方程化为gcd(a,b)(a'x+b'y)=n\\
就能发现要使(a'x+b'y)为整数，gcd(a,b)一定整除n\\
\\
先求ax'+by'=gcd(a,b)的解\\
则原方程的解x=x'n/gcd(a,b),y=y;n/gcd(a,b)
$$
实际写代码可以优化一下，每次递归交换x，y的位置
$$
将方程变形为by'+(a-a/b*b)x'=gcd(b,a\%b)\\
即ax'+b(y'-a/b*x')=gcd(b,a\%b)\\
对比ax+by=gcd(a,b)\\
有x=x',y=y'-a/b*x
$$
时间复杂度是对数级的

***代码：***

```c++
//求的是是方程ax+by=gcd(a,b)的解，注意等号右边！
//返回的是a,b的最大公约数，记为g
//有解的条件：n%g==0

//原版
int exgcd(int a,int b,int &x,int &y){
    if(!b){
        x=1; y=0;
        return;
    }
    int d=exgcd(b,a%b,x,y);
    int t=x;
    x=y;
    y=t-a/b*y;
    return d;
}

//交换版
int exgcd(int a,int b,int &x,int &y){
    if(!b){
        x=1; y=0;
        return;
    }
    int d=exgcd(b,a%b,y,x);
    y-=a/b*x;
    return d;
}

x=x/g*n,y=y/g*n;
```





### 5.2 线性同余方程

问题描述：
$$
求x使得ax\equiv b\;(mod\;m)，a,b,m都是整数
$$

$$
这就是一元线性同余方程的形式\\
它其实等价于ax+my=b\\
这就转化成了二元一次不定方程
$$

$$
方程有解的条件为：gcd(a,m)|b
$$



### 5.3 求逆元

逆元的定义：
$$
若整数 a，m 互质，并且对于任意的整数 b，如果满足 a|b，则存在一个整数 x，使得\\
b/a≡b×x\;\;(mod\;m)，则称 x 为 a 模m的乘法逆元，记作 a^{−1}\;\;(mod\;m)。
$$

$$
b/a*ax≡b×x\;\;(mod\;m)\\
由上式，x满足的方程为\;ax≡1\;\;(mod\;m)\\
这就要求gcd(a,m)|1\\
因此乘法逆元存在的充要条件是 a 与模数 m 互质
$$

$$
余数可加可减可乘\\
对于除法,要求b/a\;mod\;m\\
只要求出a模m的逆元x\\
就可以用乘法代替除法，有b/a≡b×x\;\;(mod\;m)
$$



问题描述：给定 n 组 **a~i~,p~i~**，求 **a~i~** 模 **p~i~** 的乘法逆元，若逆元不存在则输出 `impossible`。

**注意：**请返回在 **0∼p−1** 之间的逆元



#### 5.3.1 扩展欧几里得求逆元

用扩展欧几里得算法解线性同余方程 **ax ≡ 1 (mod p)**即可，要求**a,p互质**



***代码：***

```c++
#include<iostream>
using namespace std;

int exgcd(int a,int b,long long &x,long long &y){
    if(!b){
        x=1; y=0;
        return a;
    }
    int d=exgcd(b,a%b,y,x);
    y-=a/b*x;
    return d;
} 

int main(){
        
    int n;
    scanf("%d",&n);
    
    while(n--){
        int a,p;
        long long x,y;    //x，y用long long,可能会溢出
        scanf("%d%d",&a,&p);
        if(exgcd(a,p,x,y)!=1) printf("impossible\n");
        else printf("%d\n",(x%p+p)%p);
    }
    
    return 0;
}
```







#### 5.3.2 快速幂求逆元

$$
由费马小定理，当m是质数时,a^{m-1}≡1\;\;(mod\;m)\\
若a不是m的整数倍，只需 x=a^{m-2}
$$

**注意**：数据保证 **pi** 是质数。



只要 a 不是 p 的整数倍即存在逆元，快速幂求 **a^p-2^ mod p** 即可



***代码：***

```c++
#include<iostream>
using namespace std;

int qpow(int a,int k,int p){
    int ans=1;
    while(k){
        if(k&1) ans=1LL*ans*a%p;
        a=1LL*a*a%p; k>>=1;
    }
    return ans;
}

int main(){
        
    int n;
    scanf("%d",&n);
    
    while(n--){
        int a,p;
        scanf("%d%d",&a,&p);
        if(a%p==0) printf("impossible\n");
        else printf("%d\n",qpow(a,p-2,p));    //a^p-2%p
    }
    
    return 0;
}
```

______



  ## 6.中国剩余定理

问题描述：
$$
给定 2n 个整数 a_1,a_2,…,a_n 和 m_1,m_2,…,m_n，求一个最小的非负整数 x，满足 ∀i∈[1,n],x≡m_i\;(mod \;a_i)。\\
输出最小非负整数 x，如果 x 不存在，则输出 −1。
$$
问题的本质是求解**线性同余方程组**
$$
\begin{cases}
x≡m_1\;(mod \;a_1) \\
x≡m_2\;(mod \;a_2) \\
\;\;\;\;\;\;\;\;\;\;\;\;\vdots\\
x≡m_n\;(mod \;a_n)
\end{cases}
$$
转化成常规等式
$$
\begin{cases}
k_1*a_1+m_1=x \\
k_2*a_2+m_2=x \\
\;\;\;\;\;\;\;\;\;\;\;\;\vdots\\
k_n*a_n+m_n=x
\end{cases}
$$
联立先两个方程
$$
k_1*a_1-k_2*a_2=m_2-m_1\\
如果无解，说明不存在x\\
否则扩展欧几里得算法可以求出特解k_1^*,k_2^*以及a_1,a_2的最大公约数g\\
按照题目要求k_1^*或k_2^*应该是最小正整数解\\
进而可以写出k_1,k_2的通解\\
\begin{cases}
k_1=k_1^*+k*\frac{a_2}{g} \\
k_2=k_2^*+k*\frac{a_1}{g} \\
\end{cases}\\
进而x=k*\frac{a_1a_2}{g}+k_1^*a_1+m_1
$$
发现两个方程可以合并，新的 a 是 **原 a~1~，a~2~** **的最小公倍数**，新的 m 是**只考虑前两个方程 x 的特解**

循环处理这n个方程，两两合并，合并过程中出现矛盾即说明无解



***代码：***

```c++
#include<iostream>
using namespace std;

typedef long long LL;

LL exgcd(LL a,LL b,LL &x,LL &y){
    if(!b){
        x=1; y=0;
        return a;
    }
    LL d=exgcd(b,a%b,y,x);
    y-=a/b*x;
    return d;
}

int main()
{
    int n;
    scanf("%d",&n);
    
    LL a1,a2,m1,m2,k1,k2;
    bool flag=true;
    
    scanf("%lld%lld",&a1,&m1);
    for(int i=2;i<=n;i++){
        scanf("%lld%lld",&a2,&m2);
        
        LL g=exgcd(a1,a2,k1,k2);
        if((m2-m1)%g){
            flag=false;
            break;
        }
        
        k1*=(m2-m1)/g;
        LL t=a2/g;
        k1=(k1%t+t)%t;     
        
        m1=k1*a1+m1;
        a1=a1/g*a2;
    }
    
    if(flag) printf("%lld",m1);
    else printf("-1");
    
    return 0;
}
```

______



## 7.高斯消元

问题描述：**求解线性方程组**

输入一个包含 **n 个方程 n 个未知数**的线性方程组。

方程组中的系数为实数，求解这个方程组。以下是一个包含m个方程n个未知数的线性方程组示例
$$
\begin{cases}
a_{11}x_1+a_{12}x_2+\cdots+a_{1n}x_n=b_1  \\
a_{21}x_1+a_{22}x_2+\cdots+a_{2n}x_n=b_2  \\
\qquad\qquad\qquad\;\;\;\;\vdots\\
a_{m1}x_1+a_{m2}x_2+\cdots+a_{mn}x_n=b_m
\end{cases}
$$
 将系数矩阵化为行最简型，再化为等价标准型

步骤：

>1.枚举列，记录行，找到当前列的绝对值的最大值，最大值为0则跳过
>
>2.将这行与未被处理的最上面一行交换
>
>3.将交换后的这行同除以一个数，使首个非零数为1
>
>4.用这行倍加，将下面所有行的首位消成0
>
>5.这行记为已经被处理

最后判断系数矩阵的秩和增广矩阵的秩



>秩相等
>
>>小于n：无穷组解
>>
>>等于n：唯一解

>秩不等：无解



***代码：O(n^3^)***

```c++
#include<iostream>
#include<cmath>
using namespace std;

const int N=110;
const double eps=1e-8;

double a[N][N];

int gauss(int n){
    int c,r;
    for(c=r=1;c<=n;c++){
        int t=r; 
        for(int i=r;i<=n;i++)
            if(fabs(a[i][c])>fabs(a[t][c]))
                t=i; 
        if(fabs(a[t][c])<eps) continue; 
        for(int i=c;i<=n+1;i++) swap(a[r][i],a[t][i]);
        for(int i=n+1;i>=c;i--) a[r][i]/=a[r][c];
        for(int i=r+1;i<=n;i++)
            if(fabs(a[i][c])>eps)
                for(int j=n+1;j>=c;j--)
                    a[i][j]-=a[r][j]*a[i][c];
        r++; 
    }
    
    if(r<=n){
        for(int i=r;i<=n;i++)
            if(fabs(a[i][n+1])>eps)
                return 2;
        return 1;
    }
    
    for(int i=n;i>=1;i--){
        for(int j=i+1;j<=n;j++)
             a[i][n+1]-=a[i][j]*a[j][n+1],a[i][j]=0;
    }
    return 0;
}

int main(){
    
    int n;
    scanf("%d",&n);
    
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n+1;j++)
            scanf("%lf",&a[i][j]);
    
    int flag=gauss(n);
    
    if(flag==0){
        for(int i=1;i<=n;i++) 
            if(fabs(a[i][n+1])>eps) printf("%.2lf\n",a[i][n+1]);
            else printf("%.2lf\n",0.00);
    }
    else if(flag==1) printf("Infinite group solutions");
    else printf("No solution");
    
    return 0;
}
```





变形：**解线性异或方程组**

输入一个包含 n 个方程 n 个未知数的异或线性方程组。

方程组中的系数和常数为 0 或 1，每个未知数的取值也为 0 或 1。

求解这个方程组。

异或线性方程组示例如下：

```
M[1][1]x[1] ^ M[1][2]x[2] ^ … ^ M[1][n]x[n] = B[1]
M[2][1]x[1] ^ M[2][2]x[2] ^ … ^ M[2][n]x[n] = B[2]
…
M[n][1]x[1] ^ M[n][2]x[2] ^ … ^ M[n][n]x[n] = B[n]
```



步骤大体上相似，用位运算代替加减乘除即可，两个方程异或在一起即可消元



***代码：***

```c++
#include<iostream>
#include<cmath>
using namespace std;

const int N=110;

int a[N][N];

int gauss(int n){
    int c,r;
    for(c=r=1;c<=n;c++){
        int t=r; 
        for(int i=r;i<=n;i++)
            if(a[i][c]==1)
                t=i;
        if(!a[t][c]) continue; 
        for(int i=c;i<=n+1;i++) swap(a[r][i],a[t][i]);
        for(int i=r+1;i<=n;i++)
            if(a[i][c]==1)
                for(int j=c;j<=n+1;j++)
                    a[i][j]^=a[r][j];
        r++; 
    }
    
    if(r<=n){
        for(int i=r;i<=n;i++)
            if(a[i][n+1]==1)
                return 2;
        return 1;
    }
    
    for(int i=n;i>=1;i--){
        for(int j=i+1;j<=n;j++)
            a[i][n+1]^=a[i][j]*a[j][n+1];
    }
    return 0;
}

int main(){
    
    int n;
    scanf("%d",&n);
    
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n+1;j++)
            scanf("%d",&a[i][j]);
    
    int flag=gauss(n);
    
    if(flag==0){
        for(int i=1;i<=n;i++) 
            printf("%d\n",a[i][n+1]);
    }
    else if(flag==1) printf("Multiple sets of solutions");
    else printf("No solution");
    
    return 0;
}
```

______



## 8.组合数

$$
组合数C_a^b=\frac{a!}{b!(a-b)!}
$$

### 8.1 求组合数

根据数据规模，有以下四种类型

|     a,b的范围      | 询问次数 |        模数         |    方法    |
| :----------------: | :------: | :-----------------: | :--------: |
|  1 ≤ b ≤ a ≤ 2000  |  10^4^   |       10^9^+7       | 预处理答案 |
| 1 ≤ b ≤ a ≤ 10^5^  |  10^4^   |       10^9^+7       | 预处理阶乘 |
| 1 ≤ b ≤ a ≤ 10^18^ |    20    | 1 ≤ p ≤ 10^5^，质数 | 卢卡斯定理 |
|  1 ≤ b ≤ a ≤ 5000  |    1     |         无          |   高精度   |



#### 8.1.1 递推预处理答案

问题描述：
$$
给定 n 组询问，每组询问给定两个整数 a，b，请你输出 C_a^b\;mod\;(10^9+7) 的值。\\
1≤n≤10000,
1≤b≤a≤2000
$$

利用递推式预处理
$$
C_a^b=C_{a-1}^b+C_{a-1}^{b-1}
$$
***代码：***

```c++
#include<iostream>
using namespace std;

const int N=2010,mod=1e9+7;

int c[N][N];

void init(){
    for(int i=0;i<N;i++)
        for(int j=0;j<=i;j++)
            if(!j) c[i][j]=1;
            else c[i][j]=1LL*(c[i-1][j-1]+c[i-1][j])%mod;
}

int main(){
    
    int n;
    scanf("%d",&n);
    
    init();
    
    int a,b;
    while(n--){
        scanf("%d%d",&a,&b);
        printf("%d\n",c[a][b]);
    }
    
    return 0;
}
```



#### 8.1.2 预处理阶乘

问题描述：
$$
给定 n 组询问，每组询问给定两个整数 a，b，请你输出 C_a^b\;mod\;(10^9+7) 的值。\\
1≤n≤10000,
1≤b≤a≤10^5
$$



a，b的范围比较大，直接预处理出全部答案会超时

这就需要利用组合数计算的定义式，预处理阶乘的余数及其逆元



***代码：***

```c++
#include<iostream>
using namespace std;

const int N=100010,mod=1e9+7;

int fact[N]={1},infact[N]={1};

int qpow(int a,int k,int p){
    int ans=1;
    while(k){
        if(k&1) ans=1LL*ans*a%p;
        k>>=1; a=1LL*a*a%p;
    } 
    return ans; 
}

int main()
{
    for(int i=1;i<N;i++){
        fact[i]=1LL*fact[i-1]*i%mod;  
        infact[i]=1LL*infact[i-1]*qpow(i,mod-2,mod)%mod;
    }
 
    int n;
    scanf("%d",&n);
       
    while(n--){
        int a,b;
        scanf("%d%d",&a,&b);
        printf("%lld\n",1LL*fact[a]*infact[a-b]%mod*infact[b]%mod);    //防止爆long long
    }

    return 0;
}
```



#### 8.1.3 卢卡斯定理

$$
如果p是质数\\
则C_a^b\;\equiv\;C_{a/p}^{b/p}\;·\;C_{a\;mod\;p}^{b\;mod\;p}\;\;(mod\;p)
$$

证明：
$$
将a,b写成p进制数\\
a=a_0+a_1p+a_2p^2+……+a_kp^k\\
b=b_0+b_1p+b_2p^2+……+b_kp^k\\
定理：(1+x)^{p^k}\;\equiv\;1+x^{p^k}\;(mod\;p)\\

\begin{aligned}
考虑（1+x)^a&=(1+x)^{a_0+a_1p+a_2p^2+……+a_kp^k}\\
&=((1+x)^{p^k})^{a_k}·((1+x)^{p^{k-1}})^{a_{k-1}}·……·(1+x)^{a_0}\\
&\equiv(1+x^{p^k})^{a_k}·(1+x^{{p^{k-1}}})^{a_{k-1}}·……·(1+x)^{a_0}\;\;(mod\;p)\\
\end{aligned}\\
对任意x成立\\
对比左右两边x^b的系数，可得：\\
C_a^b\;\equiv\;C_{a^k}^{b^k}\;·\;C_{a^{k-1}}^{b^{k-1}}·……·\;C_{a^0}^{b^0}\;\;(mod\;p)\\
再由C_{a/p}^{b/p}\;\equiv\;C_{a^k}^{b^k}\;·\;C_{a^{k-1}}^{b^{k-1}}·……·\;C_{a^1}^{b^1}\;\;(mod\;p)\\
得C_a^b\;\equiv\;C_{a/p}^{b/p}\;·\;C_{a\;mod\;p}^{b\;mod\;p}\;\;(mod\;p)
$$


据此可以递归处理组合数，当a，b都小于p时直接循环用定义求组合数即可



问题描述：
$$
给定 n 组询问，每组询问给定三个整数 a,b,p，其中 p 是质数，请你输出 C_b^a\;mod\;p 的值。\\
1≤n≤20,
1≤b≤a≤10^{18},
1≤p≤10^5,
$$


***代码：***

```c++
#include<iostream>
using namespace std;

typedef long long LL;

int qpow(int a,int k,int p){
    int ans=1;
    while(k){
        if(k&1) ans=1LL*ans*a%p;
        k>>=1; a=1LL*a*a%p;
    }
    return ans;
}
int c(int a,int b,int p){
    int ans=1;
    for(int i=a,j=1;j<=b;i--,j++){
        ans=1LL*ans*i%p;
        ans=1LL*ans*qpow(j,p-2,p)%p;
    }
    return ans;
}
int lucas(LL a,LL b,int p){
    if(a<p&&b<p) return c(a,b,p);
    return 1LL*lucas(a/p,b/p,p)*lucas(a%p,b%p,p)%p; 
}

int main(){
    
    int n;
    scanf("%d",&n);
    
    LL a,b;  int p;
    while(n--){
        scanf("%lld%lld%d",&a,&b,&p);
        printf("%d\n",lucas(a,b,p));
    }
    
    return 0;
}
```



#### 8.1.4 高精度组合数

问题描述：
$$
输入 a,b，求 C_b^a 的值。\\
注意结果可能很大，需要使用高精度计算。\\
1≤b≤a≤5000
$$
以上三种求组合数类型都要求取余数，如果要求不取余，输出准确结果,就要用高精度

避免除法的方式是分解质因数，统计各个质因子出现的次数，出现在分子上加，出现在分母上减

最后只需高精乘以低精即可



***代码：***

```c++
#include<iostream>
#include<vector>
using namespace std;

const int N=5010;

int cnt,p[N],num[N];
bool vis[N];

void get_primes(int n){
    for(int i=2;i<=n;i++){
        if(!vis[i]) p[++cnt]=i;
        for(int j=1;p[j]<=n/i;j++){
            vis[i*p[j]]=true;
            if(i%p[j]==0) break;
        }
    } 
}
int get(int n,int p){    //1~n中共有几个因子p
    int res=0; 
    while(n){
        res+=n/p;
        n/=p;
    }
    return res;
}
vector<int> mul(vector<int> a,int b){
    vector<int> ans;
    int t=0;
    for(int i=0;i<a.size();i++){
        t+=a[i]*b;
        ans.push_back(t%10);
        t/=10;
    }
    while(t){
        ans.push_back(t%10);
        t/=10;
    }
    return ans;
}

int main(){
    
    int a,b;
    scanf("%d%d",&a,&b);
    
    get_primes(a);
     
    for(int i=1;i<=cnt;i++){
        num[i]+=get(a,p[i]);
        num[i]-=get(a-b,p[i]);
        num[i]-=get(b,p[i]);
    } 
    
    vector<int> ans; 
    ans.push_back(1);
    for(int i=1;i<=cnt;i++)
        for(int j=1;j<=num[i];j++)
            ans=mul(ans,p[i]);      
            
    for(int i=ans.size()-1;i>=0;i--) 
        printf("%d",ans[i]);
    
    return 0;
}
```





### 8.2 Catalan数

$$
Catalan数是一个数列\\
\begin{aligned}
C_n&=\frac{1}{n+1}C_{2n}^n\\
&=C_{2n}^n-C_{2n}^{n-1}\\
&=C_{2n}^n-C_{2n}^{n+1}\\\\
C_n&=C_0C_{n-1}+C_1C_{n-2}+……+C_{n-2}C_{1}+C_{n-1}C_{0}\\
&=\frac{4n-2}{n+1}C_{n-1}\;\;\;\;n=0,1,2…
\end{aligned}\\
$$



很多经典的问题都与卡特兰数有关：

>棋盘：给定n行n列的棋盘，从左下角走到右上角，**不经过对角线**，共有几种走法
>
>括号：n个左括号和n个右括号组成的字符串序列，共有多少种合法的组合
>
>出栈：给定一个字符串表示入栈序列，共有多少种可能的**出栈顺序**
>
>二叉树：n个节点构成的**二叉树**共有多少种情况



问题描述：

给定 n 个 0 和 n 个 1，它们将按照某种顺序排成长度为 2n 的序列，求它们能排列成的所有序列中，能够满足任意前缀序列中 0 的个数都不少于 1 的个数的序列有多少个。

输出的答案对 10^9^+7 取模。



这个问题等价于棋盘问题，用0，1分别表示棋盘上向右，向上走即可

将左下到右上的对角线向左上平移一格，以这条线为界，所有合法路线都在该线右下，绝对不接触该线

对于所有不合法的路线，其一定会经过该线，从第一次经过该线的位置开始，将后方的的路线沿该线对称

这些不合法的路线终点将变成（n-1，n+1），而合法的路线无法被映射

由此得到合法路线数量为**C(2n,n) - C(2n,n-1)**

任选卡特兰数计算公式即可



***代码：***

```c++
#include<iostream>
using namespace std;

const int mod=1e9+7;

int qpow(int a,int k,int p){
    int ans=1;
    while(k){
        if(k&1) ans=1LL*ans*a%p;
        a=1LL*a*a%p;
        k>>=1;
    }
    return ans;
}

int main(){
    
    int n;
    scanf("%d",&n);
    
    int ans=1;
    for(int i=n+1;i<=2*n;i++) ans=1LL*ans*i%mod;
    for(int i=1;i<=n;i++) ans=1LL*ans*qpow(i,mod-2,mod)%mod;
    ans=1LL*ans*qpow(n+1,mod-2,mod)%mod;
    
    printf("%d",ans);
    
    return 0;
}
```

****





## 9.容斥原理

问题描述：

给定一个整数 n 和 m 个不同的质数p1,p2,…,pm。

请你求出 1∼n 中能被p1,p2,…,pm 中的至少一个数整除的整数有多少个。



对于每个质数 **p~i~** ，**1~n** 的整数中能被它整除的数共有 **n/p~i~** 个

这样处理所有质数后，发现计数重复了，即是质数两两乘积的倍数的这些数被计算了两次，需要减掉 **n/p~i~p~j~**

如此，又发现即是三个质数乘积的倍数的这些数被多减了一次，又需要加回来 **n/p~i~p~j~p~k~** 

即是奇数个质数乘积的倍数的这些数，每次被计算到时，次数加1

是偶数个质数乘积的倍数的这些数，每次被计算到时，次数减1



用二进制数枚举所有质数乘积的组合情况



***代码：***

```c++
#include<iostream>
using namespace std;

const int N=20;

int p[N];


int main(){
    
    int n,m;
    scanf("%d%d",&n,&m);
    for(int i=0;i<m;i++) scanf("%d",&p[i]);
    
    int ans=0;
    for(int i=1;i<(1<<m);i++){
        long long t=1; int cnt=0;
        for(int j=0;j<m;j++){
            if(i>>j&1){
                t*=p[j];
                cnt++;
                if(t>n) break;
            }
        }
        if(t<=n){
            if(cnt&1) ans+=n/t;
            else ans-=n/t;
        }
    }
    
    printf("%d",ans);
    
    return 0;
}
```

______



## 10.博弈论

**公平组合游戏（ICG）：**

>1.有两个玩家，游戏规则对两个人是公平的
>
>2.游戏的状态有限，能走的步数也有限
>
>3.两人轮流走步，当一个玩家不能走步时游戏结束
>
>4.游戏的局势不能区分玩家身份，或者说可执行的合法行动与轮到哪位玩家行动无关（围棋不属于ICG）

这样的游戏有一个特征：***存在必胜策略***

即**给定初始局势**，**指定先手玩家**，如果**双方都采取最优策略**，那么**获胜者就已经确定了**



###10.1 尼姆游戏

问题描述：
$$
给定 n 堆石子，两位玩家轮流操作\\
每次操作可以从任意一堆石子中拿走任意数量的石子（可以拿完，但不能不拿），最后无法进行操作的人视为失败\\
问如果两人都采用最优策略，先手是否必胜\\
$$



 答案：
$$
假设有n堆石子，石子数目分别是a_1,a_2,…,a_n,\\
如果a_1⊕a_2⊕…⊕a_n≠0，先手必胜；否则，先手必败
$$



证明：在操作过程中，如果 **a~1~⊕a~2~⊕…⊕a~n~=x≠0**。那么玩家必然可以通过**拿走某一堆若干个石子将异或结果变为0**。
$$
不妨设x的二进制表示中最高一位1在第k位\\
那么在a_1,a_2,…,a_n中，必然有一个数a_i，它的第k位是1\\
由于a_i⊕x会将a_i第k位的1置成0，更高位不变，而更低位无论怎样变，都一定有a_i⊕x<a_i\\
那么从第i堆石子中拿走a_i−a_i⊕x（一定大于0）个石子，第i堆石子还剩a_i−(a_i−a_i⊕x)=a_i⊕x个\\
此时异或和为a_1⊕a_2⊕…⊕a_i⊕x⊕…⊕a_n=x⊕x=0
$$
证明：在操作过程中，如果 **a~1~⊕a~2~⊕…⊕a~n~=0**，那么**无论玩家怎么拿**，必然会**导致异或结果变为非0**。
$$
假设玩家从第i堆石子拿走若干个，异或和仍是0\\
不妨设还剩下a_i'个，因为不能不拿，所以0≤a_i'<a_i\\
根据假设，此时a_1⊕a_2⊕…⊕a_i'⊕…⊕a_n=0\\
而(a_1⊕a_2⊕…⊕a_i⊕…a_n)⊕(a_1⊕a_2⊕…⊕a_i'⊕…⊕a_n)=a_i⊕a_i'=0\\
则 a_i=a'，与假设0≤a_i'<a_i矛盾
$$
证明：操作到最后，即游戏的**终止局面**，**异或和也为0**，并且这个局面**一定会遇到**
$$
根据规则，每次操作一定会使石子总数减少，最终一定会把所有石子拿光\\
此时所有堆石子数都为0，异或和也为0
$$


根据这三个特性：

>1.如果初始局面各堆石子数量异或和**不为0**，**先手**一定可以从某堆石子**拿走若干个**，**使后手处在异或和为0的局面**
>
>此时后手**无论怎么拿**，都会使异或和**变为非0**，循环往复，后手一定会遇到终止局面，**先手必胜**

>2.如果初始局面各堆石子数量异或和**为0**，**先手无论怎么拿**，都会使异或和**变为非0**
>
>那么后手总可以通过拿走某堆若干个石子，**使先手处在异或和为0的局面**，循环往复，先手一定会遇到终止局面，**先手必败**



***代码：***

```c++
#include<iostream>
using namespace std;

int main(){
    
    int n;
    scanf("%d",&n);
    
    int x,res=0;
    for(int i=1;i<=n;i++){
        scanf("%d",&x);
        res^=x;
    }
    
    if(res) printf("Yes");
    else printf("No");
    
    return 0;
}
```

  



问题描述：**尼姆游戏的变形（台阶）**
$$
现在，有一个 n 级台阶的楼梯，每级台阶上都有若干个石子，其中第 i 级台阶上有 a_i 个石子(i≥1)。\\

两位玩家轮流操作，每次操作可以从任意一级台阶上拿若干个石子放到下一级台阶中（不能不拿）。\\

已经拿到地面上的石子不能再拿，最后无法进行操作的人视为失败。\\

问如果两人都采用最优策略，先手是否必胜。\\
$$
答案：
$$
对所有奇数级台阶上的石子数求异或和\\
如果a_1⊕a_3⊕…⊕a_{2k-1}≠0，先手必胜；否则，先手必败
$$


已知：
$$
如果奇数级台阶石子数异或和不为0，必然可以通过一步操作，使异或和变为0\\
如果异或和为0，无论怎样操作，必然使异或和变为非0\\
另外，经过有限步，所有石子都一定会被拿到地上（第0级台阶）\\
此时为终止局面，异或和为0
$$
必胜策略：

>1.初始状态下异或和不为0，先手第一步将异或和变为0
>
>> 若后手拿奇数级台阶上的石子，先手**也从奇数级上拿**，保持拿完后**异或和变回0**
>>
>> 若后手拿偶数级台阶上的石子，先手将**刚刚拿下**的石子**再往下拿一级**，**保持奇数级台阶上石子数不变**
>
>先手必胜



>2.初始状态下异或和为0
>
>> 若先手拿奇数级台阶上的石子，后手**也从奇数级上拿**，保持拿完后**异或和变回0**
>>
>> 若先手拿偶数级台阶上的石子，后手将**刚刚拿下**的石子**再往下拿一级**，**保持奇数级台阶上石子数不变**
>
>先手必败



***代码：***

```c++
#include<iostream>
using namespace std;

int main(){
    
    int n;
    scanf("%d",&n);
    
    int x,res=0;
    for(int i=1;i<=n;i++){
        scanf("%d",&x);
        if(i&1) res^=x;
    }
    
    if(res) printf("Yes");
    else printf("No");
    
    return 0;
}
```





###10.2 Sprague-Grundy函数

**有向图游戏：**
$$
一个有向无环图G(X,F)，X是点（局势）的非空集合，F是X上的函数\\
对于x∈X，有F(x)\subset X\\
对于给定的x∈X，F(x)表示玩家从x出发能够移动到的位置\\
如果F(x)为空，说明无法继续移动，称x是终点位置\\\\
一个玩家先走，起点是x_0,两人交替走步\\
在位置x，玩家可以选择移动到y点，y∈F(x)\\
位于终点位置的玩家，判负
$$

**SG函数**
$$
一个有向无环图G(X,F)中，把节点x的Sprague-Grundy函数值定义为sg(x)\\
sg(x)满足两个条件：\\
1.不能与它后继节点的sg值相同\\
2.sg(x)是满足条件1的最小非负整数\\
显然终点的sg值为0，据此可以求得所有节点的sg值
$$

由SG函数的定义可知：

>若**sg(x)≠0**，其**后继节点的sg值**必然包含所有从**0到sg(x)-1的自然数**，即只需一步就能走到sg(x)=0的节点上
>
> 若**sg(x)=0**，该节点**可能为终点**，如果不是终点，一定**存在sg值不为0的后继节点**，即只需要一步就能走到sg(x)≠0的节点上

_____



问题描述：**尼姆游戏的变形（集合）**
$$
给定 n 堆石子以及一个由 k 个不同正整数构成的数字集合 S。\\

现在有两位玩家轮流操作，每次操作可以从任意一堆石子中拿取石子，\\
每次拿取的石子数量必须包含于集合 S，最后无法进行操作的人视为失败。\\

问如果两人都采用最优策略，先手是否必胜。\\
$$
答案：
$$
设n堆石子的数量分别为，x_1,x_2,……,x_n\\
若sg(x_1)⊕sg(x_2)⊕……⊕sg(x_n)≠0，先手必胜；否则，先手必败
$$



**将n堆石子看作n个有向图**，记第k堆石子当前的数量为 x~k~ ，每个数量对应图中的一个节点

证明：在操作过程中，如果 **sg(x~1~)⊕sg(x~2~)⊕…⊕sg(x~n~)=x≠0**。那么玩家必然可以通过**在某个有向图中行走一步将异或结果变为0**。
$$
不妨设x的二进制表示中最高一位1在第k位\newline
那么在sg(x_1),sg(x_2),…,sg(x_n)中，必然有一个数sg(x_i)，它的第k位是1\newline
由于sg(x_i)⊕x会将sg(x_i)第k位的1置成0，更高位不变，而更低位无论怎样变，都一定有sg(x_i)⊕x<sg(x_i)\newline
而sg值为sg(x_i)⊕x的点一定是节点x_i的后继节点\newline
从点x_i走到该点\newline
此时异或和为sg(x_1)⊕sg(x_2)⊕…⊕sg(x_i)⊕x⊕…⊕sg(x_n)=x⊕x=0
$$
证明：在操作过程中，**如果  sg(x~1~)⊕sg(x~2~)⊕…⊕sg(x~n~)=x=0** ，那么**无论在图中怎么走 **，必然会**导致异或结果变为非0**。
$$
假设玩家在第i个有向图中走了一步，异或和仍是0\\
不妨设还新节点为x_i'，节点与其后继节点的sg值一定不同，必然有sg(x_i)≠sg(x_i')\\
根据假设，此时sg(x_1)⊕sg(x_2)⊕…⊕sg(x_i')⊕…⊕sg(x_n)=0\\
而(sg(x_1)⊕sg(x_2)⊕…⊕sg(x_i)⊕…⊕sg(x_n)⊕sg(x_1)⊕sg(x_2)⊕…⊕sg(x_i')⊕…⊕sg(x_n)=sg(x_i)⊕sg(x_i')=0\\
则 sg(x_i)=sg(x_i')，与假设sg(x_i)≠sg(x_i')矛盾
$$
证明：操作到最后，即游戏的**终止局面**，**异或和也为0**，并且这个局面**一定会遇到**
$$
最终在全部的有向图上都会走到终点，终点的sg值为0，因此异或和也为0
$$


根据这三个特性：

>1.如果初始局面sg值异或和**不为0**，**先手**一定可以**在某个有向图中走一步**，**使后手处在异或和为0的局面**
>
>此时后手**无论怎么拿，都会使异或和**变为非0**，循环往复，后手一定会遇到终止局面，**先手必胜

>2.如果初始局面sg值异或和**为0**，**先手无论怎么走**，都会使异或和**变为非0**
>
>那么后手总可以**在某个有向图中走一步**，**使先手处在异或和为0的局面**，循环往复，先手一定会遇到终止局面，**先手必败**



求sg函数值有两种方法，一种是递推，另一种是记忆化搜索

***代码（递推）：***

```c++
#include<iostream>
#include<cstring>
using namespace std;

const int N=10010;


int n,k,a[N];
int s[N],sg[N];

void get_sg(int n){
    for(int i=0;i<=n;i++){
        memset(s,0,sizeof(s));
        for(int j=1;j<=k;j++)
            if(i-a[j]>=0) s[sg[i-a[j]]]=1;
        for(int j=0;j<=i;j++)
            if(!s[j]) { sg[i]=j; break; }
    }
}

int main(){
    
   
    scanf("%d",&k);
    for(int i=1;i<=k;i++) scanf("%d",&a[i]);
    
    get_sg(10000);
    
    int res=0;
    scanf("%d",&n);
    for(int i=1,x;i<=n;i++){
        scanf("%d",&x); 
        res^=sg[x];
    }
    
    if(res) printf("Yes");
    else printf("No");
    
    return 0;
}
```



***代码（记忆化搜索）：***

```c++
#include<iostream>
#include<cstring>
#include<unordered_set>
using namespace std;

const int N=10010;

int n,k;
int a[N],f[N];

int sg(int x){
    if(f[x]!=-1) return f[x];
    
    unordered_set<int> s;
    for(int i=1;i<=k;i++)
        if(x-a[i]>=0) 
            s.insert(sg(x-a[i]));
        
    for(int i=0;;i++)
        if(!s.count(i))
            return f[x]=i;
}

int main()
{
    scanf("%d",&k);
    for(int i=1;i<=k;i++) scanf("%d",&a[i]);
    
    memset(f,-1,sizeof(f));
    int res=0;
    
    scanf("%d",&n);
    for(int i=1,x;i<=n;i++){
        scanf("%d",&x);
        res^=sg(x);
    }
    
    if(res) printf("Yes");
    else printf("No");
    
    return 0;
}
```



问题描述：**尼姆游戏的变形（拆分）**
$$
给定 n 堆石子，两位玩家轮流操作，每次操作可以取走其中的一堆石子，然后放入两堆规模更小的石子\\
（新堆规模可以为 0，且两个新堆的石子总数可以大于取走的那堆石子数），最后无法进行操作的人视为失败。\\
问如果两人都采用最优策略，先手是否必胜。
$$
**多个独立局面的SG值，等于这些局面SG值的异或和**

 

***代码：***

```c++
#include<iostream>
#include<cstring>
#include<unordered_set>
using namespace std;

const int N=110;

int n,f[N];

int sg(int x){
    if(f[x]!=-1) return f[x];
    
    unordered_set<int> s;
    for(int i=0;i<x;i++)
        for(int j=0;j<=i;j++)
            s.insert(sg(i)^sg(j));
            
    for(int i=0;;i++)
        if(!s.count(i)) 
            return f[x]=i;
}

int main(){
    
    scanf("%d",&n);
    
    memset(f,-1,sizeof(f));
    
    int res=0;
    for(int i=1,x;i<=n;i++){
        scanf("%d",&x);
        res^=sg(x);
    }
    
    if(res) printf("Yes");
    else printf("No");
    
    return 0;
}
```

______

