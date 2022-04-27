## strange_classic 

### 出题idea

> 最早的出题思路来源于很早之前复现图灵攻击Enigma密码机时的思路，故题目名称为strange_classic，而hint名则为记录了这一事件的图灵传记电影“模仿游戏”。

​		Enigma密码机的加密可以抽象为$P = F^{-1}RF$, 图灵攻击Enigma密码机时的思路是考虑在明密文对中寻找能够构成首位相同的链，从而得到部分$P^n(x) = x$ 的情况，此时即有$R^n(F(x)) = F(x)$，这样只需要爆破 $R$ 函数相关的设定，大大降低了爆破的复杂度。

​		参考这一思路，于是我想把Enigma密码机中的 $F, R$ 换一换。首先整个变换是在 $GF(N)$ 上的，这里的 $F$ 直接使用一个可逆变换即可，这里为了方便使用一个可逆矩阵，$R$ 则使用的是基于多项式的流密码，本身也可以表示为多项式对应友矩阵的幂。考虑如果找到足够多的向量 $x_i, i = 1, 2, \dots$， 使得 $R^n(F(x_i)) = k_i*F(x_i)$， 那么向量 $F(x_i)$ 则为矩阵 $R^n$ 的特征向量，$k_i$ 即为特征值，得到了n组线性无关且符合上述等式的 $x_i$ 后，即获得了 $R^n$ 的特征值。我们可以由此构造出与 $R^n$ 相似的对角矩阵 $T^n$ ，$R^n \sim T^n$， 有 $R \sim T$，故我们只需要把 $R^n$ 的特征值在模$N$上开$n$ 次根即可得到 $R$ 的特征值，从而得到 $R$ 的特征方程即得到strange_classic的flag。

​		但是有个问题，对于矩阵 $R$ 的特征向量 $r_i$，一定也为 $R^n$ 的特征向量，如果想要找到 $x_i$，则需要能区分其特征向量的特点。而想要区分这一点，需要将加密n次的向量和加密其余次数特征向量区分开来， 需要构造一个与 $R$ 矩阵可交换的矩阵 $G$，$G$ 会改变 $R$ 的特征向量，但 $G^n$ 不改变 $R$ 的特征向量，即有 $G \neq I, G^n = I$。

​		最简单的情况是直接给一个数量矩阵 $g * I$，只要 $g^n = 1$ 即可，但这样并不能起到区分的作用，所以不能为数量矩阵。经过推导发现，由于 $R$ 矩阵为友矩阵，故任意矩阵 $G$，设置好 $G$ 最后一列的数据，其余数据由 $R$ 中的 $c_i$ 约束直接得出。由此，可以构造出一个非数量矩阵的与 $R$ 矩阵可交换的矩阵 $G$，随机选取 $G$ 的最后一列数据，得到 $G$ 有 $\frac{N-1}{N}$ 的概率不为0。而由于 $G$ 与 $R$ 可交换，故 $G$ 的任意次幂与 $R$ 也可交换。接下来需要研究有限域上矩阵的阶。如果 $G_{i\times i}$ 在 $GF(N)$ 上有 $i$ 个不同的特征值，则存在 $m$ 使得 $G^m = I$。而由于 $G$ 只由最后一列，即 $i$ 个元素决定，可以藉由这 $i$ 个变量表示出特征方程。给定一个有 $i$ 个不同根的多项式，通过系数联立，则可以得到 $i$ 方程，由此解出矩阵 $G$。但由于是有限域上多元高次方程，本身就是困难问题，所以考虑将模数设的小一些，然后爆破部分变量寻找符合的矩阵，这样矩阵的阶为 $\phi(N)*N$ 的因子，最终选取 $G' = G^{\frac{N*\phi N}{n}}$，则有 $G'^n = I$。为了不要在生成密文上耗费太多时间，由于$G'$ 的设置，存在 $n\ | \ \phi(N)$ ，故在模 $N$ 上开 $n$ 次根时，会出现 $n-1$ 个假根，为了降低爆破复杂度，选取特征值时，需注意不同特征值的 $n$ 次根不相等，最终只需要进行一个 $GCD(m*n, \phi(N))^m$ 的爆破即可，最后选取 $N = 43, m = 9, n = 3$。

​		最终的加密为 $P(x) = F^{-1}((F\sdot x)\sdot R^m \sdot G')$，梳理一遍整体的破解思路：

1. 寻找明文密文对中能构成长度为n的链
2. 得到 $m$ 个 $x_j$ 满足 $((F\sdot x_j) \sdot R^{m*n} \sdot G'^n) = k_j(F\sdot x_j) \Rightarrow ((F \sdot x_j) \sdot R^{m*n}) = k_j(F\sdot x_j)$
3. 将 $m$ 个 $k_j$ 在模 $N$ 上开 $m*n$ 次根，爆破 $n^m$ 种 $R$ 特征值的可能，从而得到 flag。

​		之后开始思考怎么让选手会往预期的思路上去靠，因为不可能直接暴露 $G'$，所以只能提示 $G$ 与 $R$ 可交换，且 $G$ 的特征值均落在 $GF(N)$ 上，并且使用的 $G' = G^{\frac{N*(N-1)}{3}}$，也提示了 $n$ 为3。

​		由于给出了flag的hash值方便爆破，考虑选手直接爆破特征值的复杂度为 $43^9 > 2^{48}$，36h的比赛中大概是够了。同时密文给出了4倍的冗余数据，并且没出成交互题限制时间，希望能看到有意思的非预期。

### exp

代入数据即可

```python
from sage.all import *
from Crypto.Util.number import *
from hashlib import sha1


def is_k(v0, v1):
    judge = v0[0] * inverse(v1[0], N) % N
    for w in range(1, len(v1)):
        if v0[w] * inverse(v1[w], N) % N != judge:
            return False
    return judge


def find_k(plain_texts, cipher_texts):
    k = []
    for i in range(len(plain_texts)):
        if plain_texts[i] in cipher_texts and cipher_texts[i] in plain_texts:
            xi = plain_texts[cipher_texts.index(plain_texts[i])]
            xj = cipher_texts[plain_texts.index(cipher_texts[i])]
            tmp_k = is_k(xj, xi)
            if tmp_k:
                k.append(tmp_k)
    return k


def get_t(roots):
    PR = PolynomialRing(GF(N), 'x')
    x = PR.gens()[0]
    fx = 1
    for i in range(m):
        fx *= (x - roots[i])
    t = [int(-i % N) for i in fx.coefficients()[:-1]]
    return t


def get_possible_roots(y, e):
    roots = []
    for i in range(N):
        if i ** e % N == y:
            roots.append(i)
    return roots


N = 43
m = 9
n = 3
plain = 
cipher = 
hash = 'fd1f241a4d3ff9fc25d1e2480baa8b0c3b5a4559'
all_k = list(set(find_k(plain, cipher)))
assert len(all_k) == m
all_roots = [get_possible_roots(ki, n * m) for ki in all_k]

all_possible_t = []
for j in range(n**m):
    index = []
    tmp = j
    for _ in range(m):
        index.append(tmp % n)
        tmp //= n
    all_possible_t.append([all_roots[i][index[i]] for i in range(len(all_roots))])
for j in range(len(all_possible_t)):
    flag = 'MRCTF{%s}' % sha1(str(get_t(all_possible_t[j])).encode()).hexdigest()
    if sha1(flag.encode()).hexdigest() == hash:
        print('Got the flag:', flag)
        break
```

