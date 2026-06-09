## 數學原理與矩陣推導

為了模擬程式碼內部的每一道四則運算，我們設定以下參數進行完全展開：
* **子載波數量**: $M = 2$ （欲估計係數為 $a_1, a_2$）
* **基頻頻率**: $f_0 = 1 \text{ Hz}$
* **採樣點數**: $L = 3$ ，時間向量

$$
f_t =
\begin{bmatrix}
0 \\
0.25 \\
0.5
\end{bmatrix}
$$
* **載波振幅**: $A_c = 1$

---

### 步驟一：建立傅立葉基底矩陣 $A$ (每個元素逐一展開)
程式公式為 $A(n, k) = \exp(j \cdot 2\pi \cdot k \cdot f_0 \cdot t_n) = \cos(2\pi k f_0 t_n) + j\sin(2\pi k f_0 t_n)$。
矩陣大小為 $3 \times 2$：

* **第 1 欄 ($k=1$)，對應頻率 $1 \cdot f_0 = 1 \text{ Hz}$**:
  * $A(1,1) = \exp(j \cdot 2\pi \cdot 1 \cdot 1 \cdot 0) = \exp(0) = 1 + j0 = 1$
  * $A(2,1) = \exp(j \cdot 2\pi \cdot 1 \cdot 1 \cdot 0.25) = \exp(j\frac{\pi}{2}) = \cos(\frac{\pi}{2}) + j\sin(\frac{\pi}{2}) = 0 + j1 = j$
  * $A(3,1) = \exp(j \cdot 2\pi \cdot 1 \cdot 1 \cdot 0.5) = \exp(j\pi) = \cos(\pi) + j\sin(\pi) = -1 + j0 = -1$

* **第 2 欄 ($k=2$)，對應頻率 $2 \cdot f_0 = 2 \text{ Hz}$**:
  * $A(1,2) = \exp(j \cdot 2\pi \cdot 2 \cdot 1 \cdot 0) = \exp(0) = 1 + j0 = 1$
  * $A(2,2) = \exp(j \cdot 2\pi \cdot 2 \cdot 1 \cdot 0.25) = \exp(j\pi) = \cos(\pi) + j\sin(\pi) = -1 + j0 = -1$
  * $A(3,2) = \exp(j \cdot 2\pi \cdot 2 \cdot 1 \cdot 0.5) = \exp(j2\pi) = \cos(2\pi) + j\sin(2\pi) = 1 + j0 = 1$

**觀測矩陣 $A$ 的最終結果：**

$$
A=
\begin{bmatrix}
1 & 1 \\
j & -1 \\
-1 & 1
\end{bmatrix}
$$

---

### 步驟二：計算並移除載波

#### 為什麼要移除載波？

想像寄了一封信，這封信裡面裝著我們真正要估計的多載波資料（基頻訊號 $y_{\text{clean}}(t)$）：
$$y_{\text{clean}}(t) = \sum_{k=1}^{M} a_k e^{j 2\pi k f_0 t}$$

為了保護信件並送到遠方，郵局把它裝進了一個厚重的防護快遞箱裡，這個箱子就是位於第 $M+1$ 個頻率的高頻實數載波 $\text{carrier}(t)$：
$$\text{carrier}(t) = 2 \cdot \text{real}\left( A_c \cdot e^{j 2\pi (M+1) f_0 t} \right) = A_c e^{j 2 \pi (M+1) f_0 t} + A_c e^{-j 2 \pi (M+1) f_0 t} $$

當天線在接收端收到整個包裹（`rx_signal`）時，它同時包含了「信件」與「快遞箱」的疊加：
$$\text{rx}\_\text{signal}(t) = y_{\text{clean}}(t) + \text{carrier}(t)$$

如果我們不拆箱，直接用最小平方法 `A \ rx_signal` 去讀取，厚重的箱子（高頻載波）會造成嚴重的干擾，使能量產生**頻譜洩漏**，導致解出來的數據 $a_k$ 嚴重失真。

因此，**我們必須先用剪刀把快遞箱拆掉，也就是在接收端接收到信號後，先計算強載波，並將他從原信號中扣除，去除強載波的干擾：**
$$y(t) = \text{rx}\_\text{signal}(t) - \text{carrier}(t)$$

我們可以利用「正交性」來推回強載波：

兩個不同頻率的正弦波（例如 $m \cdot f_0$ 與 $n \cdot f_0$），只要在一個完整的週期 $T = 1/f_0$ 內相乘並積分（加總），只要 $m \neq n$，其結果就絕對等於 0。

只有當兩者頻率完全相同（ $m = n$ ）時，積分才不會是 0，而是會彈出一個正比於振幅的數值：

$$ 
\frac{1}{T} \int_0^T e^{j 2 \pi m f_0 t} \cdot e^{-j 2 \pi n f_0 t} dt = 
\begin{cases} 
1, & \text{if } m = n \\
0,& \text{if } m \neq n 
\end{cases} 
$$

並且我們知道強載波被刻意放在第 $(M+1) \cdot f_0$ 的頻率上。

接收端首先生成一個頻率同樣為 $(M+1)f_0$ 的反向旋轉向量： $e^{-j 2 \pi (M+1) f_0 t}$。

接收端再把收到的總訊號 $s_{\text{tx}}(t)$ 與這個內部向量相乘。這個動作在頻域上相當於「**把整個頻譜往左移動了 $(M+1)f_0$ 的距離**」。

這時，原本位於 $(M+1)f_0$ 的強載波，會被剛好移到 $0\text{ Hz}$ 的位置：

$$ A_c e^{j 2 \pi (M+1) f_0 t} \cdot e^{-j 2 \pi (M+1) f_0 t} = A_c e^0 = A_c $$

為了取得 $A_c$ ，接收再端對一個週期 $T$ 內的所有訊號進行求平均（即程式碼中的 `mean`）：
* **訊號去哪了？** 原本在 $1f_0 \sim Mf_0$ 的數據，被左移後變成了 $-Mf_0 \sim -1f_0$。因為它們依然是 $f_0$ 的整數倍頻率，在一個週期內求平均後，通通變成 0。
* **強載波的負頻率成分呢？** 原本在 $-(M+1)f_0$ 的成分，左移後變成 $-2(M+1)f_0$，求平均後也變成 0。
* **剩下的強載波：** 只有剛剛被移到 $0\text{ Hz}$ 的那個 $A_c$，因為它變成了穩定的直流訊號（不受時間 $t$ 影響），求平均後完好無損地保留了下來。

經過上述的乘法與平均，理論上我們就已經得到了強載波振幅 $A_c$ 了。

這樣就能求出加上強載波前的信號 $y(t) \approx y_{\text{clean}}(t)$，接下來交給最小平方法（Least-Squares）時，才能順利還原 $a_k$。

現在載波頻率定義在第 $M+1 = 3 \text{ Hz}$。
$$\text{carrier} = 2 \cdot \text{real}\left( 1 \cdot \exp(j \cdot 2\pi \cdot 3 \cdot 1 \cdot \mathbf{t}) \right) = 2 \cdot \cos(6\pi \mathbf{t})$$

1. $t_1 = 0 \implies 2 \cdot \cos(6\pi \cdot 0) = 2 \cdot \cos(0) = 2 \cdot 1 = 2$
2. $t_2 = 0.25 \implies 2 \cdot \cos(6\pi \cdot 0.25) = 2 \cdot \cos(1.5\pi) = 2 \cdot 0 = 0$
3. $t_3 = 0.5 \implies 2 \cdot \cos(6\pi \cdot 0.5) = 2 \cdot \cos(3\pi) = 2 \cdot (-1) = -2$

$$
\text{carrier}=
\begin{bmatrix}
2 \\
0 \\
-2
\end{bmatrix}
$$


#### 假設發射端原本的真實數據為 $a_1 = 3, a_2 = 5$，則天線收到的 `rx_signal` 應為：

$$
\text{rx}_{\text{signal}}=
A
\begin{bmatrix}
3 \\
5
\end{bmatrix}
+
\text{carrier}=
\begin{bmatrix}
1\cdot3+1\cdot5 \\
j\cdot3+(-1)\cdot5 \\
(-1)\cdot3+1\cdot5
\end{bmatrix}
+
\begin{bmatrix}
2 \\
0 \\
-2
\end{bmatrix}=
\begin{bmatrix}
8 \\
-5+3j \\
2
\end{bmatrix}
+
\begin{bmatrix}
2 \\
0 \\
-2
\end{bmatrix}=
\begin{bmatrix}
10 \\
-5+3j \\
0
\end{bmatrix}
$$



#### 程式執行減法，剝離載波，恢復基頻訊號 $y$：

$$
y=
\begin{bmatrix}
10 \\
-5+3j \\
0
\end{bmatrix}-
\begin{bmatrix}
2 \\
0 \\
-2
\end{bmatrix}=
\begin{bmatrix}
8 \\
-5+3j \\
2
\end{bmatrix}
$$

---

### 步驟三：最小平方法內部計算 (`A \ y` 的每一步加減乘除)
MATLAB 執行 `ak_rec = A \ y` 時，因為列數大於欄數（ $3 > 2$ ），內部會解標準方程（Normal Equations）： $\mathbf{a_k} = (A^H A)^{-1} A^H y$。

#### 1. 計算 $A^H$（共軛轉置矩陣：行列對調，虛部變號）

$$
A^H=
\begin{bmatrix}
1 & -j & -1 \\
1 & -1 & 1
\end{bmatrix}
$$

#### 2. 計算 $A^H A$
$$A^H A = \begin{bmatrix} 1 & -j & -1 \\ 1 & -1 & 1 \end{bmatrix} \begin{bmatrix} 1 & 1 \\ j & -1 \\ -1 & 1 \end{bmatrix}$$
* **左上元素**: $(1 \cdot 1) + (-j \cdot j) + (-1 \cdot -1) = 1 + (-j^2) + 1 = 1 + (-(-1)) + 1 = 1 + 1 + 1 = 3$
* **右上元素**: $(1 \cdot 1) + (-j \cdot -1) + (-1 \cdot 1) = 1 + j - 1 = j$
* **左下元素**: $(1 \cdot 1) + (-1 \cdot j) + (1 \cdot -1) = 1 - j - 1 = -j$
* **右下元素**: $(1 \cdot 1) + (-1 \cdot -1) + (1 \cdot 1) = 1 + 1 + 1 = 3$

$$
A^H A=
\begin{bmatrix}
3 & j \\
-j & 3
\end{bmatrix}
$$

#### 3. 計算 $A^H y$
$$A^H y = \begin{bmatrix} 1 & -j & -1 \\ 1 & -1 & 1 \end{bmatrix} \begin{bmatrix} 8 \\ -5+3j \\ 2 \end{bmatrix}$$
* **第 1 個元素**: 
  $$1 \cdot 8 + (-j) \cdot (-5+3j) + (-1) \cdot 2$$
  $$= 8 + (5j - 3j^2) - 2$$
  $$= 8 + 5j - 3(-1) - 2$$
  $$= 8 + 5j + 3 - 2 = 9 + 5j$$
* **第 2 個元素**: 
  $$1 \cdot 8 + (-1) \cdot (-5+3j) + 1 \cdot 2$$
  $$= 8 + 5 - 3j + 2 = 15 - 3j$$

$$
A^H y=
\begin{bmatrix}
9+5j \\
15-3j
\end{bmatrix}
$$

#### 4. 計算 $(A^H A)^{-1}$ 反矩陣
利用 $2 \times 2$ 矩陣求逆公式 

$$
\begin{bmatrix} 
a & b \\
c & d 
\end{bmatrix}^{-1} = 
\frac{1}{ad-bc}
\begin{bmatrix} 
d & -b \\
-c & a 
\end{bmatrix}
$$

* **行列式值 $\Delta$**: $3 \cdot 3 - (j \cdot -j) = 9 - (-j^2) = 9 - (-(-1)) = 9 - 1 = 8$
* **反矩陣結果**:

$$
(A^H A)^{-1}=
\frac{1}{8}
\begin{bmatrix}
3 & -j \\
j & 3
\end{bmatrix}
$$

#### 5. 最終解出 $\mathbf{a_k} = (A^H A)^{-1} A^H y$
$$
\mathbf{a_k} = 
\frac{1}{8} 
\begin{bmatrix} 
3 & -j \\ 
j & 3 
\end{bmatrix} 
\begin{bmatrix} 
9 + 5j \\
15 - 3j 
\end{bmatrix}
$$

* **計算 $a_1$**:
  $$a_1 = \frac{1}{8} \left[ 3 \cdot (9+5j) + (-j) \cdot (15-3j) \right]$$
  $$a_1 = \frac{1}{8} \left[ (27 + 15j) + (-15j + 3j^2) \right]$$
  $$a_1 = \frac{1}{8} \left[ 27 + 15j - 15j + 3(-1) \right]$$
  $$a_1 = \frac{1}{8} \left[ 27 - 3 \right] = \frac{1}{8} \cdot 24 = 3$$
* **計算 $a_2$**:
  $$a_2 = \frac{1}{8} \left[ j \cdot (9+5j) + 3 \cdot (15-3j) \right]$$
  $$a_2 = \frac{1}{8} \left[ (9j + 5j^2) + (45 - 9j) \right]$$
  $$a_2 = \frac{1}{8} \left[ 9j + 5(-1) + 45 - 9j \right]$$
  $$a_2 = \frac{1}{8} \left[ -5 + 45 \right] = \frac{1}{8} \cdot 40 = 5$$

$$
\mathbf{ak\_rec}=
\begin{bmatrix}
3 \\
5
\end{bmatrix}
$$

經由完二元一次正規方程組求解，演算法還原了初始發射係數 `[3, 5]`

---
