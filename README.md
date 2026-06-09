# OTFS-ZCD MATLAB 模擬

本專案使用 MATLAB 實作：

```text
OTFS
+
ZCD
```

模擬以下：

```text
發射 → 傳輸 → 接收 → 還原
```

# 本專案基本流程如下：

```text

→ 隨機產生資料 
→ OTFS 調變 (ISFFT)
→ ZCD 波形 (時域波形生成)
→ 傳輸 (理想/無雜訊通道)
→ ZCD 接收 (載波移除 + 最小平方法估計)
→ OTFS 解調 (SFFT)
→ QAM 解調 
→ 還原原始資料
```

---

# 2. 模擬參數設定

| 參數符號 | 參數名稱 | 設定值 |
| :--- | :--- | :--- |
| $M$ | 子載波數量 (延遲網格數) | 16 |
| $N$ | 時間槽數量 (都卜勒網格數) | 8 |
| $M \times N$ | 每個 OTFS 訊框的總符號數 | 128 |
| $f_0$ | 基頻頻率 (Subcarrier Spacing) | 20 kHz |
| $T = 1/f_0$ | 符號週期時間 | 50 $\mu$ s |
| $f_s$ | 系統取樣頻率 | $1000 \times f_0 \times (M+1)$ |

---

# 3. OTFS 二維網格變換

### 3.1 延遲-都卜勒域資料生成
輸入資料在延遲-都卜勒網格 $D$ 上進行 QAM 映射，網格座標為 $k \in [0, M-1]$ (延遲軸) 與 $l \in [0, N-1]$ (都卜勒軸)，生成二維離散信號 $X_{DD}[k, l]$。

### 3.2 逆辛傅立葉轉換 (ISFFT)
ISFFT 將訊號從延遲-都卜勒域映射到時頻域網格 $X_{TF}[m, n]$，其中 $m$ 代表子載波索引， $n$ 代表時間槽索引。其二維離散數學公式定義為：

$$X_{TF}[m, n] = \frac{1}{\sqrt{MN}} \sum_{k=0}^{M-1} \sum_{l=0}^{N-1} X_{DD}[k, l] e^{j2\pi \left( \frac{nk}{M} - \frac{ml}{N} \right)}$$

在 MATLAB 中，此雙重求和轉換可透過對行做 IFFT、對列做 FFT 高效實現：
```matlab
X_TF = sqrt(M/N) * fft(ifft(X_DD, [], 1), [], 2);
```


---

# 4. ZCD 傳輸

在特定的時間槽 $n$ 中，將該時段的時頻符號向量 $\mathbf{a}_n = [a_1, a_2, \dots, a_M]^T$ 

 (即 $X_{TF}[ :, n]$ ) 映射為連續時間波形。傳輸訊號中引入了一個位於第 $M+1$ 個子載波位置的強正弦保護載波，振幅為 $A_c$。
 
### 4.1 發射訊號
合成的連續時域信號 $s_{tx}(t)$ 定義如下：

$$s_{tx}(t) = \sum_{k=1}^{M} 2\Re [ a_k e^{j2\pi k f_0 t} ] + 2\Re [ A_c e^{j2\pi(M+1)f_0 t} ]$$

展開實部後可表示為：
$$s_{tx}(t) = 2 \sum_{k=1}^{M} \left[ \Re\{a_k\}\cos(2\pi k f_0 t) - \Im\{a_k\}\sin(2\pi k f_0 t) \right] + 2A_c\cos(2\pi(M+1)f_0 t)$$

### 4.2 MATLAB 實作波形合成
```matlab
s_tx = zeros(size(t));
for k = 1:M
    s_tx = s_tx + 2 * real(ak(k) * exp(1j * 2 * pi * k * f0 * t));
end

% 疊加第 M+1 個子載波作為參考強載波
s_tx = s_tx + 2 * real(carrier_amp * exp(1j * 2 * pi * (M+1) * f0 * t));
```

---

# 5. ZCD 接收機

本模擬假設理想通道，接收訊號 $r(t) = s_{tx}(t)$。接收端利用子載波之間的**正交性 (Orthogonality)** 進行載波的分離與資料恢復。

### 5.1 強載波正交投影估計
由於子載波之間滿足正交條件：

$$
\frac{1}{T}\int_{0}^{T} e^{j2\pi k f_0 t} e^{-j2\pi m f_0 t} dt = \delta_{km} = \begin{cases}
1, & k=m \\
0, & k \neq m 
\end{cases}
$$


因此，將接收訊號 $r(t)$ 投影至第 $M+1$ 子載波空間，可消去前 $M$ 個資料項的干擾，提取出強載波振幅 $A_c$：

$$\hat{A}_c = \Re \left( \frac{1}{T} \int_{0}^{T} r(t) e^{-j2\pi(M+1)f_0 t} dt \right)$$


**MATLAB 離散時間近似（均值運算）：**
```matlab
carrier_est = real(mean(rx_signal(:) .* exp(-1j * 2 * pi * (M+1) * f0 * t(:))));
```

### 5.2 強載波消去 
自接收信號中扣除估計出的強載波成分，得到純資料波形 $y(t)$：
$$y(t) = r(t) - 2\hat{A}_c\cos(2\pi(M+1)f_0 t)$$

```matlab
y = rx_signal(:) - 2 * carrier_est * real(exp(1j * 2 * pi * (M+1) * f0 * t(:)));
```

### 5.3 傅立葉基底構建與最小平方法 (LS) 估計
將離散觀測時間點 $t_n$ 帶入，建立大小為 $K \times M$ 的傅立葉基底複數矩陣 $\mathbf{A}$（其中 $K$ 為取樣點數）：

$$\mathbf{A} = \begin{bmatrix} 
e^{j2\pi (1) f_0 t_1} & e^{j2\pi (2) f_0 t_1} & \dots & e^{j2\pi M f_0 t_1} \\
e^{j2\pi (1) f_0 t_2} & e^{j2\pi (2) f_0 t_2} & \dots & e^{j2\pi M f_0 t_2} \\
\vdots & \vdots & \ddots & \vdots \\
e^{j2\pi (1) f_0 t_K} & e^{j2\pi (2) f_0 t_K} & \dots & e^{j2\pi M f_0 t_K}
\end{bmatrix}$$

```matlab
A = zeros(length(t), M);
for k = 1:M
    A(:, k) = exp(1j * 2 * pi * k * f0 * t(:));
end
```

由於發射信號包含了正反頻率成分（即 $2\Re\{x\} = x + x^*$），接收端重建本質上是求解一個過定線性方程組（Overdetermined System）。我們利用廣義逆矩陣執行**最小平方法 (Least Squares)** 估計，還原出原始時頻符號
$\mathbf{\hat{a}}$：

$$\mathbf{\hat{a}} = \arg\min_{\mathbf{a}} \|\mathbf{y} - 2\Re\{\mathbf{A}\mathbf{a}\}\|^2 = (\mathbf{A}^H \mathbf{A})^{-1} \mathbf{A}^H \mathbf{y}$$

註：在無雜訊且取樣率足夠時， $\mathbf{A}^H \mathbf{A}$ 趨近於對角矩陣。*

**MATLAB 矩陣左除實作：**
```matlab
ak_rec = A \ y;
```

### 5.4 辛傅立葉轉換 (SFFT)
使用 SFFT 將還原的時頻符號 $X_{TF, rec}[m, n]$ 逆轉回延遲-都卜勒域：

$$X_{DD, rec}[k, l] = \frac{1}{\sqrt{MN}} \sum_{m=0}^{M-1} \sum_{n=0}^{N-1} X_{TF, rec}[m, n] e^{-j2\pi \left( \frac{nk}{M} - \frac{ml}{N} \right)}$$

**MATLAB 實作程式碼：**
```matlab
X_DD_rec = sqrt(N/M) * fft(ifft(X_TF_rec, [], 2), [], 1);
```

最後將 $$X_{DD, rec}[k, l]$$ 裡的每個資料經由QAM解調，還原原始資料 $X_int$ 

---

# 6. 實驗結果驗證

### 輸出：

```text
Trials         : 200
Mean BER       : 0.00000000
Median BER     : 0.00000000
Min BER        : 0.00000000
Max BER        : 0.00000000
Success Rate   : 100%
```

![image](output.png)


# 7. 改進方向

### 1. AWGN channel

加入雜訊：

```text
BER vs SNR
```

---

### 2. PAPR

比較：

```text
CP-OTFS
vs
ZCD-OTFS
```

---

### 3. 真實無線通道

加入：

- multipath
- delay spread
- Doppler spread

---


