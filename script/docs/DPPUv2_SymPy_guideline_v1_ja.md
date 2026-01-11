# DPPUv2 実装ガイドライン (Engineering & Theoretical Guidelines)

- **Version:** 1.0 
- **Target:** AI Coding Assistants
- **Context:** Einstein-Cartan Theory on Curved Spacetime (S³×S¹, etc.)

本ドキュメントは、DPPUv2プロジェクトにおける数式処理（SymPy）の実装指針、および理論物理学的規約（Convention）を定めるものである。以下のルールを厳守すること。

⇒ [English version](DPPUv2_SymPy_guideline_v1.md)

-----

## 1\. 高速化エンジニアリング指針 (Optimization Rules)

SymPyを用いた積分計算における「処理時間爆発（2時間以上）」を防ぎ、数秒で完了させるための鉄則。

### Rule 1.1: `expand()` + `cancel()` 戦略の徹底

積分（`integrate`）を行う直前に、被積分関数に対して高コストな `simplify()` を使用してはならない。代わりに「展開と約分」を行うこと。

  * **Don't:**
    ```python
    density = simplify(density)  # NG 禁止：巨大な式の因数分解は計算コストが極大
    result = integrate(density, x)
    ```
  * **Do:**
    ```python
    density = cancel(expand(density))  # OK 推奨：多項式の和にすることで項別積分を誘発
    result = integrate(density, x)
    ```

### Rule 1.2: 積分後の中間簡約の抑制

多重積分（例：$\phi$ 積分 → $\theta$ 積分）の際、中間結果に対して過度な `simplify` を行わない。`cancel` 程度に留め、最終的な積分直前で再度 `expand` する方が高速である。

-----

## 2\. 理論実装指針 (Theoretical Implementation Rules)

曲がった時空（Curved Spacetime）において、テンソル演算の整合性を保つための鉄則。

### Rule 2.1: 配列インデックス操作の禁止 (Robust Method)

非対角計量や $g_{\mu\nu} \neq 1$ の環境下では、`T[mu, nu, lam]` のように配列のインデックスを入れ替える操作は、物理的なテンソルの添字操作（上げ下げ）と等価ではない。

  * **Don't:**
    ```python
    # NG 禁止：物理的に誤った成分を参照する恐れがある
    term = T_tensor[mu, nu, lam]
    ```
  * **Do:** 必ず計量 $g_{\mu\nu}$ を介して操作する。
    1.  全ての添字を下げて完全共変形 $T_{\lambda\mu\nu}$ を作る。
    2.  添字の入れ替え（Permutation）を行う。
    3.  必要に応じて計量で添字を上げる。

#### 正規直交フレーム基底での最適化

正規直交フレーム基底においては、計量が単位行列（$g_{ab} = \eta_{ab} = \text{diag}(1,1,1,1)$ または $\text{diag}(-1,1,1,1)$）であるため、添字の上げ下げ計算を省略し、直接成分計算を行うことで高速化する。

ただし、物理的定義（Hehl 1976）の符号パターン $(+1, +1, -1)$ は厳守すること。

```python
# ============================================================
# Golden Logic for Contortion (Frame Basis / v3 Standard)
# ============================================================
# Assumption: Metric is diagonal/identity (Orthonormal Frame)
# Therefore T^a_bc and T_abc behave identically in code logic.

K_tensor = MutableDenseNDimArray.zeros(dim, dim, dim)

for a in range(dim):
    for b in range(dim):
        for c in range(dim):
            # Formula: K_abc = (1/2)(T_abc + T_bca - T_cab)
            # Note: Using T[a,b,c] directly as T_abc
            
            term = (T_tensor[a, b, c] + T_tensor[b, c, a] - T_tensor[c, a, b])
            
            val = term * Rational(1, 2)
            
            if val != 0:
                K_tensor[a, b, c] = simplify(val)
```

### Rule 2.2: 自己無撞着性の検証 (Consistency Check)

EC接続 ($\Gamma_{\text{EC}}$) を構築した後、必ず以下の検証コードを実行し、ミスマッチが **0** であることを確認しなければならない。

```python
# Torsion Consistency Check
T_verify = Gamma_EC[lam, mu, nu] - Gamma_EC[lam, nu, mu] # Hehl定義
mismatch = count(simplify(T_verify - T_original) != 0)
assert mismatch == 0
```

-----

## 3\. 理論的規約 (Standard Conventions)

論文執筆時の混乱を防ぐため、**Hehl (1976) 標準** に準拠する。

### 3.1 Torsion Definition

捩れテンソル $T^\lambda_{\ \mu\nu}$ の定義：
$$T^\lambda_{\ \mu\nu} \equiv \Gamma^\lambda_{\ \mu\nu} - \Gamma^\lambda_{\ \nu\mu}$$
（注：テトラドからの計算時は $T^a_{\mu\nu} = \partial_\nu e^a_\mu - \partial_\mu e^a_\nu$ がこれに対応する）

### 3.2 Contortion Formula

上記のTorsion定義と整合するContortion $K^\lambda_{\ \mu\nu}$ の公式（Verified Formula）：

$$K_{\lambda\mu\nu} = \frac{1}{2} \left( T_{\lambda\mu\nu} + T_{\mu\nu\lambda} - T_{\nu\lambda\mu} \right)$$

  * 符号パターン: **$(+1, +1, -1)$**
  * 注意: この式は $T_{\lambda\mu\nu}$（全ての添字を下げたもの）に対して適用すること。

### 3.3 Einstein-Cartan Connection

$${\Gamma_{\text{EC}}}^\lambda_{\ \mu\nu} = {\Gamma_{\text{LC}}}^\lambda_{\ \mu\nu} + K^\lambda_{\ \mu\nu}$$

  * $\Gamma_{\text{LC}}$: Levi-Civita接続（Christoffel記号）
  * $K$: Contortion


