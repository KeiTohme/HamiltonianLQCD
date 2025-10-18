# Lattice QCD - Hamiltonian Formalism

格子QCD（Quantum Chromodynamics）のハミルトニアン形式による実時間ダイナミクスシミュレーションプログラム

## 概要

このプログラムは、ハミルトニアン形式での格子QCDシミュレーションを実装しており、グルーオンの実時間ダイナミクスを研究するためのツールです。

### 主な特徴

1. **格子構造の選択**
   - 正方格子（ハイパーキューブ）
   - 六角格子（2D/3D）
   - 任意の時空次元に対応

2. **SU(Nc)ゲージ理論**
   - 任意の色数Ncに対応
   - Gell-Mann生成子と構造定数の自動生成
   - ゲージ不変性の保持

3. **フェルミオン定式化**
   - Wilson fermion（連続的スピン自由度）
   - Staggered fermion（格子スピン自由度）

4. **実時間発展**
   - Leapfrog（シンプレクティック）積分法
   - Runge-Kutta 4次積分法
   - ハミルトニアン（エネルギー）保存の追跡

5. **並列化**
   - OpenMPによるスレッド並列化
   - MPIによる分散メモリ並列化
   - オプションで選択可能

6. **高度な計算手法**
   - Tensor Networkインターフェース（TEBD, PEPSなど）
   - 量子コンピュータインターフェース（Qiskit, Cirq, Q#など）

## ディレクトリ構造

```
.
├── src/
│   ├── modules/              # Fortran90モジュール
│   │   ├── mod_parameters.f90
│   │   ├── mod_su_algebra.f90
│   │   ├── mod_lattice.f90
│   │   ├── mod_gauge_field.f90
│   │   ├── mod_wilson_fermion.f90
│   │   ├── mod_staggered_fermion.f90
│   │   ├── mod_hamiltonian.f90
│   │   ├── mod_time_evolution.f90
│   │   ├── mod_parallel.f90
│   │   ├── mod_tensor_network.f90
│   │   └── mod_quantum_interface.f90
│   └── main.f90              # メインプログラム
├── input/
│   └── parameters.inp        # 入力パラメータファイル
├── output/                   # 出力ディレクトリ
├── Makefile                  # ビルドシステム
└── README.md                 # このファイル
```

## コンパイル

### 必要な環境

- Fortran90/95対応コンパイラ（gfortran, ifortなど）
- (オプション) OpenMP対応コンパイラ
- (オプション) MPI実装（OpenMPI, MPICHなど）

### ビルド方法

基本的なビルド（OpenMP有効）：
```bash
make
```

オプション指定：
```bash
# OpenMPなし
make USE_OPENMP=no

# MPI有効
make USE_MPI=yes

# 両方有効
make USE_OPENMP=yes USE_MPI=yes
```

クリーンビルド：
```bash
make clean
make
```

ヘルプ表示：
```bash
make help
```

## 使用方法

### 1. パラメータ設定

`input/parameters.inp`ファイルを編集して、シミュレーションパラメータを設定します。

```fortran
&lattice_config
  lattice_type = 'square'      ! 'square' or 'hexagonal'
  d_euclid = 4                 ! 時空次元
  nsize = 8, 8, 8, 8          ! 各次元の格子サイズ
/

&physical_params
  gauge_coupling = 1.0         ! ゲージ結合定数 g
  fermion_mass = 0.1          ! フェルミオン質量 m
  n_colors = 3                ! カラー数 Nc
  beta = 6.0                  ! 逆結合定数 β = 2Nc/g²
/

&fermion_config
  fermion_type = 'wilson'      ! 'wilson' or 'staggered'
  n_flavors = 2               ! フレーバー数
/

&time_evolution
  time_step = 0.01            ! 時間刻み dt
  n_steps = 1000              ! ステップ数
  output_interval = 10        ! 出力間隔
/

&parallel_config
  use_openmp = .true.         ! OpenMP使用
  use_mpi = .false.           ! MPI使用
  n_threads = 4               ! OpenMPスレッド数
/

&computation_method
  method = 'classical'         ! 'classical', 'tensor_network', 'quantum'
  tensor_chi = 50             ! Tensor Networkボンド次元
  quantum_backend = 'qiskit'  ! 量子計算バックエンド
/
```

### 2. 実行

```bash
./lattice_qcd
```

または、カスタム入力ファイルを指定：
```bash
./lattice_qcd my_input.inp
```

OpenMP並列実行（スレッド数指定）：
```bash
export OMP_NUM_THREADS=8
./lattice_qcd
```

MPI並列実行：
```bash
mpirun -np 4 ./lattice_qcd
```

### 3. 出力

プログラムは以下のファイルを`output/`ディレクトリに生成します：

- `observables.dat`: 時間発展中の物理量（エネルギー、プラケットなど）
- `final_config.dat`: 最終的なゲージ場配置（バイナリ形式）

## 物理的背景

### ハミルトニアン形式

格子QCDのハミルトニアンは以下のように表されます：

```
H = H_E + H_B + H_F

H_E = (g²/2) Σ_{x,μ,a} [E^a_μ(x)]²       (電場エネルギー)
H_B = (1/g²) Σ_{x,μ<ν} [1 - Re Tr U_P]   (磁場エネルギー)
H_F = Σ_x ψ̄(x) D_W ψ(x)                  (フェルミオン項)
```

ここで：
- `E^a_μ(x)`: 電場（正準運動量）
- `U_P`: プラケット（格子上の最小ループ）
- `D_W`: Wilson-Diracオペレータ
- `g`: ゲージ結合定数

### 時間発展

ハミルトンの運動方程式：
```
dU/dt = δH/δE
dE/dt = -δH/δU
```

これらをLeapfrog法やRunge-Kutta法で数値的に解きます。

### SU(Nc)群

SU(Nc)群の生成子T^aは以下を満たします：
```
[T^a, T^b] = i f^abc T^c
Tr(T^a T^b) = (1/2) δ^ab
```

リンク変数は群元素：
```
U_μ(x) ∈ SU(Nc)
U_μ(x) = exp(i Σ_a α^a_μ(x) T^a)
```

## 技術的詳細

### Wilsonフェルミオン

Wilson-Diracオペレータ：
```
D_W = m + (4 + m) - (1/2) Σ_μ [(1-γ_μ)U_μ(x)δ_{x+μ,y} + (1+γ_μ)U†_μ(x-μ)δ_{x-μ,y}]
```

γ行列はDirac表現を使用。

### Staggeredフェルミオン

位相因子η_μ(x) = (-1)^(x₁+...+x_{μ-1})を用いた簡約表現：
```
D_stag = m χ(x) + (1/2) Σ_μ η_μ(x) [U_μ(x)χ(x+μ) - U†_μ(x-μ)χ(x-μ)]
```

スピン自由度が格子に埋め込まれ、1サイトあたり1成分。

### Tensor Network法

- ボンド次元χで精度を制御
- TEBD（Time Evolution Block Decimation）によるシンプレクティック時間発展
- PEPSやMPSによる効率的な状態表現

### 量子計算インターフェース

- SU(Nc)群元素を量子ゲートに変換
- Trotterization: `exp(-iHt) ≈ exp(-iH_E t/2) exp(-iH_B t) exp(-iH_E t/2)`
- バックエンド: Qiskit (IBM), Cirq (Google), Q# (Microsoft)

## 拡張可能性

### カスタム観測量の追加

`mod_hamiltonian.f90`にカスタム観測量関数を追加：

```fortran
function my_observable() result(obs)
  real(dp) :: obs
  ! 計算コード
end function my_observable
```

### 新しい格子構造

`mod_lattice.f90`に新しい格子初期化ルーチンを追加：

```fortran
subroutine initialize_custom_lattice()
  ! カスタム格子の定義
end subroutine initialize_custom_lattice
```

### 異なる境界条件

周期的、開放、反周期的境界条件などを`mod_lattice.f90`で実装可能。

## 参考文献

1. Kogut, J., & Susskind, L. (1975). "Hamiltonian formulation of Wilson's lattice gauge theories." Physical Review D, 11(2), 395.

2. Wilson, K. G. (1974). "Confinement of quarks." Physical review D, 10(8), 2445.

3. Stacey, R. (1982). "Eliminating lattice fermion doubling." Physical Review D, 26(2), 468.

4. Troyer, M., & Wiese, U. J. (2005). "Computational complexity and fundamental limitations to fermionic quantum Monte Carlo simulations." Physical review letters, 94(17), 170201.

5. Martinez, E. A., et al. (2016). "Real-time dynamics of lattice gauge theories with a few-qubit quantum computer." Nature, 534(7608), 516-519.

## ライセンス

このプログラムは教育・研究目的で自由に使用できます。

## 作者

Lattice QCD Hamiltonian Formalism Project

## 問い合わせ

バグ報告や機能要望は、GitHubのIssueトラッカーまでお願いします。

---

**注意**: Tensor NetworkおよびQuantum Computingの完全な実装には外部ライブラリ（ITensor, TeNPy, Qiskit等）との連携が必要です。現在の実装はインターフェース層のみを提供しています。
