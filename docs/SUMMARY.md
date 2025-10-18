# プロジェクト概要

## 格子QCDハミルトニアン形式プログラム

### プロジェクト構成

```
lattice-qcd-hamiltonian/
├── src/                          # ソースコード
│   ├── modules/                  # Fortran90モジュール群
│   │   ├── mod_parameters.f90           # パラメータ管理
│   │   ├── mod_su_algebra.f90           # SU(Nc)群代数
│   │   ├── mod_lattice.f90              # 格子構造
│   │   ├── mod_gauge_field.f90          # ゲージ場
│   │   ├── mod_wilson_fermion.f90       # Wilsonフェルミオン
│   │   ├── mod_staggered_fermion.f90    # Staggeredフェルミオン
│   │   ├── mod_hamiltonian.f90          # ハミルトニアン
│   │   ├── mod_time_evolution.f90       # 時間発展
│   │   ├── mod_parallel.f90             # 並列化
│   │   ├── mod_tensor_network.f90       # TensorNetwork IF
│   │   └── mod_quantum_interface.f90    # 量子計算IF
│   └── main.f90                  # メインプログラム
│
├── input/                        # 入力ファイル
│   └── parameters.inp                   # デフォルトパラメータ
│
├── examples/                     # サンプル入力
│   ├── small_test.inp                   # 小規模テスト
│   ├── su3_test.inp                     # SU(3)テスト
│   └── hexagonal_test.inp               # 六角格子テスト
│
├── output/                       # 出力ディレクトリ
│
├── docs/                         # ドキュメント
│   ├── THEORY.md                        # 理論的背景
│   ├── PHYSICS.md                       # 物理的応用
│   └── SUMMARY.md                       # 本ファイル
│
├── Makefile                      # ビルドシステム
├── README.md                     # プロジェクトREADME
└── INSTALL.md                    # インストールガイド
```

### 実装された機能

#### 1. 格子構造（mod_lattice.f90）

✅ **実装完了**
- 正方格子（ハイパーキューブ）：任意次元対応
- 六角格子：2D/3D対応
- 周期境界条件
- サイト・リンク・プラケットの管理
- 近傍探索アルゴリズム

**主要関数**：
- `initialize_lattice()`: 格子構造の初期化
- `initialize_square_lattice()`: 正方格子生成
- `initialize_hexagonal_lattice()`: 六角格子生成
- `coords_to_index()`: 座標↔線形インデックス変換
- `get_plaquette_sites()`: プラケット構成サイト取得

#### 2. SU(Nc)群代数（mod_su_algebra.f90）

✅ **実装完了**
- 任意のカラー数Ncに対応
- Gell-Mann生成子の自動生成
- 構造定数の計算
- 群元素の演算（積、随伴、トレース）
- 行列指数関数

**特殊化**：
- SU(2): Pauli行列
- SU(3): Gell-Mann行列
- SU(N): 一般化生成子

**主要関数**：
- `initialize_su_algebra()`: 群構造初期化
- `generate_su2_matrices()`: SU(2)生成子
- `generate_su3_matrices()`: SU(3)生成子
- `calculate_structure_constants()`: f^abc計算
- `su_matrix_exp()`: exp(iαT^a)計算

#### 3. ゲージ場（mod_gauge_field.f90）

✅ **実装完了**
- リンク変数 U_μ(x) ∈ SU(Nc)
- 電場 E^a_μ(x)（正準運動量）
- Cold/Hot スタート
- プラケット計算
- 場の強さテンソル

**主要関数**：
- `initialize_gauge_field()`: 場の初期化
- `random_su_matrix()`: ランダムSU(Nc)生成
- `calculate_plaquette()`: プラケット計算
- `plaquette_average()`: 平均プラケット
- `field_strength_tensor()`: F_μν計算

#### 4. Wilsonフェルミオン（mod_wilson_fermion.f90）

✅ **実装完了**
- 4成分Diracスピノル
- Wilson-Diracオペレータ
- γ行列（Dirac表現）
- フェルミオン作用
- バイリニア演算子

**実装内容**：
- フェルミオン場：ψ(x, spin, color)
- Wilson parameter r = 1
- ホッピング項の計算
- Doublerの除去

**主要関数**：
- `initialize_wilson_fermion()`: 場の初期化
- `setup_gamma_matrices()`: γ行列定義
- `wilson_dirac_operator()`: D_W適用
- `wilson_fermion_action()`: ψ̄Dψ計算
- `fermion_bilinear()`: ψ̄Γψ計算

#### 5. Staggeredフェルミオン（mod_staggered_fermion.f90）

✅ **実装完了**
- 1サイト1成分
- Staggered位相因子
- 簡約されたDiracオペレータ
- カイラル凝縮
- 相関関数

**実装内容**：
- フェルミオン場：χ(x, color)
- 位相因子：η_μ(x) = (-1)^(x₁+...+x_{μ-1})
- 効率的な計算
- カイラル対称性

**主要関数**：
- `initialize_staggered_fermion()`: 場の初期化
- `setup_staggered_phases()`: 位相因子計算
- `staggered_dirac_operator()`: D_stag適用
- `staggered_condensate()`: ⟨χ̄χ⟩計算
- `staggered_pion_correlator()`: パイオン相関

#### 6. ハミルトニアン（mod_hamiltonian.f90）

✅ **実装完了**
- H = H_E + H_B + H_F
- 電場エネルギー：(g²/2)ΣE²
- 磁場エネルギー：プラケット項
- フェルミオンエネルギー
- ハミルトンの運動方程式

**実装内容**：
- エネルギー計算
- 力の計算（δH/δU）
- Staple計算
- エネルギー密度分布

**主要関数**：
- `calculate_hamiltonian()`: 全ハミルトニアン
- `electric_field_energy()`: H_E計算
- `magnetic_field_energy()`: H_B計算
- `calculate_equations_of_motion()`: dU/dt, dE/dt
- `calculate_gauge_force()`: ゲージ力

#### 7. 時間発展（mod_time_evolution.f90）

✅ **実装完了**
- Leapfrog積分法（シンプレクティック）
- Runge-Kutta 4次積分法
- エネルギー保存の追跡
- 観測量の記録
- 配置の保存/読込

**実装内容**：
- 時間ステップ管理
- 観測量履歴
- SU(Nc)への射影
- ファイル入出力

**主要関数**：
- `initialize_time_evolution()`: 初期化
- `evolve_leapfrog()`: Leapfrog発展
- `evolve_rk4()`: RK4発展
- `record_observables()`: 観測量記録
- `write_observables()`: データ出力

#### 8. 並列化（mod_parallel.f90）

✅ **実装完了**
- OpenMPスレッド並列化
- MPI分散メモリ並列化
- オプション選択可能
- 並列リダクション
- 同期バリア

**実装内容**：
- 条件付きコンパイル（#ifdef）
- ランク管理
- 通信プリミティブ

**主要関数**：
- `initialize_parallel()`: 並列環境初期化
- `finalize_parallel()`: 終了処理
- `parallel_sum_real()`: 並列和
- `parallel_barrier()`: 同期

#### 9. Tensor Networkインターフェース（mod_tensor_network.f90）

✅ **インターフェース実装**
- テンソル構造定義
- ボンド次元管理
- TEBD概念実装
- MPO構築準備

**注意**：完全な実装には外部ライブラリ（ITensor, TeNPyなど）が必要

**実装内容**：
- テンソル型定義
- サイトテンソル管理
- ゲージ場↔テンソル変換インターフェース

**主要関数**：
- `initialize_tensor_network()`: TN初期化
- `gauge_to_tensor()`: 変換（インターフェース）
- `tebd_evolution_step()`: TEBD概念
- `apply_two_site_gate()`: 2サイトゲート

#### 10. 量子計算インターフェース（mod_quantum_interface.f90）

✅ **インターフェース実装**
- 量子回路構造
- キュービットマッピング
- ゲート分解
- バックエンド選択（Qiskit/Cirq/Q#）

**注意**：実際の実行には量子計算ライブラリとの連携が必要

**実装内容**：
- 量子ゲート型定義
- 回路構築
- Trotterization
- SU(Nc)→量子ゲート変換

**主要関数**：
- `initialize_quantum_interface()`: 初期化
- `gauge_to_quantum()`: エンコード
- `encode_sun_matrix()`: SU(Nc)分解
- `quantum_time_evolution()`: 量子発展
- `execute_quantum_circuit()`: 実行（インターフェース）

#### 11. パラメータ管理（mod_parameters.f90）

✅ **実装完了**
- Namelist形式の入力ファイル
- 全パラメータの一元管理
- デフォルト値設定
- パラメータ表示

**管理項目**：
- 格子設定（種類、次元、サイズ）
- 物理パラメータ（g, m, Nc, β）
- フェルミオン設定
- 時間発展設定
- 並列化オプション
- 計算手法選択

#### 12. メインプログラム（main.f90）

✅ **実装完了**
- 全モジュールの統合
- 実行フロー制御
- 計算手法の分岐
- 結果の出力

**フロー**：
1. 初期化（パラメータ、格子、場）
2. 初期観測量の計算
3. 時間発展（Classical/TN/Quantum）
4. 最終観測量の計算
5. 結果の出力とクリーンアップ

### ビルドシステム

#### Makefile

✅ **実装完了**

**機能**：
- 依存関係の自動解決
- モジュールコンパイル順序の管理
- OpenMP/MPIオプション
- クリーンビルド
- ヘルプ表示

**使用例**：
```bash
make                              # 標準ビルド
make USE_OPENMP=yes USE_MPI=yes  # 完全並列版
make clean                        # クリーン
make help                         # ヘルプ
```

### ドキュメント

#### README.md
- プロジェクト概要
- 使用方法
- 出力説明

#### INSTALL.md
- 環境構築
- コンパイラ設定
- トラブルシューティング
- 最適化オプション

#### docs/THEORY.md
- 理論的背景
- 数式の詳細
- 格子QCDの定式化
- 群論の基礎

#### docs/PHYSICS.md
- 物理的応用
- 観測量の解釈
- 実験との比較
- 相図と相転移

#### docs/SUMMARY.md（本ファイル）
- プロジェクト全体の概要
- 実装状況
- ファイル構成

### テスト例

#### examples/small_test.inp
- 小規模4⁴格子
- SU(2)、短時間
- 動作確認用

#### examples/su3_test.inp
- 標準SU(3) QCD
- 8⁴格子
- Staggeredフェルミオン

#### examples/hexagonal_test.inp
- 六角格子
- 2次元
- 代替幾何のテスト

### 実行方法

```bash
# ビルド
make

# 小規模テスト実行
./lattice_qcd examples/small_test.inp

# 標準実行
./lattice_qcd input/parameters.inp

# 並列実行（OpenMP）
export OMP_NUM_THREADS=8
./lattice_qcd

# 並列実行（MPI）
mpirun -np 4 ./lattice_qcd
```

### 出力ファイル

#### output/observables.dat
時間発展中の物理量：
- ステップ番号
- 時刻
- 全エネルギー
- 平均プラケット

#### output/final_config.dat
最終配置（バイナリ）：
- リンク変数
- 電場
- メタデータ

### 拡張性

#### 新しい観測量の追加
1. `mod_hamiltonian.f90`に関数追加
2. `mod_time_evolution.f90`で記録
3. 出力ファイルに書き込み

#### 新しい格子構造
1. `mod_lattice.f90`に初期化ルーチン追加
2. `initialize_lattice()`で分岐

#### カスタムフェルミオン
1. 新しいモジュール作成
2. Diracオペレータ定義
3. `mod_hamiltonian.f90`に統合

#### 境界条件の変更
`mod_lattice.f90`の近傍計算を修正

### 性能

#### 計算量

**格子サイズ**: N⁴ (4次元の場合)
**リンク数**: 4N⁴
**計算量**: O(N⁴) per time step

#### メモリ使用量

**SU(3), 8⁴格子の場合**：
- リンク変数: 8⁴ × 4 × 9 × 16 bytes ≈ 1.9 MB
- 電場: 8⁴ × 4 × 8 × 16 bytes ≈ 1.7 MB
- フェルミオン（Wilson）: 8⁴ × 4 × 3 × 16 bytes ≈ 0.8 MB
- 合計: ~5-10 MB

**16⁴格子の場合**：
- 約16倍（~100 MB）

**32⁴格子の場合**：
- 約256倍（~1.5 GB）

#### 並列化効率

**OpenMP**：
- スレッド並列（共有メモリ）
- スケーラビリティ: ~8-16スレッド

**MPI**：
- プロセス並列（分散メモリ）
- 格子分割が必要
- より大規模な計算に適用

### 制限事項と今後の課題

#### 現在の制限

1. **Tensor Network**：完全実装にはITensorなど外部ライブラリが必要
2. **量子計算**：実際の実行にはQiskit等との連携が必要
3. **格子サイズ**：メモリとCPU時間の制約
4. **フェルミオン行列式**：動的フェルミオンは未実装

#### 今後の拡張

1. **動的フェルミオン**：
   - HMCアルゴリズム
   - Rational hybrid Monte Carlo
   
2. **改良作用**：
   - Symanzik改良
   - Stout smearing
   - Twisted mass fermion

3. **有限温度・密度**：
   - 時間方向の境界条件
   - 化学ポテンシャル
   
4. **トポロジー**：
   - インスタントン測定
   - トポロジカル感受率

5. **観測量の拡充**：
   - Wilson/Polyakovループ
   - グルーボールスペクトル
   - ハドロン質量

6. **最適化**：
   - GPU対応（CUDA/OpenCL）
   - より効率的なアルゴリズム
   - キャッシュ最適化

### まとめ

本プロジェクトは、格子QCDのハミルトニアン形式による実時間ダイナミクスシミュレーションの完全な実装を提供します。

**完成度**：
- コア機能：✅ 100%実装
- 並列化：✅ 100%実装
- 高度な手法：🔄 インターフェースのみ（外部ライブラリ連携が必要）

**用途**：
- 教育・研究
- グルーオンダイナミクスの研究
- 非平衡QCDの解析
- 新しいアルゴリズムのテストベッド

**コードの質**：
- モジュール化された設計
- 明確なインターフェース
- 拡張可能な構造
- 詳細なドキュメント

このプログラムは、格子QCD研究の出発点として、また教育ツールとして有用です。
