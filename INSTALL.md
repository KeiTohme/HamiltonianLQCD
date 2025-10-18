# インストールガイド

## 必要な環境

### Fortranコンパイラ

以下のいずれかのFortran90/95対応コンパイラが必要です：

#### GFortran（推奨）
```bash
# Ubuntu/Debian
sudo apt-get install gfortran

# CentOS/RHEL
sudo yum install gcc-gfortran

# macOS (Homebrew)
brew install gcc

# macOS (MacPorts)
sudo port install gcc12
```

#### Intel Fortran (ifort)
Intel oneAPI HPC Toolkitをインストール：
https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html

### オプション：並列化ライブラリ

#### OpenMP
- GFortran 4.2以降で自動的にサポート
- コンパイラフラグ: `-fopenmp`

#### MPI（オプション）
```bash
# Ubuntu/Debian
sudo apt-get install libopenmpi-dev openmpi-bin

# CentOS/RHEL
sudo yum install openmpi openmpi-devel

# macOS (Homebrew)
brew install open-mpi
```

## ビルド手順

### 1. ソースコードの取得

```bash
cd /path/to/lattice_qcd
```

### 2. コンパイル

基本的なビルド（OpenMP有効）：
```bash
make
```

並列化オプション付きビルド：
```bash
# OpenMPのみ
make USE_OPENMP=yes

# MPIのみ
make USE_MPI=yes

# 両方
make USE_OPENMP=yes USE_MPI=yes
```

### 3. コンパイラの変更

Makefileの先頭でコンパイラを変更できます：

```makefile
# GFortranの場合
FC = gfortran

# Intel Fortranの場合
FC = ifort

# NAG Fortranの場合
FC = nagfor
```

または、コマンドラインから：
```bash
make FC=ifort
```

### 4. 最適化オプション

デバッグビルド：
```bash
make FFLAGS="-g -O0 -fcheck=all -fbacktrace"
```

最適化ビルド：
```bash
make FFLAGS="-O3 -march=native -funroll-loops"
```

Intel Fortranでの高度な最適化：
```bash
make FC=ifort FFLAGS="-O3 -xHost -ipo -parallel"
```

## ビルドの確認

コンパイルが成功すると、実行ファイル`lattice_qcd`が作成されます：

```bash
./lattice_qcd --help  # ヘルプ表示（実装されている場合）
```

## トラブルシューティング

### エラー: "gfortran: command not found"

**原因**: Fortranコンパイラがインストールされていない

**解決策**:
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install gfortran

# または他のパッケージマネージャーを使用
```

### エラー: "mpif90: command not found"

**原因**: MPIがインストールされていないか、パスが通っていない

**解決策**:
```bash
# MPIのインストール
sudo apt-get install libopenmpi-dev

# または、MPI無効でビルド
make USE_MPI=no
```

### エラー: モジュール依存関係のエラー

**原因**: 不完全なビルド

**解決策**:
```bash
make clean
make
```

### 警告: "warning: unused variable"

これらは通常、動作には影響しません。警告を無視する場合：
```bash
make FFLAGS="-O3 -w"
```

### メモリ不足エラー

大規模な格子での計算時にメモリ不足になる場合：

1. 格子サイズを小さくする
2. より多くのRAMを持つマシンを使用
3. MPIで分散メモリを使用

```bash
# MPIで実行
mpirun -np 4 ./lattice_qcd
```

## 環境変数

### OpenMPスレッド数の設定

```bash
export OMP_NUM_THREADS=8
./lattice_qcd
```

### MPIプロセス数の設定

```bash
mpirun -np 16 ./lattice_qcd
```

### スタックサイズの増加（大規模計算時）

```bash
# bashの場合
ulimit -s unlimited

# または環境変数で
export OMP_STACKSIZE=2G
```

## 性能チューニング

### コンパイラ最適化フラグ

#### GFortran
```bash
make FFLAGS="-O3 -march=native -funroll-loops -ftree-vectorize"
```

#### Intel Fortran
```bash
make FC=ifort FFLAGS="-O3 -xHost -ipo -parallel -qopt-report"
```

### プロファイリング

gprof用：
```bash
make FFLAGS="-O3 -pg" LDFLAGS="-pg"
./lattice_qcd
gprof ./lattice_qcd gmon.out > profile.txt
```

## Docker環境（オプション）

Dockerfileの例：

```dockerfile
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    gfortran \
    make \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /lattice_qcd
COPY . .

RUN make USE_OPENMP=yes

CMD ["./lattice_qcd"]
```

ビルドと実行：
```bash
docker build -t lattice-qcd .
docker run -v $(pwd)/output:/lattice_qcd/output lattice-qcd
```

## クラスタ環境での実行

### SLURM例

```bash
#!/bin/bash
#SBATCH --job-name=lattice_qcd
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00

module load gcc/11.2.0
module load openmpi/4.1.1

mpirun ./lattice_qcd input/parameters.inp
```

### PBS例

```bash
#!/bin/bash
#PBS -N lattice_qcd
#PBS -l nodes=4:ppn=16
#PBS -l walltime=24:00:00

cd $PBS_O_WORKDIR
module load gcc openmpi

mpirun ./lattice_qcd
```

## 動作確認

小規模テスト実行：

```bash
# input/parameters.inpを編集して小さい格子に設定
# nsize = 4, 4, 4, 4
# n_steps = 10

./lattice_qcd

# 正常に動作すれば output/ ディレクトリに結果が出力される
ls -lh output/
```

## サポート

問題が解決しない場合は、以下の情報と共に報告してください：

1. OSとバージョン
2. コンパイラとバージョン (`gfortran --version`)
3. エラーメッセージ全文
4. 使用したコンパイルコマンド
