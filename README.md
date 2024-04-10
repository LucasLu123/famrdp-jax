## 系统要求
- Ubuntu 22.04

## 安装工具链
```
sudo apt install cmake gfortran mpich
```

## 克隆仓库
```
git clone http://10.136.126.255/xxx/famrdp.git
```

## 编译
```
cd famrdp
cmake . -B build
cmake --build build -j 16
```

## 运行
```
cp param.inp test3.grd build/src/
./build/src/HOSTA.mpi
```