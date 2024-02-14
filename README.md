Thanks for your attention. This is the source code of our VLDB 2024 submission "Influence Maximization via Vertex Countering".

# Experimental Environment

The experiments are performed on a CentOS Linux serve (Release 7.5.1804) with Quad-Core Intel Xeon CPU (E5-2640 v4 @ 2.20GHz) and 128G memory. All the algorithms are implemented in C++17. The source code is compiled by GCC (10.2.1) under the -O3 optimization.

# Compile & Run

## Input

The input data should be stored in `./data`, we give a graph in the data for example.

The data for experiments can be found in [SNAP](http://snap.stanford.edu).

## Output

The experiment results will be saved in `./result`.

## MIC

### Compile


```shell
cd ./static/
g++ -o MIC MIC.cpp -std=c++17 -O3
```

### Run

```shell
./MIC dataset_file seed_file budget 𝜀 l
```

For example

```shell
./MIC example.txt ../data/seed.txt 1 0.3 1
```

## MIC+

### Compile

```shell
cd ./static/
g++ -o MIC+ MIC+.cpp -std=c++17 -O3
```

### Run

```shell
./MIC+ dataset_file seed_file budget 𝜀 l
```

For example

```shell
./MIC+ example.txt ../data/seed.txt 1 0.3 1
```


## MICP for Dynamic Graphs

### Compile


```shell
cd ./dynamic/
make
```

### Run

```shell
./main dataset_file seed_file budget 𝜀 l
```

For example

```shell
./main example.txt ../data/seed.txt 1 0.3 1
```

