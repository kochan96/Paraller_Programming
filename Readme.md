# Build
1. go to directory with Makefile
2. run `make`

# Usage
```bash
Usage ./main -n size [-m max_number] [-d] [-o file]
Options are:
    -h: display what you are reading now
    -n size: size of matrix
    -m max_number: maximum number of cell in generated matrix (default 100)
    -v: display A,Q,R, A=Q*R matrices (default false)
    -o file: append size and execution time to file (in csv format)
```

#### Usage examples
```
Calculate QR decomposition of 4x4 matrix and display size and time
./main -n 4
```
```
Calculate QR decomposition of 5x5 matrix and append size and time to test.csv file
./main -n 5 -o test.csv
```
```
Calculate QR decomposition of 10x10 matrix and display size, time, Q matrix and R matrix
./main -n 10 -v
```