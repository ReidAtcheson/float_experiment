#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=272
numactl -m1 ./main 1
numactl -m1 ./main 10
numactl -m1 ./main 100
numactl -m1 ./main 1000
numactl -m1 ./main 10000
numactl -m1 ./main 100000
numactl -m1 ./main 1000000
