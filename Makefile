
CXX := icpc
INCLUDE := -I $(NAG_ROOT)/include 
CFLAGS := -O0 --std=c++11 $(INCLUDE)
LDFLAGS := -L $(NAG_ROOT)/lib -lnagc_mkl 

main : main.cpp myrand.o
	$(CXX) $(CFLAGS) main.cpp -o main myrand.o $(LDFLAGS)

myrand.o : myrand.cpp
	$(CXX) $(CFLAGS) -c myrand.cpp -o myrand.o




.PHONY : clean



clean :
	rm -rf ./main
	rm -rf ./*.o
