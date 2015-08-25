import os 

os.system("g++ -c -Wall derive.cpp")
os.system("g++ -o derive derive.o -larmadillo")
os.system("./derive 10")
os.system("gnuplot plot.gp")
