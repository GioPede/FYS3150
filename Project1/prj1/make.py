import os 

print "compiling...\n"
os.system("g++ -c -Wall -fPIC -I. -O3 -I/local/include/python2.4 -c  main.cpp")
print "linking...\n"
os.system("g++ -o  prj1 lib.o main.o  ")
print "running...\n"
os.system("./prj1")

print "plotting functions...\n"
os.system("gnuplot plot.gp")
print "plotting errors...\n"
os.system("gnuplot plot_err.gp")
os.system("gnuplot error.gp")
