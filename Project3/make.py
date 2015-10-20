import os 

print "compiling...\n"
os.system("g++ -c -Wall -fPIC -I. -I/local/include/python2.4 -fopenmp -std=c++11 main.cpp")
print "linking...\n"
os.system("g++ -o prj3 -fopenmp  lib.o main.o  ")
print "running...\n"
os.system("./prj3")

print ("cleaning up...\n")
os.system("rm main.o")
#os.system("rm Data/*.txt")
