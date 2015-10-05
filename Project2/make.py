import os 

print "compiling...\n"
os.system("g++ -c -Wall -fPIC -I. -I/local/include/python2.4 -c  main.cpp")
print "linking...\n"
os.system("g++ -o  prj2 lib.o main.o  ")
print "running...\n"
os.system("./prj2")

print "plotting functions...\n"
os.system("gnuplot plot.gp")

print ("cleaning up...\n")
os.system("rm main.o")
#os.system("rm Data/*.txt")
