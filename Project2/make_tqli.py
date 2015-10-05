import os 

print "compiling...\n"
os.system("g++ -c -Wall -fPIC -I. -I/local/include/python2.4 -c  main_tqli.cpp")
print "linking...\n"
os.system("g++ -o  prj2 lib.o main_tqli.o  ")
print "running...\n"
os.system("./prj2")

print "plotting functions...\n"
os.system("gnuplot plot.gp")

print ("cleaning up...\n")
os.system("rm main_tqli.o")
#os.system("rm Data/*.txt")
