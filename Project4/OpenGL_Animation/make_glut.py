import os 

print "compiling...\n"
os.system("g++ -c -Wall -fPIC -I. -I/local/include/python2.4 -fopenmp -std=c++11 glut.cpp -lGL -lGLU -lglut")
print "linking...\n"
os.system("g++ -o prj4 -fopenmp   glut.o -lGL -lGLU -lglut ")
print "running...\n"
os.system("./prj4")

print ("cleaning up...\n")
os.system("rm glut.o")
#os.system("rm Data/*.txt")
