#ストリング法を実装する前に溶媒和自由エネルギーの地形を見てみるためのプログラム

import makeInputFile
import subprocess
import variables as var

numberOfGrid = var.numberOfGrit#numberOfGrid**2がメッシュの数になる
contour.makePoints(numberOfGrid)
print(contour.points)
fileContour = open("forContour.dat","w") #ファイルをまっさらにする
# fileContour.close()
# fileContour = open("forContour.dat","a")
numberOfPoints = (numberOfGrid + 1)*(numberOfGrid + 1) #gridの点の数
for i in range(numberOfPoints):
    makeInputFile.makeInputFile(contour.points[i][0],contour.points[i][1])
    subprocess.call(["zsh","start_program.sh"])
    fileFreeEnergy = open("freeEnergy.dat","r")
    freeEnergy = float(fileFreeEnergy.read())
    fileFreeEnergy.close()
    fileContour.write("{0:4.3f}".format(contour.points[i][0]))
    fileContour.write("   " + "{0:4.3f}".format(contour.points[i][1]))
    fileContour.write("   " + str(freeEnergy) + "\n")
fileContour.close()


