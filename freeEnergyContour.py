zeropoint = 5.582 + (-75.98)

f = open("freeEnergyContour.dat","w")

fInput = open("forContour.dat","r")

while True:
    readLine = fInput.readline()
    if readLine == "":
        break
    else:
        x, y, z = map(float, readLine.split())
        f.write("{0:4.3f}".format(x))
        f.write("   " + "{0:4.3f}".format(y))
        f.write("   " + str(z - zeropoint) + "\n")
#x,y,z = map(float, fInput.readline().split())
