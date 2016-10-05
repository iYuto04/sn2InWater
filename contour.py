def makePoints(numberOfGrid):
#x軸もy軸も同じスケールのものをつくる
    startPosition = 1.5
    endPosition = 4.0
    gridOfRow = numberOfGrid
    gridOfLine = numberOfGrid
    dtOfLine = (endPosition - startPosition) / gridOfRow
    dtOfRow = dtOfLine

    global points
    points = []

    for i in range(gridOfLine + 1): #pythonのrangeの仕様上+1している
        for j in range(gridOfRow + 1):
            x = startPosition + dtOfLine*i
            y = startPosition + dtOfRow*j
            points.append([x,y])