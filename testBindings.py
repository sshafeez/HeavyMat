from HeavyMat import matrix

def main():
    #contructors and data
    myMat = matrix(4,3,True)
    myMat.print()
    print(myMat.grid[3][2])
    print(myMat.at(3,2))
    
    myMat.set(0,0,-9999)
    myMat.set(2,2,0)

    myMat2 = myMat
    print(myMat2)

if __name__ == "__main__":
    main()
