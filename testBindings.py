from HeavyMat import *

def main():
    #contructors and data
    myMat = heavy_matrix(6,3,True)
    myMat.print()

    myMat.cache()
    myMat.writeback()
    print(myMat)

    myMat.cache(70)
    myMat.writeback()
    print(myMat)



if __name__ == "__main__":
    main()
