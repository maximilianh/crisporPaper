# test the Wang implementation

from sklearn import svm
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier

from annotateOffs import *
from scipy.stats import pearsonr
import pickle

sys.path.append("libsvm-260/python")
from svm import *

#vecOrder = {"A":0, "C":1, "T":2,"G":3}
vecOrder = {"A":0, "C":1, "T":2,"G":3}

def main():
    startPos = 0
    endPos = 20

    m = svm_model("wangSabatiniSvm/wang.model")

    ofh = open("svmTraining/wang.comparison.txt", "w")
    testSeqScores = parseSvmOut("wangSabatiniSvm/output.txt")


    testPreds = []
    testVecs = []
    testScores = []
    testSeqs = []
    for seq, score in testSeqScores.iteritems():
        vec = seqToVec(seq, offsets=vecOrder)
        testVecs.append(vec)
        testScores.append(score)
        #predScore = clf.predict(vec)
        #testPreds.append(predScore)
        testSeqs.append(seq)
        #print seq, score, predScore
        probs = m.predict_probability(vec)
        print seq, score, probs, vec


main()
