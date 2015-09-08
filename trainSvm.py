#train a scikit learn SVM onto the Wang data

from sklearn import svm
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier

from annotateOffs import *
from scipy.stats import pearsonr
import pickle

sys.path.append("libsvm-3.20/python")
from svmutil import *

def main():
    startPos = 0
    endPos = 20

    #mixedScores = parseScores("svmTraining/chari.tab")
    wangScores = parseScores("svmTraining/wang.tab")

    #random.seed(0)
    #random.shuffle(wangScores)
    #mixedScores.extend(wangScores[:len(mixedScores)])
    #mixedScores.extend(wangScores[:100])
    #mixedScores.extend(wangScores[-100:])

    vecOrder = {"A":0, "C":1, "T":2,"G":3}

    vecs = []
    scores = []
    for seq, score in wangScores:
        seq = seq[startPos:endPos]
        vec = seqToVec(seq, offsets=vecOrder)
        vecs.append(vec)
        scores.append(score)
        #vecStrs = ["%d:%d" % (i+1, x) for i,x in enumerate(vec)]
        #ofh.write("%d %s\n" % (score, " ".join(vecStrs)))
    #ofh.close()
    #print "wrote", ofh.name

    #clf = svm.SVC()
    #clf = svm.SVC(kernel="rbf", probability=True)
    #clf = svm.LinearSVC()
    #clf = svm.SVC(kernel="rbf")
    #clf = svm.SVC(kernel="rbf", probability=True)
    #clf = RandomForestRegressor()
    #clf = RandomForestClassifier()
    #clf.fit(vecs, scores)
    #print vecs
    #print scores
    #m = svm_train(scores, vecs, '-b 1 -t 0')
    m = svm_load_model("wangSabatiniSvm/wang.model")

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

    print testVecs[-1]
    print len(testVecs[-1])
    p_label, p_acc, p_val = svm_predict(testScores, testVecs, m, "-b 1")
    for seq, trueScore, predScore in zip(testSeqs, testScores, p_val):
        print seq, trueScore, predScore
    #ofh = open("svmTraining/chari.test.svmlight", "w")
    #ofh2 = open("svmTraining/chari.test.txt", "w")
    #testSeqScores = parseSvmOut("svmTraining/chariScoreCheck.tab")
    #testPreds = []
    #testScores = []
    #for seq, score in testSeqScores.iteritems():
        #seq = seq[startPos:endPos]
        ##assert(len(seq)==21)
        #vec = seqToVec(seq)
        #vecStrs = ["%d:%d" % (i+1, x) for i,x in enumerate(vec)]
        #ofh.write("%d %s\n" % (score, " ".join(vecStrs)))
        #ofh2.write("%s %f\n" % (seq, score))
        ##print clf.predict(vec), score
        #print clf.predict_log_proba(vec), score
        #testPreds.append(clf.predict_proba(vec)[0])
        ##testPreds.append(clf.predict(vec)[0])
        ##testScores.append(float(score))
    #print pearsonr(testScores, testPreds)

    #pickle.dump(clf, open("out/svm.pickle","w"))
    #print "wrote", ofh.name
    #print "wrote", ofh2.name
    print "wrote out/svm.pickle"

main()
