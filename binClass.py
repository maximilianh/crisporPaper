# test the scores and rules on various datasets using a binary classification
# the basis is that we consider the top25 of any dataset as the target
# all scores use 75-percentiles as cutoffs to predict the targets

import operator

import numpy as np
from sklearn import svm
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree
from sklearn import dummy
from sklearn import cross_validation
from sklearn.cross_validation import LeaveOneOut, KFold
from sklearn.metrics import recall_score, precision_score, f1_score

from annotateOffs import *
from scipy.stats import pearsonr
import pickle
#import pydot
from sklearn.externals.six import StringIO

# start position of sequence to take into account
startPos = 14
# end position to take into account for classifier
endPos = 20

# datasets to evaluate with loo cross validation
evalSets = ["farboud2015", "ren2015", "gagnon2014", "varshney2015"]

# datasets to train dec tree classifiers on
#trainSets = ["farboud2015", "ren2015", "gagnon2014", "varshney2015"]
trainSets = ["gagnon2014"]

# datasets to apply dec tree classifier on
testSets = ["schoenig", "xu2015TrainHl60", "chari2015Train", "eschstruth", "varshney2015", "ren2015", "farboud2015", "doench2014-Hs", "housden2015", "morenoMateos2015", "alenaAll"]
#testSets = ["eschstruth"]
# removed:"museumIC50", 
# removed: "xu2015FOX-AR", "xu2015AAVS1", 'chari2015Valid_293T', 
testSets.extend(trainSets)
# training datasets
#inData = ["ren2015", "farboud2015"]
#inData = ["ren2015"]
inData = ["farboud2015"]
#inData = ["ren2015"]
#inData = ["varshney2015"]
#inData = ["gagnon2014"]
#inData = ["chari2015Train"]

def parseIn(inData):
    allSeqs, allScores = [], []
    for dataset in inData:
        seqs, scores = parseSeqScores(dataset)
        if dataset.startswith("xu2015Train"):
            scores = [-x for x in scores]
        scores = useRanks(scores, doQuart=True)
        allSeqs.extend(seqs)
        allScores.extend(scores)
    return np.array(allSeqs), np.array(allScores)

def scorePreds(y, yPred):
    rec= recall_score(y, yPred)
    prec= precision_score(y, yPred)
    f1= f1_score(y, yPred)
    posCount = len([y for y in yPred if y==1.0])
    return posCount, rec, prec, f1

def highFinalGc(vec):
    cgCount = 0
    vecLen = len(vec)
    for i in range(0,6):
        cgCount += vec[vecLen-((i*4)+1)] + vec[vecLen-((i*4)+2)]

    if cgCount >=4:
        return 1
    return 0

def parseAndBinarize(inData):
    seqs, scores = parseIn([inData])
    vecs = seqsToVecs(seqs, startPos=startPos, endPos=endPos)
    vecs = np.array(vecs)

    convScores = {1:0, 2:0, 3:0, 4:1}
    # these only have three levels, so consider everything in the top two quartiles as good
    if inData=="schoenig" or inData=="concordet2":
        convScores = {1:0, 2:0, 3:1, 4:1}
    if inData=="eschstruth":
        convScores = {1:0, 2:0, 3:0, 4:1}
    # museum scores are inverted, so use the lower two quartiles as "top"
    if inData.startswith("museum"):
        convScores = {1:1, 2:1, 3:0, 4:0}

    labels = np.array([convScores[s] for s in scores])
    return seqs, vecs, labels

def rulePredScores(vecs, labels):
    ruleGGPreds = [] # ends with GG ?
    ruleGCPreds = [] # no of GC in final 6bp > 4?
    for vec in vecs:
        ruleGCPreds.append(highFinalGc(vec))

        if vec[-2]==1 and vec[-6]==1:
            ruleGGVal = 1
        else:
            ruleGGVal = 0
        ruleGGPreds.append(ruleGGVal)

    ruleGGRec, ruleGGPrec, ruleGGF1 = scorePreds(labels, ruleGGPreds)
    ruleGCRec, ruleGCPrec, ruleGCF1 = scorePreds(labels, ruleGCPreds)
    return ruleGGRec, ruleGGPrec, ruleGGF1, ruleGCRec, ruleGCPrec, ruleGCF1

def evalDatasets(trainSets):
    for inData in evalSets:
        seqs, vecs, labels = parseAndBinarize(inData)
        ggRec, ggPrec, ggF1, gcRec, gcPrec, gcF1 = rulePredScores(vecs, labels)

        cVal = LeaveOneOut(len(labels))
        cValPreds = []
        cValTests = []
        for train, test in cVal:
            X_train, X_test = vecs[train], vecs[test]
            y_train, y_test  = labels[train], labels[test]
            #clf = svm.SVC(kernel="linear", probability=True)
            clf = tree.DecisionTreeClassifier(min_samples_leaf=4, max_depth=4)
            #clf = RandomForestClassifier()
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            cValPreds.append(y_pred[0])
            cValTests.append(y_test[0])

        clfRec, clfPrec, clfF1 = scorePreds(cValTests, cValPreds)


        row = [ggRec, ggPrec, gcRec, gcPrec, clfRec, clfPrec, clfF1]
        row = ["%0.2f" % x for x in row]
        row.insert(0, inData)
        print "\t".join(row)
    
        classifiers.append( ("DecTree_"+inData, clf) )

        #clf = dummy.DummyClassifier()
        #clf.fit(vecs, labels)
        #classifiers.append( ("Dummy_"+inData, clf) )

    return classifiers

def trainClassifiers(trainSets):
    " train classifiers on trainSets and return as a list (name, clfObj) "
    classifiers = {}
    for inData in trainSets:
        seqs, vecs, labels = parseAndBinarize(inData)
        ggRec, ggPrec, ggF1, gcRec, gcPrec, gcF1 = rulePredScores(vecs, labels)

        clf = tree.DecisionTreeClassifier(min_samples_leaf=4, max_depth=4)
        clf.fit(vecs, labels)
        classifiers["DecTree_"+inData] = clf

        #clf = dummy.DummyClassifier()
        #clf.fit(vecs, labels)
        #classifiers["Dummy_"+inData] = clf

    return classifiers

def writeTree(clf, outFname):
    " write a decision tree to a pdf file "
    f = open("out/temp.dot", "w")
    names = []
    for i in range(startPos, endPos):
        for c in "ACGT":
            names.append("%d%s" % (i, c))
    tree.export_graphviz(clf, out_file=f, feature_names=names)
    #graph = pydot.graph_from_dot_data(dot_data.getvalue())
    #graph.write_pdf("temp.pdf")
    cmd = "dot out/temp.dot -Tpng -o %s" % outFname
    print cmd
    assert(os.system(cmd)==0)
    print "wrote tree to %s" % outFname

def map23To34():
    " return a map from 23mer to 34mer "
    shortToLong = {}
    for fname in glob.glob("effData/*.ext.tab"):
        #print fname
        for row in iterTsvRows(fname):
            long = row.extSeq
            short = row.seq
            shortToLong[short] = long
    return shortToLong

def evalAllScores_takeBestX(datasetName, seqs, labels):
    """ create posCount/rec/prec/f1 for all main scores, defining positive as 'best X'
    with X being the number of TPs in the target
    """
    #shortToLong = map23To34()
    #longSeqs = [shortToLong[s] for s in seqs]
    #seqScores = calcEffScores(longSeqs, skipOof=True)
    seqScores, freqs = parseEffScores(datasetName)
    #seqScoreDict = dict(zip(seqs, seqScores))

    posCount = len([x for x in labels if x==1]) # number of positives

    typeEval = {}
    #scoreTypes = ('doench', 'svm', 'chariRank', 'ssc')
    scoreTypes = getScoreTypes()
    scoreTypes.append("finalGc6")
    scoreTypes.append("finalGg")

    for scoreType in scoreTypes:
        # first create a list (seq, score) 
        seqScoreList = []
        for seq in seqs:
            scoreDict = seqScores[str(seq)]
            score = scoreDict[scoreType]
            seqScoreList.append ( (score, seq) )
        seqScoreList.sort(reverse=True)

        # get the best sequences, up to posCount
        posSeqs = set()
        for i, (score, seq) in enumerate(seqScoreList):
            if i==posCount:
                break
            posSeqs.add(seq)

        # now create the labels 1/0 and save them
        predLabels = []
        for seq in seqs:
            if seq in posSeqs:
                predLabels.append(1)
            else:
                predLabels.append(0)

        typeEval[scoreType] = scorePreds(labels, predLabels)
    return typeEval

def evalAllScores_75Perc(datasetName, seqs, labels):
    " create rec/prec/f1 for all main scores, defining TP as 'over the 75 percentile' "
    # run scoreCutoffs.py to get these
    #cutoffs = {'doench': 0.32398638678651903, 'wangOrig': 0.77276957430125004, 'chariRank': 79.365247038449994, 'ssc': 0.41655549999999997, "crisprScan" :0.61, "fusi":0.62}
    #cutoffs = {'oof': 68.0, 'wang': 0.7466315, 'drsc': 6.83643, 'finalGg': 0.0, 'wangOrig': 0.77927021985800005, 'chariRaw': 0.53029115000000004, 'finalGc6': 1.0, 'chariRank': 79.0, 'mh': 5447.0, 'doench': 32.0, 'crisprScan': 60.0, 'ssc': 0.420796, 'fusi': 0.62621093132700001}
    cutoffs = {'oof': 68.0, 'wang': 0.76407499999999995, 'drsc': 6.8228999999999997, 'finalGg': 0.0, 'wangOrig': 0.78028833768625006, 'chariRaw': 0.52830615000000003, 'finalGc6': 1.0, 'chariRank': 79.0, 'mh': 5424.25, 'doench': 33.0, 'crisprScan': 60.0, 'ssc': 0.41735925000000001, 'fusi': 0.62641138794899998}


    #shortToLong = map23To34()
    #longSeqs = [shortToLong[s] for s in seqs]
    #seqScores = calcEffScores(longSeqs, skipOof=True)
    seqScores, freqs = parseEffScores(datasetName)

    typeEval = {}
    #for scoreType, cutoff in cutoffs.items():
    scoreTypes = getScoreTypes()
    scoreTypes.append("finalGc6")
    scoreTypes.append("finalGg")
    for scoreType in scoreTypes:
        cutoff = cutoffs[scoreType]
        # the two rules already have binary values 0 or 1
        if scoreType.startswith("final"):
            cutoff = 0.5

        predLabels = []
        for seq in seqs:
            scoreDict = seqScores[seq]
            score = scoreDict[scoreType]
            if score > cutoff:
                predLabels.append(1)
            else:
                predLabels.append(0)
        typeEval[scoreType] = scorePreds(labels, predLabels)
    return typeEval

def evalClassifiers(testSets):
    """ get recall, precision, f1 for a list of classifiers. Also get these
    metrics for for efficiency scores (doench, etc) and the two heuristics
    """
    rows = []
    for dataset in testSets:
        seqs, vecs, labels = parseAndBinarize(dataset)
        posCount = len([x for x in labels if x==1]) # number of positives
        dataSize = len(seqs)

        #for clfName, clf in classifiers.iteritems():
            #y_pred = clf.predict(vecs)
            #clfRec, clfPrec, clfF1 = scorePreds(labels, y_pred)
            #row = [clfName, dataset, dataSize, posCount, clfRec, clfPrec, clfF1]
            #rows.append(row)

        #ggRec, ggPrec, ggF1, gcRec, gcPrec, gcF1 = rulePredScores(vecs, labels)
        #row = ["finalGg", dataset, dataSize, posCount, ggRec, ggPrec, ggF1]
        #rows.append(row)

        #row = ["finalGc6", dataset, dataSize, posCount, gcRec, gcPrec, gcF1]
        #rows.append(row)

        # add the metrics for all efficiency scores
        scoreTypeEvalsBestX = evalAllScores_takeBestX(dataset, seqs, labels)
        scoreTypeEvalsGt75 = evalAllScores_75Perc(dataset, seqs, labels)
        for scoreType, (predCount, rec, prec, f1) in scoreTypeEvalsGt75.iteritems():
            bestXPredCount, bestXRec, bestXPrec, bestXF1 = scoreTypeEvalsBestX[scoreType]
            #assert(bestXPredCount==predCount)
            assert(bestXRec==bestXPrec)
            row = [scoreType, dataset, dataSize, posCount, predCount, rec, prec, f1, bestXPredCount, bestXPrec]
            rows.append(row)
    return rows

def main():
    random.seed(0)
    ofh = open("out/binClassMetrics.tsv", "w")

    #classifiers = trainClassifiers(trainSets)
    #classifiers = {}

    headers = ["classifierName", "dataset", "size", "posCount", "predCount", "recall", "precision", "f1", "bestXPredCount", "bestXAcc"]
    ofh.write ( "\t".join(headers)+"\n")

    rows = evalClassifiers(testSets)
    rows.sort(key=operator.itemgetter(-1), reverse=True)

    for row in rows:
        row[2:] = ["%0.2f" % x for x in row[2:]]
        ofh.write( "\t".join(row)+"\n")

    print "wrote to %s" % ofh.name
    #writeTree(classifiers["DecTree_gagnon2014"], "out/%s.png" % "gagnon")

    #loo = LeaveOneOut(len(labels), vecs)
    #for train, test in loo:
        #print "train", train, "test", test
        #print scores[train], scores[test]
        #print scores[np.ix_(train)]
        #X_train, X_test = vecs[train], vecs[test]
        #y_train, y_test  = labels[train], labels[test]
        #clf = tree.DecisionTreeClassifier()
        #clf.fit(X_train, y_train)
        #y_pred = clf.predict(X_test)
        #print y_test, y_pred



    #clf = svm.LinearSVR()
    #clf = svm.SVC(kernel="rbf", probability=True)
    #clf = RandomForestClassifier()
    #scores = cross_validation.cross_val_score(clf, vecs, scores, cv=10, scoring='precision')
    #scores = cross_validation.StratifiedKFold(clf, vecs, scores, cv=10, scoring='precision')
    #print scores

    #clf = svm.SVC()
    #clf = svm.SVR(kernel="rbf")
    #clf = RandomForestRegressor()
    #clf = RandomForestClassifier()
    #clf.fit(vecs, scores)


main()
