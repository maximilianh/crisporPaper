# use 5-fold cross validation to quantify the "predictability" of modification frequency
import numpy as np
from sklearn import cross_validation
from sklearn import datasets
from sklearn import svm
from sklearn import linear_model
from annotateOffs import *
from scipy.stats import linregress, pearsonr, spearmanr, mannwhitneyu, rankdata
import sys
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.mixture import *
from sklearn.linear_model import LinearRegression

def parseTab(fname):
    xList = []
    yList = []
    seqs = []
    for row in iterTsvRows(open(fname)):
        seq = row.seq[17:20]
        seqs.append(seq)
        score = float(row.modFreq)
        #vec = seqToVec(seq, offsets={"G":0, "other":1})
        vec = seqToVec(seq)
        xList.append(vec)
        yList.append(score)
    return xList, yList

def splitData(xList, yList, ratio):
    " given two lists, split into four lists, ratio indicates size of first "
    xTest, yTest = [], []
    xTrain, yTrain = [], []
    assert(len(xList)==len(yList))
    for i in range(0, len(xList)):
        if random.random()<ratio:
            xTrain.append(xList[i])
            yTrain.append(yList[i])
        else:
            xTest.append(xList[i])
            yTest.append(yList[i])
    return xTrain, yTrain, xTest, yTest

def parseMany(mask, names, maxSize=None):
    " parse many tab files and return as merged lists "
    xList, yList = [], []
    for name in names:
        xAdd, yAdd = parseTab(mask % name)
        yAdd = useRanks(yAdd, doQuart=True)
        if maxSize is not None:
            xAdd, yAdd = shuffleOrder(xAdd, yAdd)
            xAdd = xAdd[:maxSize]
            yAdd = yAdd[:maxSize]
        xList.extend(xAdd)
        yList.extend(yAdd)

    return list(xList), list(yList)

def shuffleOrder(xList, yList):
    mixList = zip(xList, yList)
    random.shuffle(mixList)
    return zip(*mixList)

models = [ \
    #("DPGMM", DPGMM(n_components=5, covariance_type='diag', alpha=100, n_iter=100)), \
    #("DPGMM", DPGMM(n_components=5, covariance_type='diag', alpha=100, n_iter=100)), \
    #("GMM", GMM(n_components=10, covariance_type="tied", init_params='wc', n_iter=20)),
    ("lasso", linear_model.Lasso(alpha=0.01)), \
    ("SVM", svm.SVR(kernel="poly")), \
    ("RF-regression", RandomForestRegressor()), \
    ("Ridge-Regression", linear_model.Ridge(alpha=0.1)) \
    ]
models = [ \
    #("DPGMM", DPGMM(n_components=5, covariance_type='diag', alpha=100, n_iter=100)), \
    #("DPGMM", DPGMM(n_components=5, covariance_type='diag', alpha=100, n_iter=100)), \
    #("GMM", GMM(n_components=10, covariance_type="tied", init_params='wc', n_iter=20)),
    ("lasso", linear_model.Lasso(alpha=0.01)) \
    ]

#xList, yList = parseTab("effData/varshney2015.tab")
#xList, yList = parseTab("effData/gagnon2014.tab")
datasets = ["meta", "doench2014-Hs", "xu2015Train", "chari2015Train", "varshney2015", "gagnon2014"]
#, "schoenig", "museumT7", "farboud2015", "ren2015"]
for dataset in datasets:
    if dataset=="meta":
        xList, yList = parseMany("effData/%s.ext.tab", ["doench2014-Hs", "xu2015Train", "chari2015Train"], 1200)
        xTrain, yTrain, metaTestX, metaTestY = splitData(xList, yList, 0.2)
        #print "xl", xList
        #print "yl", yList
    else:
        xList, yList = parseTab("effData/%s.ext.tab" % dataset)
        if dataset=="doench2014-Hs":
            xAdd, yAdd = parseTab("effData/doench2014-Mm.tab")
            xList.extend(xAdd)
            yList.extend(yAdd)
    xList, yList = shuffleOrder(xList, yList)
    #xList = xList[:1200]
    #yList = yList[:1200]
    print "dataset %s, size %d=%d:" % (dataset,  len(xList), len(yList))
    if dataset=="meta":
        #metaRegr = svm.SVR(kernel="poly")
        #metaRegr = RandomForestRegressor()
        metaRegr = RandomForestClassifier()
        #metaRegr = linear_model.LassoCV()
        #metaRegr = svm.LinearSVR()
        #metaRegr = linear_model.Lasso(alpha=0.01)
        #metaRegr = linear_model.Lasso(alpha=0.01)
        #metaRegr = linear_model.BayesianRidge()
        # metaRegr = linear_model.LogisticRegressionCV() too slow
        #metaRegr = linear_model.Ridge (alpha = .05)
        #metaRegr = linear_model.RANSACRegressor(linear_model.LinearRegression())
        #metaRegr = RandomForestRegressor()
        #yTrain = useRanks(yTrain, doQuart=True)
        # data is already quartiled, see parseMany()
        metaRegr.fit(xTrain, yTrain)

    xTrain, yTrain, xTest, yTest = splitData(xList, yList, 0.2)
    for modelName, regr in models:
        rVals = []
        for i in range(0, 5):
            #regr = linear_model.Lasso(alpha=0.01)
            #regr = svm.SVR()
            #regr = RandomForestRegressor()
            #regr = linear_model.Ridge(alpha=0.1) # -4.5
            regr.fit(xTrain, yTrain)
            testPreds = [regr.predict(x)[0] for x in xTest]
            pearR, pearP = pearsonr(testPreds, yTest)
            #print "on 1/5 of input: Pearson R %0.3f (P %0.3f)" % (pearR, pearP)
            rVals.append(pearR)
        print "%s: Avg Pearson R %0.3f" % (modelName, np.mean(pearR))


    #for x, y in zip(xTest, yTest):
        #testPreds.append(regr.predict(x)[0])
#for seq, x in zip(seqs, xList, yList):
    #print seq, x, y
#x_train = xList[:-20]
#x_test = xList[-20:]
#y_train = yList[:-20]
#y_test = yList[-20:]
#diabetes = datasets.load_diabetes()
#diabetes_X_train = diabetes.data[:-20]
#diabetes_X_test  = diabetes.data[-20:]
#diabetes_y_train = diabetes.target[:-20]
#diabetes_y_test  = diabetes.target[-20:]

#iris = datasets.load_iris()
#print iris.data
#print iris.target

#regr = linear_model.LinearRegression() # -4.73
#regr = linear_model.Ridge(alpha=0.1) # -4.5
#regr = linear_model.Lasso(alpha=0.01)
#regr = linear_model.Lasso(alpha=0.01)
#regr = RandomForestRegressor()
#scores = cross_validation.cross_val_score(regr, xList, yList, cv=5, scoring="mean_absolute_error")
#scores = cross_validation.cross_val_score(regr, xList, yList, cv=5, scoring="r2")
#print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

# best score: -4.5
#alphas = np.logspace(-4, -1, 6)
#alphas = [0.0, 0.001, 0.01, 0.1, 0.25, 0.5, 1.0]
#for alpha in alphas:
    #regr = linear_model.Lasso(alpha=alpha)
    #scores = cross_validation.cross_val_score(regr, xList, yList, cv=5)
    #print ("alpha %f" % alpha)
    #print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

#regr = svm.SVR()
#regr = svm.SVC(kernel='linear', C=1)
#regr.fit(xTrain, yTrain)
#coef = regr.coef_
#printCoef(coef)
#print "coeff", (coef)
#print "mean square error", np.mean((regr.predict(xList)-yList)**2)
#print "expl variance score", regr.score(xList, yList) 
#print "mean error", np.mean((regr.predict(xTest)-yTest))
#print "expl variance score", regr.score(xTest, yTest)
#testPreds = []
#for x, y in zip(xTest, yTest):
    #testPreds.append(regr.predict(x)[0])
#pearR, pearP = pearsonr(testPreds, yTest)
#print "1/2 of input: Pearson R %0.3f (P %0.3f)" % (pearR, pearP)

testSets = ["mixTest", "doench2014-Hs", "xu2015Train", "chari2015Train", "varshney2015", "gagnon2014", "schoenig", "museumT7", "farboud2015", "ren2015", "xu2015"]

print
print "application to full datasets"
for testSet in testSets:
    if testSet=="mixTest":
        xTest, trueY = metaTestX, metaTestY
    else:
        xTest, trueY = parseTab("effData/%s.tab" % testSet)
        trueY = useRanks(trueY, doQuart=True)

    testPreds = []
    for x, y in zip(xTest, trueY):
        testPreds.append(metaRegr.predict(x)[0])
    bestIdx = [i for i in range(0, len(trueY)) if testPreds[i]==4]
    if len(bestIdx)==0:
        print testSet, "no single 'high' prediction"
        continue
    okCount = [trueY[i]==4 for i in bestIdx].count(True)
    print testSet, "len", len(bestIdx), "ok", okCount, "perc", float(okCount)/len(bestIdx)
