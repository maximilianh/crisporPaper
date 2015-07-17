# train an offtarget model based on examples
from annotateOffs import *
from sklearn import linear_model, svm

def parseOts(fname, data, y):
    for row in iterTsvRows(fname):
        vec = seqToVec(row.guideSeq[:20], row.otSeq[:20])
        data.append(vec)
        y.append(float(row.readFraction))

regr = linear_model.LinearRegression()
#regr = linear_model.Ridge(alpha=.1)
# regr = svm.SVR(kernel='sigmoid')

data = []
y = []
parseOts("hsu2013/hsuSingle.tab", data, y)

print "fitting..."
regr.fit(data, y)
#print regr.coef_
print "report:"
#coef = regr.coef_
#print "coefficients:"
#print "\n".join([str(x) for x in coef])
print "explained variance score on train data", regr.score(data, y)

testData = []
testY = []
parseOts("out/annotFiltOfftargets.tsv", testData, testY)
print "explained variance score on real data", regr.score(testData, testY)
