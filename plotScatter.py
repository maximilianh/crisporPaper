import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy
import matplotlib.backends.backend_pdf as pltBack

pdf = pltBack.PdfPages("gcCont" + '.pdf')
fig = plt.figure(figsize=(5,5),
               dpi=300, facecolor='w')
data = numpy.genfromtxt("table2.tsv", names=True)
fig = plt.figure()
plt.scatter(data["gcContent"], data["offtargetCount"], alpha=.5, marker="o", s=50)
plt.xlabel("GC content")
plt.ylabel("Number of off-targets")
fig.savefig(pdf, format = 'pdf')

pdf.close()

