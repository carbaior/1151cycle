#!/usr/bin/python3

#horoscycles is a program to search for planetary cycles of the solar system of the ancients:
#Sun, Moon, Saturn, Jupiter, Mars, Venus, Mercury
#copyleft 2024 Carlos Baiget Orts (asinfreedom@gmail.com)
#https://github.com/carbaior/horoscycles


import lzma, numpy as np
from skyfield.api import Star, load, GREGORIAN_START
from joblib import Parallel, delayed
	
def desvest(l):
    global x,pos_jd_ref
    i = x[2] + l  # Día juliano del día a comparar
    j = pos_jd_ref + l  # Día juliano del día de referencia
    angs = np.mod(np.abs(a[i] - a[j]), 3600)
    return np.std(np.where(angs > 1800, 3600 - angs, angs))

print ()    
print ("Loading planets positions...")
data_list = []
with lzma.open('planetpos.dat.lzma', 'rt') as f:
	for line in f:
	    data_list.append([int(value) for value in line.strip().split('\t')])
plposdat = np.array(data_list)
print("Done.")    
print()

print ("Generating list of series beginning dates...")
ts = load.timescale()
ts.julian_calendar_cutoff = GREGORIAN_START
data_list=[]
r_ini=-1500
r_fin=1500
for i in range(r_ini, r_fin):
	t = ts.utc(i,12,25,12)
	jd_comp=int(t.tdb)
	data_list.append([jd_comp,i,np.searchsorted(plposdat[:, 0], jd_comp)])
print ("Done.")
print ()
listanav = np.array(data_list)
del data_list

#Reference date: (any date is as good as any other)
day=25
month=12
year=1
###################################################
jd_ref=ts.utc(year,month,day,12)
jd_ref=int(jd_ref.tdb)
pos_jd_ref = np.searchsorted(plposdat[:, 0], jd_ref)

a = plposdat[:, 1:8]
serie=[]
listaseries=[]
centuries=100 #timespan of series (greater values wont give other solutions)
series_long=int(centuries*365.25) # four centuries: 400*365.25=146100

print (f"Beginning date of the reference series: ** {year}/{month}/{day} (Julian Day {jd_ref}) **")
print ()
print (f"Series of {series_long} days long. ({centuries} centuries.)")
print ()
print (f"Search range: from {r_ini} to {r_fin}")
print ()
tot=len(listanav)
act=0
for x in listanav:
	act+=1
	pc=round(((act*100)/tot),1)
	print (f"Computing series beggining on year {x[1]} ({pc}%)          ", end="\r")
	serie = Parallel(n_jobs=-1)(delayed(desvest)(l) for l in range(series_long))
	avgserie = np.mean(serie)/10
	stdserie = np.std(serie)/10
	linea=[x[0]-jd_ref,round(avgserie,1)+round(stdserie,1),x[1]]
	listaseries.append(linea)

listaseries=np.array(listaseries)

srtd = listaseries[np.argsort(listaseries[:, 1])]

print()
print()
print ("The THREE best cycles found are:")
print()
print("Days before/after\tYears before/after\tAverage deviation")
for elem in srtd[:3]:
    print("{:>18}\t{:>18}\t{:>16}º".format(int(elem[0]), int(elem[0]/365.25), elem[1]))
	
#write to file	
filename = "cycles.csv"
with open(filename, "w") as f:
    for elem in listaseries:
        f.write(f"{elem[2]};{elem[1]}\n")

print()
print(f"Full output in '{filename}'")
print("(Load to spreadsheet and draw a scatter plot to compare each series affinity to the reference series)")
