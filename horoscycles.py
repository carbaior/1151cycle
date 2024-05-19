#!/usr/bin/python3

#horoscycles is a program to search for planetary cycles of the solar system of the ancients:
#Sun, Moon, Saturn, Jupiter, Mars, Venus, Mercury

#copyleft 2024 Carlos Baiget Orts (asinfreedom@gmail.com)
#https://github.com/carbaior/horoscycles

import lzma, numpy as np
from skyfield.api import Star, load, GREGORIAN_START

#Reference date: (any date is as good as any other)
day=25
month=12
year=1

#timespan of series (greater values wont give other solutions)
centuries=1 

ts = load.timescale()
ts.julian_calendar_cutoff = GREGORIAN_START

jd_ref=ts.utc(year,month,day,12)
jd_ref=int(jd_ref.tdb)
###################################################

print ("Generating list of series beginning dates...")
listanav=[]
r_ini=-1500
r_fin=1500
for i in range(r_ini, r_fin):
	t = ts.utc(i,month,day,12)
	jd_comp=int(t.tdb)
	listanav.append([int(t.tdb),i])
print ("Done.")
print ()
print ("Loading planets positions...")
data_list = []
j=0
k=0
tam=len(listanav)
a=[]
with lzma.open('planetpos.dat.lzma', 'rt') as f:
	for line in f:
		reg=[int(value) for value in line.strip().split('\t')]
		a.append(reg)
		if (reg[0]==jd_ref):
			pos_jd_ref=k
		if j<tam and reg[0]==listanav[j][0]:
			listanav[j].append(k)
			j+=1
		k+=1	
print("Done.")    
print()

serie=[]
listaseries=[]
series_long=int(100*centuries*365.25) # four centuries: 400*365.25=146100

print (f"Beginning date of the reference series: ** {year}/{month}/{day} (Julian Day {jd_ref}) **")
print ()
print (f"Series of {series_long} days long. ({centuries} centuries.)")
print ()
print (f"Search range: from {r_ini} to {r_fin}")
print ()
tot=len(listanav)
act=0

a=np.array(a)
a1 = a[pos_jd_ref:pos_jd_ref+series_long, 1:8]
for x in listanav:
	act+=1
	pc=round(((act*100)/tot),1)
	print (f"Computing series beginning on year {x[1]} ({pc}%)          ", end="\r")
	a2 = a[x[2]:x[2]+series_long, 1:8]
	b = np.abs(a2 - a1)
	b = np.where(b > 1800, 3600 - b, b)
	mean_b = np.mean(b, axis=1).reshape(-1, 1)
	m1 = np.mean(mean_b)
	m2 = np.std(mean_b)
	linea=[x[0]-jd_ref,round((m1+m2)/10,1),x[1]]
	listaseries.append(linea)
	serie=[]

listaseries=np.array(listaseries)
srtd = listaseries[np.argsort(listaseries[:, 1])]

print()
print()
print ("The TWO best cycles found are:")
print()
print("Days before/after\tYears before/after\tAverage deviation")
for elem in srtd[1:3]:
    print("{:>18}\t{:>18}\t{:>16}º".format(int(elem[0]), int(elem[0]/365.25), round(elem[1],1)))
	
#write to file	
filename = "cycles.csv"
with open(filename, "w") as f:
    for elem in listaseries:
        f.write(f"{elem[2]};{elem[1]}\n")

print()
print(f"Full output in '{filename}'")
print("(Load to spreadsheet and draw a scatter plot to compare each series affinity to the reference series)")
