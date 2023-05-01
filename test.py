
import sys, os
import math
import random

try:
	coordfile=sys.argv[1] # coordsNICO.txt
except:
	sys.exit('<coordinates file>')

#N=430597 # dimensioni del genoma # 430597
rdiff=20 # spessore dei geni
rtext=2 # spazio tra testo e geni
r=250 # raggio principale
r1= r+rdiff # raggio esterno (geni plus)
r2= r-rdiff # raggio interno (geni minus)
r3=50 # raggio nome geni in overlap
r4=3 # spazio tra nome gene in overlap e linea di connessione, e' lo stesso tra la linea di connessione e il gene
r5=30 # raggio per la curvatura nella linea di connessione
scp=0.5 # spazio tra i nomi dei geni plus in overlap 
scm=2 # spazio tra i nomi dei geni minus in overlap

pi=math.pi

def cos(x):
	return math.cos(x)

def sin(x):
	return math.sin(x)

def bp2gr(x):
	"""
	gradi=360 * x / N
	"""
	gr=((360*x)/float(N))+270
	return gr

def getColor():
	return '#'+hex(random.randint(0, 16777215))[2:].upper()

def bp2rad(x):
	"""
	radianti = gradi * pi /180
	"""
	gr=bp2gr(x) #((360*x)/float(N))+270
	return (gr*pi)/180

def gr2rad(x):
	return (x*pi)/180

def getPos(r,ar): # r: raggio ar: angolo radianti 
	x=r*cos(ar)
	y=r*sin(ar)
	return str(x),str(y) 

def getPath(start,end,r,rx,color,id): # #e2001a
	# M ax,ay A r,r 0 0,1 bx,by L cx,cy A rx,rx 0 0,0 dx,dy Z
	addcolor='fill:%s' %(color)
	style=['stroke:none', 'stroke-dasharray:none', 'stroke-opacity:1']
	style.append(addcolor)
	a=getPos(r, bp2rad(start))
	b=getPos(r,bp2rad(end))
	c=getPos(rx, bp2rad(end))
	d=getPos(rx, bp2rad(start))
	p="M %s,%s A %s,%s 0 0,1 %s,%s L %s,%s A %s,%s 0 0,0 %s,%s Z" %(a[0],a[1],str(r),str(r),b[0],b[1],c[0],c[1],str(rx),str(rx),d[0],d[1])	
	path='<path id="%s" d="%s" style="%s" />' %(id,p,';'.join(style))
	return path

def getText(point,r,testo,direction):
	style='font-variant:normal;font-weight:normal;font-size:9.5px;font-family:Helvetica;-inkscape-font-specification:Helvetica;writing-mode:lr-tb;fill:#000000;fill-opacity:1;fill-rule:nonzero;stroke:none'
	recm=bp2gr(point)
	recu=recm-0.5
	recd=recm+0.5
	apos=bp2gr(point) # angolo in gradi
	if abs((360+270)-apos) < 180:
		pos= getPos(r,gr2rad(recu))
		style=style+';text-anchor:end'
		apos+=180
	else: pos= getPos(r,gr2rad(recd))
	t='<text x="%s" y="%s" transform="rotate(%s %s,%s)" style="%s" direction="%s" >%s</text>' %(pos[0],pos[1],str(apos),pos[0],pos[1],style,direction,testo)
	return t

def getText2(start,end,r,testo,direction):
	style='font-variant:normal;font-weight:normal;font-size:9.5px;font-family:Helvetica;-inkscape-font-specification:Helvetica;writing-mode:lr-tb;fill:#000000;fill-opacity:1;fill-rule:nonzero;stroke:none'
	recu=start
	recd=end
	apos=start
	if abs((360+270)-apos) < 180:
		pos= getPos(r,gr2rad(recu))
		style=style+';text-anchor:end'
		apos+=180
	else: pos= getPos(r,gr2rad(recd))
	t='<text x="%s" y="%s" transform="rotate(%s %s,%s)" style="%s" direction="%s" >%s</text>' %(pos[0],pos[1],str(apos),pos[0],pos[1],style,direction,testo)
	return t

def readCoords(infile):
	plus,minus=[],[]
	gname='' #nome del genoma
	size=0
	f=open(infile)
	for i in f:
		if i.startswith('>'):
			gname=(i.strip()).replace('>','')
			continue
		if i.startswith('#'):
			size=int((i.strip()).replace('#',''))*1.0
			continue		
		l=(i.strip()).split('\t')
		if l[-1]=='-': minus.append((int(l[0]),int(l[1]),l[2]))
		elif l[-1]=='+': plus.append((int(l[0]),int(l[1]),l[2]))
	f.close()
	plus.sort()
	minus.sort()
	return plus,minus,gname,size

def makeCluster(allcoord):
	cluster=[]
	remaining=[]
	c1=allcoord[0][0]
	c2=allcoord[0][1]
	for i in range(len(allcoord)):
		if allcoord[i]!=(c1,c2):
			if c1<=allcoord[i][0]<=c2:
				cluster.append(allcoord[i])
				if allcoord[i][1]>c2:
					c2=allcoord[i][1]
			else:
				remaining.append(allcoord[i])
		else:
			cluster.append((c1,c2))
	return remaining,cluster

def getClusters(interval):
	clusters=[]
	interval.sort()
	while len(interval)!=0:
		interval,cluster=makeCluster(interval)
		clusters.append(cluster)
	return clusters

def writeName(name="mygenome",size="unknown"):
	t1='<text x="0" y="0" text-anchor="middle" style="font-style:oblique;font-variant:normal;font-weight:bold;font-size:20px;font-family:Helvetica;-inkscape-font-specification:Helvetica-BoldOblique;writing-mode:lr-tb;fill:#000000;fill-opacity:1;fill-rule:nonzero;stroke:none">%s</text>' %(name)
	#t2='<text x="0" y="25" text-anchor="middle" style="font-variant:normal;font-weight:normal;font-size:18px;font-family:Helvetica;-inkscape-font-specification:Helvetica;writing-mode:lr-tb;fill:#000000;fill-opacity:1;fill-rule:nonzero;stroke:none">type</text>'
	t3='<text x="0" y="25" text-anchor="middle" style="font-variant:normal;font-weight:normal;font-size:18px;font-family:Helvetica;-inkscape-font-specification:Helvetica;writing-mode:lr-tb;fill:#000000;fill-opacity:1;fill-rule:nonzero;stroke:none">%s bp</text>' %(size.replace('.',','))
	return t1+t3

plus,minus,gname,N=readCoords(coordfile)

SVG=['<svg version="1.1" baseProfile="full" width="800" height="800" viewBox="0 0 800 800" preserveAspectratio="xMinYMin" xmlns="http://www.w3.org/2000/svg">']
SVG.append('<g transform="translate(400,400)">')

plusAcoords,minusAcoords=[],[]
for i in plus:
	path=getPath(i[0],i[1],r,r1,getColor(),i[2]) #'#e2001a'
	SVG.append(path)
	midpoint=(i[0]+i[1])/2.0
	recm=bp2gr(midpoint)
	recu=recm-0.7
	recd=recm+0.7
	plusAcoords.append((recu,recd,i,recm))

for i in minus:
	path=getPath(i[0],i[1],r2,r,getColor(),i[2])
	SVG.append(path)
	midpoint=(i[0]+i[1])/2.0
	recm=bp2gr(midpoint)
	recu=recm-0.7
	recd=recm+0.7
	minusAcoords.append((recu,recd,i,recm))
	
plusClusters=getClusters(plusAcoords)
minusClusters=getClusters(minusAcoords)

for i in plusClusters:
	if len(i)>1:
		sd=[j[1]-j[0] for j in i]
		areal=i[-1][1]-i[0][0]
		sc = ((sum(sd)+(scp*(len(sd)-1)))-areal)/2
		start=i[0][0]-sc
		for j in range(len(i)):
			text=getText2(start,start+sd[j],r1+r3,i[j][2][2],"ltr")
			SVG.append(text)
			midt=((start+sd[j])+start)/2.0
			midg=i[j][3]
			a=getPos(r1+r4,gr2rad(midg))
			b=getPos(r1+(r3-r4),gr2rad(midt))
			c=getPos(r1+r5,gr2rad(midg))
			pr='<path d="M %s %s S %s %s, %s %s" fill="transparent" stroke="black"/>' %(a[0],a[1],c[0],c[1],b[0],b[1])
			SVG.append(pr)
			start+=sd[j]+scp
	else:
		text=getText((i[0][2][0]+i[0][2][1])/2,r1+rtext,i[0][2][2],"ltr")
		SVG.append(text)

for i in minusClusters:
	if len(i)>1:
		sd=[j[1]-j[0] for j in i]
		areal=i[-1][1]-i[0][0]
		sc = ((sum(sd)+(scm*(len(sd)-1)))-areal)/2
		start=i[0][0]-sc
		for j in range(len(i)):
			text=getText2(start,start+sd[j],r2-r3,i[j][2][2],"rtl")
			SVG.append(text)
			midt=((start+sd[j])+start)/2.0
			midg=i[j][3]
			a=getPos(r2-r4,gr2rad(midg))
			b=getPos(r2-(r3-r4),gr2rad(midt))
			c=getPos(r2-r5,gr2rad(midg))
			pr='<path d="M %s %s S %s %s, %s %s" fill="transparent" stroke="black"/>' %(a[0],a[1],c[0],c[1],b[0],b[1])
			SVG.append(pr)			
			start+=sd[j]+scm
	else:
		text=getText((i[0][2][0]+i[0][2][1])/2,r2-rtext,i[0][2][2],"rtl")
		SVG.append(text)


SVG.append('<circle cx="0" cy="0" r="%i" style="stroke: black; stroke-width:2; fill:none;"/>' %(r))
genomeName=writeName(gname,str(int(N)))
SVG.append(genomeName)
##
SVG.append('</g>')
SVG.append('</svg>')
  
print '\n'.join(SVG)
