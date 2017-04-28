

def spline(xl,yl,yp1=1.e30,ypn=1.e30):
	
	n = len(xl)
	yp1 = (yl[1]-yl[0])/(xl[1]-xl[0])
	ypn = (yl[n-2]-yl[n-1])/(xl[n-2]-xl[n-1])
	u = []
	y2 = []
	for i in range(0,n):
		u.append(0)
		y2.append(0)
	
	if yp1 < .99e30:
		y2[0] = -.5
		u[0] = (3./(xl[1]-xl[0]))*((yl[1]-yl[0])/(xl[1]-xl[0])-yp1)
	
	for i in range(0,n-1):
		sig=(xl[i]-xl[i-1])/(xl[i+1]-xl[i-1])
		p = sig*y2[i-1]+2.
		y2[i] = (sig-1.)/p
		u[i] = (yl[i+1]-yl[i])/(xl[i+1]-xl[i]) - (yl[i]-yl[i-1])/(xl[i]-xl[i-1])
		u[i] = (6.*u[i]/(xl[i+1]-xl[i-1])-sig*u[i-1])/p
	
	if ypn > .99e30:
		qn=un=0
	
	else:
		qn = 0.5
		un = (3./(xl[n-1]-xl[n-2]))*(ypn-(yl[n-1]-yl[n-2])/(xl[n-1]-xl[n-2]))
		
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-1]+1.)
	
	rng = range(n-1)
	rng.reverse
	for k in rng:
		y2[k]=y2[k]*y2[k+1]+u[k]
	
	return y2,xl,yl
	


def splint(xl,yl,y2l,x):
	#Note, this requires xlist is ordered with sequentially increasing x values
	k = 0
	while xl[k] < x:
		k += 1
	klo = k - 1
	khi = k
	if xl[khi] < x or xl[klo] > x:
		print xl[khi],xl[klo],x
		return 'your khi/klo do not bracket x'
	h = xl[khi]-xl[klo]
	if h == 0.0:
		return 'you used the same x value twice, dummy'
	a = (xl[khi]-x)/h
	b = (x-xl[klo])/h
	return a*yl[klo]+b*yl[khi]+((a*a*a-a)*y2l[klo]+(b*b*b-b)*y2l[khi])*(h*h)/6.