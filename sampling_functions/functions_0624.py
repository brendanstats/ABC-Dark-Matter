import numpy as np
import matplotlib.pyplot as plt
import random as rand
from scipy import integrate
import scipy.optimize
import time

#model function; param = [a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax]
def f(x, vr, vt, param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5]):

	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param

        rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)
        Pr = Ps * ( 1 - (np.log(1+x))/x )
        J  = x*rs * vt           #J = v * r*sin(theta)
        E  = (vt*vt + vr*vr)/2.0 + Pr # v*v/2 + Pr

        Ec   = Ec * Ps
        xlim = rlim / rs #turn rlim in unit of rs
	Plim = Ps * ( 1 - (np.log(1+xlim))/xlim ) #0.45*Ps
        Jb   = Jb * rs * (Ps**0.5) #*0.086

        if b <= 0:
                gJ = 1.0/(1 + (J/Jb)**-b)
        else:
                gJ = 1 + (J/Jb)**b

        N  = 1.0*10**1
        if E < Plim and E >= 0:
                hE = N*(E**a) * ((E**q + Ec**q)**(d/q)) * ((Plim - E)**e)
        else:
                hE = 0.0

        #return hE * gJ
        return hE


def hE(x, v, param):
	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
	rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)
        Pr = Ps * ( 1 - (np.log(1+x))/x )

	E  = (v*v)/2.0 + Pr
	Ec = Ec * Ps
	xlim = rlim / rs #turn rlim in unit of rs
        Plim = Ps * ( 1 - (np.log(1+xlim))/xlim ) #0.45*Ps

	N = 1.
	if E < Plim and E >= 0:
                h = N*(E**a) * ((E**q + Ec**q)**(d/q)) * ((Plim - E)**e)
        else:
                h = 0.0

	return h

#model probility function
def fprob(x, vr, vt, param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5]):
	return x*x*vt*f(x, vr, vt, param)



#step function g(x,vr,vt) use for sampling
def fbdvalues(param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], steps = 10):
	xi  = 10**-7
	vri = 0
	vti = 0
	
	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param

	#find the maximum point
	rs    = rmax/2.16       # rmax=2.16*rs
        fmax1 = findmax(param) * 1.1
        xm,vrm, vtm = scipy.optimize.fmin(lambda (x,vr,vt): -fprob(x,vr,vt, param), (0.1,5,5), maxiter=999999)
        fmax2 = fprob(xm, vrm, vtm, param)
        fmax  = max([fmax1, fmax2])


	#stepping over x, vr, vt
	xlim = rlim / rs
        x0   = 10**-7
        vmax = vesc(x0, [rlim, Vmax, rmax])			
	
	xarr  = np.append( np.linspace(x0, xlim, steps)  , xm)
	vrarr = np.append( np.linspace(vri, vmax, steps) , vrm)
	vtarr = np.append( np.linspace(vti, vmax, steps) , vtm)

	xarr  = xarr[np.argsort(xarr)]
	vrarr = vrarr[np.argsort(vrarr)]
	vtarr = vtarr[np.argsort(vtarr)]

	xxx, vrr, vtt = np.meshgrid(xarr, vrarr, vtarr, indexing='ij')

	#fillup the values.
	fvalues = np.zeros((steps+1, steps+1, steps+1))
	i = 0
	while i <= steps:
		j = 0
		while j <= steps:
			k = 0
			while k<= steps:
				fvalues[i,j,k] = fprob(xxx[i,j,k], vrr[i,j,k], vtt[i, j, k], param)
				k = k+1
			j=j+1
		i=i+1

	return xarr, vrarr, vtarr, fvalues

def G(param, steps):
	xarr, vrarr, vtarr, p = fbdvalues(param, steps)
	i = 1
	Gx = np.zeros(steps+1)
	Gxvr = []
	Gxvrvt = []
	pxvrvt = []
	while i <= steps:
		Gvr = np.zeros(steps+1)
		Gvrvt = []
		pvrvt = []
		j = 1
                while j <= steps:
			Gvt = np.zeros(steps+1)
			pvt = []
			k = 1
                        while k<= steps:
				gijk = max(p[i-1,j-1,k-1], p[i-1,j-1,k], p[i-1,j,k-1], p[i-1,j,k],
					   p[i  ,j-1,k-1], p[i  ,j-1,k], p[i  ,j,k-1], p[i  ,j,k] )
				gvt = gijk * (vtarr[k] - vtarr[k-1])
				Gvt[k] = Gvt[k-1] + gvt
				pvt.append(gijk)
                                k = k+1

			Gvrvt.append(Gvt/(max(Gvt)+10**-10))
			pvrvt.append(pvt)
			Gvr[j] = Gvr[j-1] + max(Gvt) * (vrarr[j] - vrarr[j-1])
                        j=j+1

		Gxvrvt.append(Gvrvt)
		pxvrvt.append(pvrvt)
		Gxvr.append(Gvr/(max(Gvr)+10**-10))
		Gx[i] = Gx[i-1] + max(Gvr)*(xarr[i]-xarr[i-1])
                i=i+1
	Gx = Gx/(max(Gx)+10**-10)
	return xarr, vrarr, vtarr, Gx, Gxvr, Gxvrvt, pxvrvt  #size steps+1 array, or list of that array,
						    	     #or list of such list


def sampleg(xarr, vrarr, vtarr, Gx, Gxvr, Gxvrvt, pxvrvt):
	ux  = rand.random()
	uvr = rand.random()
	uvt = rand.random()

	xupindex = sum(Gx<ux)
	xloindex = xupindex - 1
	xi = xarr[xloindex] + rand.random() * (xarr[xupindex] - xarr[xloindex])

	Gvr = Gxvr[xloindex]
	vrupindex = sum(Gvr < uvr)
	vrloindex = vrupindex - 1
	vri = vrarr[vrloindex] + rand.random() * (vrarr[vrupindex] - vrarr[vrloindex])

	Gvt = Gxvrvt[xloindex][vrloindex]
        vtupindex = sum(Gvt < uvt)
        vtloindex = vtupindex - 1
        vti = vtarr[vtloindex] + rand.random() * (vtarr[vtupindex] - vtarr[vtloindex])	

	gi = pxvrvt[xloindex][vrloindex][vtloindex]

	return xi, vri, vti, gi



#not efficient
def  g2(x1,x2, vr1, vr2, vt1, vt2, param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5]):
	f111 = fprob(x1,vr1, vt1, param)
	f112 = fprob(x1,vr1, vt2, param)
	f121 = fprob(x1,vr2, vt1, param)
	f122 = fprob(x1,vr2, vt2, param)
	f211 = fprob(x2,vr1, vt1, param)
	f212 = fprob(x2,vr1, vt2, param)
	f221 = fprob(x2,vr2, vt1, param)
	f222 = fprob(x2,vr2, vt2, param)
	return max(f111, f112, f121, f122, f211, f212, f221, f222)



#model function in terms of x, v, theta
def ftheta(x, v, theta, param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5]):

	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param

        rs = rmax/2.16                   # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0           # Vmax=0.465*sqrt(Ps)
        Pr = Ps * ( 1 - (np.log(1+x))/x )
        J  = x*rs * v * np.sin(theta)    #J = v * r*sin(theta)
        E  = v*v/2.0 + Pr    # v*v/2 + Pr


	xlim = rlim/rs
	Plim = Ps * ( 1 - (np.log(1+xlim))/xlim )
	Ec   = Ec * Ps
	Jb   = Jb * rs * (Ps**0.5) #*0.086


        if b <= 0:
                gJ = 1.0/(1 + (J/Jb)**-b)
        else:
                gJ = 1 + (J/Jb)**b

        N  = 1.0
        if E < Plim:
                hE = N*(E**a) * ((E**q + Ec**q)**(d/q)) * ((Plim - E)**e)
        else:
                hE = 0.0

        return gJ*hE


#define the general potential
def genphi(r, gamma = 1, beta = 3, alpha = 1, rmax = 1.5, Vmax = 21.0):
	rs = rmax/2.16
	Ps = (Vmax/0.465)**2.0
	
	def Integrand1(ri):
		xi = ri/rs
		I  = (xi**(-gamma)) * ( (1+xi**alpha)**((gamma-beta)/alpha) )
		#I  = 1.0/(xi * ((1+xi)**2))
		return -Ps * I * xi*xi

	def Integrand2(ri):
		xi = ri/rs
		I  = (xi**(-gamma)) * ( (1+xi**alpha)**((gamma-beta)/alpha) )
		#I  = 1.0/(xi * ((1+xi)**2))
		return -Ps * I * xi/rs
		

	phi1 = integrate.quad(Integrand1, 0.0, r)[0]
	phi2 = integrate.quad(Integrand2, r, np.inf)[0]
	x = r/rs
	return (phi1/r + phi2) 

'''
def Integrand1(ri):
		gamma = 1
		beta  = 3
		alpha = 1
		rs = 1.5/2.16
                xi = ri/rs
                I  = xi**(-gamma) * ( (1+xi**alpha)**((gamma-beta)/alpha) )
                return I * xi*xi

def Integrand2(ri):
                gamma = 1
                beta  = 3
                alpha = 1
		rs = 1.5/2.16
                xi = ri/rs
                I  = xi**(-gamma) * ( (1+xi**alpha)**((gamma-beta)/alpha) )
                return I * xi
'''
	
def phi(r, rmax = 1.5, Vmax = 21.0):
	rs = rmax/2.16
        Ps = (Vmax/0.465)**2.0
	x  = r/rs
	return Ps * (  - (np.log(1+x))/x )




#get escape velocity; param = [rlim, Vmax, rmax]
def vesc(x, param = [1.5, 21.0, 1.5]): #pop1 MR, pop2 MP

	rlim, Vmax, rmax = param

	rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)

	xlim = rlim/rs		 #turn rlim in unit of rs
	Pr   = Ps * ( 1 - (np.log(1+x))/x )
	Plim = Ps * ( 1 - (np.log(1+xlim))/xlim ) #lms 0.45*Ps
	
	vesc = (2 * (Plim - Pr))**0.5
	
	return vesc



#finding the max value of f(...) by evalute f spaning the whole parameter space
def findmax (param = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5]):
	flist = []
	dx    = 0.1
	dvr   = 2
	dvt   = dvr
	
	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
	rs    = rmax/2.16
	xmax  = rlim/rs
	x = 0.0001
	while x<xmax:  #loop through r/rs

        	vmax = vesc(x,[rlim, Vmax, rmax])
        	vr   = 0.00
        	while vr<vmax:  #loop through radial velocity

	                vtmax = (vmax * vmax - vr*vr)**0.5
	                vt    = 0.00
	                while vt < vtmax:
        	                fvt = fprob(x, vr, vt, param)
                	        flist.append(fvt)
                        	vt  = vt + dvt
       
       	        	vr = vr + dvr
		x = x + dx
	return max(flist)



def sample( param = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], samplesize = 3000):

	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
	rs    = rmax/2.16	# rmax=2.16*rs

	fmax1 = findmax(param) * 1.1
        xm,vrm, vtm = scipy.optimize.fmin(lambda (x,vr,vt): -fprob(x,vr,vt, param), (0.1,5,5), maxiter=999999)
        fmax2 = fprob(xm, vrm, vtm, param)
	fmax  = max([fmax1, fmax2])
        #print "fmax: ", fmax1, fmax2

	samplelist = []
	num = 1
	aux = 1
	while (num <= samplesize):

	        x = rlim * rand.random() / rs
		x0 = 10**-8
	        vmax0 = vesc(x0, [rlim, Vmax, rmax])
        	vr    = vmax0 * rand.random()
		vt    = vmax0 * rand.random()

                u  = rand.random() #new *fmax 
		fi = fprob(x, vr, vt, param)
		gi = fmax
        	if (fi/gi >= u): #new
                	samplelist.append( [x,vr,vt] )
                	num = num + 1

	        aux = aux + 1
	print "acceptance rate: ", num/(aux+0.00001)
	samplelist = np.asarray(samplelist)
	xarr  = samplelist[:,0]
	vrarr = samplelist[:,1]
	vtarr = samplelist[:,2]
	
	#tranform to the right coordinate
	r  = xarr * rs
	v  = (vrarr*vrarr + vtarr*vtarr )**0.5

	rsize = len(r)
	ru    = np.random.random_sample( rsize ) 
	theta = np.arccos(1-2*ru) #inverse sampling the distribution for theta
	phi   = np.random.random_sample( rsize ) * 2.0 * np.pi

	vsign = np.sign( np.random.random_sample( rsize ) - 0.5 )
	vphi  = np.random.random_sample( rsize ) * 2.0 * np.pi

	x = r * np.sin(theta) * np.cos(phi)
	y = r * np.sin(theta) * np.sin(phi)
	z = r * np.cos(theta)

	vz2 = vsign * vrarr
	vx2 = vtarr * np.cos(vphi)
	vy2 = vtarr * np.sin(vphi)

	#passive rotation, using rotation matrix 
	#to rotate the zhat of calculate velocity into the zhat of the spatial coordinate
	vx = np.cos(theta)*np.cos(phi)*vx2 - np.sin(phi)*vy2 + np.sin(theta)*np.cos(phi)*vz2
	vy = np.cos(theta)*np.sin(phi)*vx2 + np.cos(phi)*vy2 + np.sin(theta)*np.sin(phi)*vz2
	vz = -np.sin(theta)*vx2 + np.cos(theta)*vz2
	return x, y, z, vx, vy, vz

def sampleR( param = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], samplesize = 3000):

	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
	rs    = rmax/2.16	# rmax=2.16*rs

	fmax1 = findmax(param) * 1.1
        xm,vrm, vtm = scipy.optimize.fmin(lambda (x,vr,vt): -fprob(x,vr,vt, param), (0.1,5,5), maxiter=999999)
        fmax2 = fprob(xm, vrm, vtm, param)
	fmax  = max([fmax1, fmax2])
        #print "fmax: ", fmax1, fmax2

	samplelist = []
	num = 1
	aux = 1
	while (num <= samplesize):

	        x = rlim * rand.random() / rs
		x0 = 10**-8
	        vmax0 = vesc(x0, [rlim, Vmax, rmax])
        	vr    = vmax0 * rand.random()
		vt    = vmax0 * rand.random()

                u  = rand.random() #new *fmax 
		fi = fprob(x, vr, vt, param)
		gi = fmax
        	if (fi/gi >= u): #new
                	samplelist.append( [x,vr,vt] )
                	num = num + 1

	        aux = aux + 1
	print "acceptance rate: ", num/(aux+0.00001)
	samplelist = np.asarray(samplelist)
	xarr  = samplelist[:,0]
	vrarr = samplelist[:,1]
	vtarr = samplelist[:,2]
	
	#tranform to the right coordinate
	r  = xarr * rs
	v  = (vrarr*vrarr + vtarr*vtarr )**0.5

	rsize = len(r)
	ru    = np.random.random_sample( rsize ) 
	theta = np.arccos(1-2*ru) #inverse sampling the distribution for theta
	phi   = np.random.random_sample( rsize ) * 2.0 * np.pi

	vsign = np.sign( np.random.random_sample( rsize ) - 0.5 )
	vphi  = np.random.random_sample( rsize ) * 2.0 * np.pi

	x = r * np.sin(theta) * np.cos(phi)
	y = r * np.sin(theta) * np.sin(phi)
	z = r * np.cos(theta)

	vz2 = vsign * vrarr
	vx2 = vtarr * np.cos(vphi)
	vy2 = vtarr * np.sin(vphi)

	#passive rotation, using rotation matrix 
	#to rotate the zhat of calculate velocity into the zhat of the spatial coordinate
	vx = np.cos(theta)*np.cos(phi)*vx2 - np.sin(phi)*vy2 + np.sin(theta)*np.cos(phi)*vz2
	vy = np.cos(theta)*np.sin(phi)*vx2 + np.cos(phi)*vy2 + np.sin(theta)*np.sin(phi)*vz2
	vz = -np.sin(theta)*vx2 + np.cos(theta)*vz2
	#return x, y, z, vx, vy, vz
	return [list(x),list(y),list(z),list(vx),list(vy),list(vz)]


def sample2( param = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], steps = 10, samplesize = 3000):

        a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs    = rmax/2.16       # rmax=2.16*rs

        #fmax1 = findmax(param) * 1.1
        #xm,vrm, vtm = scipy.optimize.fmin(lambda (x,vr,vt): -fprob(x,vr,vt, param), (0.1,5,5), maxiter=999999)
        #fmax2 = fprob(xm, vrm, vtm, param)
        #fmax  = max([fmax1, fmax2])
        #print "fmax: ", fmax1, fmax2

        #new    
        #steps = 10
        xarr, vrarr, vtarr, Gx, Gxvr, Gxvrvt, pxvrvt = G(param, steps)

        samplelist = []
        num = 1
        aux = 1
        while (num <= samplesize):

                x = rlim * rand.random() / rs

                x0 = 10**-7
                #vmax0 = vesc(x0, [rlim, Vmax, rmax])
                #vr    = vmax0 * rand.random()
                #vt    = vmax0 * rand.random()
                #vtmax = (vmax * vmax - vr*vr)**0.5
                #vt    = vtmax * rand.random()

                u  = rand.random() #new *fmax
		xi, vri, vti, gi = sampleg(xarr, vrarr, vtarr, Gx, Gxvr, Gxvrvt, pxvrvt) 
                fi = fprob(xi, vri, vti, param)
                if (fi/gi >= u): #new
                        samplelist.append( [xi,vri,vti] )
                        num = num + 1
		if (fi/gi > 1.0):
			print 'Warning: f(x)/g(x) > 1'

                aux = aux + 1
        print "acceptance rate: ", num/(aux+0.00001)
        samplelist = np.asarray(samplelist)
        xarr  = samplelist[:,0]
        vrarr = samplelist[:,1]
        vtarr = samplelist[:,2]

        #tranform to the right coordinate
        r  = xarr * rs
        v  = (vrarr*vrarr + vtarr*vtarr )**0.5

        rsize = len(r)
	ru    = np.random.random_sample( rsize )
        theta = np.arccos(1-2*ru) #inverse sampling the distribution for theta
        phi   = np.random.random_sample( rsize ) * 2.0 * np.pi

        vsign = np.sign( np.random.random_sample( rsize ) - 0.5 )
        vphi  = np.random.random_sample( rsize ) * 2.0 * np.pi

        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        vz2 = vsign * vrarr
        vx2 = vtarr * np.cos(vphi)
        vy2 = vtarr * np.sin(vphi)

	#passive rotation, using rotation matrix 
        #to rotate the zhat of calculate velocity into the zhat of the spatial coordinate
        vx = np.cos(theta)*np.cos(phi)*vx2 - np.sin(phi)*vy2 + np.sin(theta)*np.cos(phi)*vz2
        vy = np.cos(theta)*np.sin(phi)*vx2 + np.cos(phi)*vy2 + np.sin(theta)*np.sin(phi)*vz2
        vz = -np.sin(theta)*vx2 + np.cos(theta)*vz2
        return x, y, z, vx, vy, vz



#sampling realistic mock data with metalicity and measurement error
def mockdata(param1, param2, N = 3000, f1 = 0.5, gamma = 1.0): 
	x1, y1, z1, vx1, vy1, vz1 = sample( param1, N*f1 )
	x2, y2, z2, vx2, vy2, vz2 = sample( param2, N*(1-f1) )
	
	#adding mean metalicity and observation error to the mock matallicity
	mstdev = (0.5**2 + 0.2**2)**0.5   #total stdev from true metalicity and observation error
	m1 = np.random.normal(1, mstdev, np.size(x1)) #mean metallicity check with matt
	m2 = np.random.normal(-2, mstdev, np.size(x2))
	
	m = np.hstack((m1,m2))
	x = np.hstack((x1,x2))
	y = np.hstack((y1,y2))
	z = np.hstack((z1,z2))
	vx = np.hstack((vx1,vx2))
	vy = np.hstack((vy1,vy2))
	vz = np.hstack((vz1,vz2))
	vz = vz + np.random.normal(0,2.0, np.size(vz))
	
	return x, y, z, vx, vy, vz, m

#mock data without error in velocity and metallicity
def mockdata_noerror(param1, param2, N = 3000, f1 = 0.5, gamma = 1.0):
        x1, y1, z1, vx1, vy1, vz1 = sample( param1, N*f1 )
        x2, y2, z2, vx2, vy2, vz2 = sample( param2, N*(1-f1) )

        #adding mean metalicity and observation error to the mock matallicity
        mstdev = 0.5   #total stdev from true metalicity and observation error
        m1 = np.random.normal(1, mstdev, np.size(x1)) #mean metallicity check with matt
        m2 = np.random.normal(-2, mstdev, np.size(x2))

        m = np.hstack((m1,m2))
        x = np.hstack((x1,x2))
        y = np.hstack((y1,y2))
        z = np.hstack((z1,z2))
        vx = np.hstack((vx1,vx2))
        vy = np.hstack((vy1,vy2))
        vz = np.hstack((vz1,vz2))
        #vz = vz + np.random.normal(0,2.0, np.size(vz))

        return x, y, z, vx, vy, vz, m




#integrating the orbit. time in Gyr 
def orbit(x,y,z, vx,vy,vz, orbittime = 2, dt = 0.001, rmax=1.5, Vmax=21.0):

	rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)

	convertF  = 1.02269032     #constant to convert force into km/s / Gyr
	convertv  = 1.02269032    	   #constand to convert km/s into kpc/Gyr


	t = 0
	while t <= orbittime:
		r = (x*x + y*y + z*z)**0.5
		drdx = x/r
        	drdy = y/r
        	drdz = z/r

		#calculate force or acceleration, mass = 1
		xi = r/rs
		#dphidr = Ps * ( -(1+xi)**-2 - 1.0/(xi+xi*xi) + np.log(1+xi)/(xi*xi) )/rs 

        	dphidr = -Ps * ( ( 1.0/(1.0*r + r*r/rs) )  -  ( rs*np.log(1.0+r/rs) / (r*r) ) )
        	Fx = -dphidr * drdx * convertF
        	Fy = -dphidr * drdy * convertF
        	Fz = -dphidr * drdz * convertF

      		x  = x + vx*convertv*dt/2.0
		y  = y + vy*convertv*dt/2.0
		z  = z + vz*convertv*dt/2.0
		vx = vx + Fx*convertF*dt
		vy = vy + Fy*convertF*dt
		vz = vz + Fz*convertF*dt
		x  = x + vx*convertv*dt/2.0
                y  = y + vy*convertv*dt/2.0
                z  = z + vz*convertv*dt/2.0
		t  = t + dt
		#print 'time: ', t, orbittime 
	return x,y,z, vx,vy,vz


#calculate projected density, using f(x,vr,vt)
def projecteddensity(param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5]):

	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)

        xlim = rlim/rs  #in unit of rs
        dx = 0.001 #*.7kpc * 10
        dX = 0.1
        projrholist = []
        Xlist = []
        Xi = 0.00
        while Xi<xlim:
            Integrantlist = []
            IntegrantNormlist = []
            xlist = []
            xi = Xi+0.00001
            while xi<xlim:  #loop through r/rs

		def rhofunc(vr,vt):
                        return f(xi, vr, vt, param) * vt


		vmax = vesc(xi,[rlim, Vmax, rmax])
                def bounds_vt():
                        return [0, vmax]

                def bounds_vr(vt):
                        lim = (vmax*vmax - vt*vt)**0.5
                        return [0, lim]


                #marginalize over theta and velocity, as a function of x
                rho     = integrate.nquad(rhofunc, [bounds_vr, bounds_vt])[0]

                #integrants to be integrated over x, as a function of X (or R/rs), 
                zz = (xi*xi - Xi*Xi)
                dz = rs*xi / np.sqrt(zz)
                Integrant = rho*dz

                Integrantlist.append(Integrant)
                xlist.append(xi)

                xi = xi + dx

            #integrate over X as a function of X
            #print 'time: ', (time.time() - start_time)/60.0
            projrho = integrate.trapz(Integrantlist,xlist)
            projrholist.append( projrho )
            Xlist.append(Xi)
            print Xi, projrho
            Xi = Xi + dX
        return np.asarray(Xlist) * rs,  np.asarray(projrholist)/max(projrholist) 
	#np.sqrt( losvsig/(norm+0.00000001) )

#calculate projected density,  by defining rho(r)
def projecteddensity2(param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], dR = 0.02):

        a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)

        xlim = rlim/rs  #in unit of rs

	dX = dR/rs
        projrholist = []
        Xlist = []
        Xi = 0.00
        while Xi<xlim:
		def Integrand(r):   #It's using r, because I define below has limits in r, not x
			xi = r/rs
			def rhofunc(vt,vr):
                        	return f(xi, vr, vt, param) * vt

                	vmax = vesc(xi,[rlim, Vmax, rmax])
                	def bounds_vr():
                        	return [0, vmax]

                	def bounds_vt(vr):
                        	lim = (vmax*vmax - vr*vr)**0.5
                        	return [0, lim]

			rho = integrate.nquad(rhofunc, [bounds_vt, bounds_vr])[0]
			zz = (xi*xi - Xi*Xi)
                	dz = rs*xi / np.sqrt(zz)
			return rho*dz


		I  = integrate.quad(Integrand, Xi*rs, rlim)[0]
		projrholist.append(I)
            	Xlist.append(Xi)
		print 'fvrvt Xi: ', Xi, I
		Xi = Xi + dX
	return np.asarray(Xlist) * rs,  np.asarray(projrholist)/max(projrholist)



#calculate projected density, by ftheta(), defining rho(r)
def projden(param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], dR = 0.02):

        a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)

        xlim = rlim/rs  #in unit of rs

        dX = dR/rs
	projrholist = []
        Xlist = []
        Xi = 0.00
	while Xi<xlim:
                def Integrand(r):
                        xi = r/rs 
			def rhof(theta,v):
	                        return ftheta(xi, v, theta, param) * v*v*np.sin(theta)
	
			vmax = vesc(xi,[rlim, Vmax, rmax])
			bounds_theta = [0, np.pi]
                	bounds_v     = [0, vmax]

			#marginalize over theta and velocity, as a function of x
	                rho = integrate.nquad(rhof, [bounds_theta, bounds_v])[0]
			zz  = (xi*xi - Xi*Xi)
                        dz  = rs*xi / np.sqrt(zz)
                        return rho*dz

		I = integrate.quad(Integrand, Xi*rs, rlim)[0]
		projrholist.append(I)
                Xlist.append(Xi)
                print "ftheta Xi: ", Xi, I
                Xi = Xi + dX
        return np.asarray(Xlist) * rs,  np.asarray(projrholist)/max(projrholist)


#it returns functions: rho(r), rho*vrdispersion^2(r), rho*vtdispersion^2(r), using ftheta()
def rhosigrt2theta(r, param = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5]):	
	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs

	xi = r/rs
        def sigrI(theta, v):
		return ftheta(xi, v, theta, param) * pow(v,4) * pow( np.cos(theta),2 ) * np.sin(theta)

        def sigtI(theta, v):
                return ftheta(xi, v, theta, param) * pow(v,4) * pow( np.sin(theta) , 3) 

        vmax = vesc(xi,[rlim, Vmax, rmax])
	bounds_theta = [0, np.pi]
        bounds_v     = [0, vmax]

        rhosigr2 = integrate.nquad(sigrI, [bounds_theta, bounds_v])[0]
        rhosigt2 = integrate.nquad(sigtI, [bounds_theta, bounds_v])[0]
        return rhosigr2, rhosigt2 

#it returns function rho(r) using ftheta()
def rhortheta(r, param = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5]):
        a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs

        xi = r/rs
        def rhorI(theta, v):
                return ftheta(xi, v, theta, param) * pow(v,2) * np.sin(theta)

        vmax = vesc(xi,[rlim, Vmax, rmax])
        bounds_theta = [0, np.pi]
        bounds_v     = [0, vmax]

        rhor = integrate.nquad(rhorI, [bounds_theta, bounds_v])[0]
        return rhor

#it calculate rhoR using ftheta()
def SigR(param):
	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs
	xlim = rlim/rs
	
	def rhoRI(theta, v, x, X):
		z   = (x*x - X*X +10**-10)**.5
		aux = 2 * ftheta(x, v, theta, param) * pow(v,2) * np.sin(theta) * 2 * np.pi  * rs
		return aux * x/z

	def bounds_theta(v,x):
                return [0, np.pi]
        def bounds_v(x):
                vmax = vesc( x, [rlim, Vmax, rmax] )                        
		return [0, vmax]
	#def bounds_r(X):
        #        return [X, xlim]

	#def result(Xi):
	#	ans = integrate.nquad(rhoRI, [bounds_theta, bounds_v, bounds_r], args=(Xi,))[0]
	#	return ans

	Xarr = np.linspace(0, xlim, 10)
	projlist = []
	for Xi in Xarr:
		def rhoRI2(theta,v,x):
			return rhoRI(theta,v,x,Xi)
		ans = integrate.nquad(rhoRI2, [bounds_theta, bounds_v, [Xi,xlim]])[0]
		projlist.append(ans)
		print 'SigR Xi: ', Xi
	proj = np.array(projlist)

	return Xarr*rs, proj/max(proj)


#calculate the projdensity, los vdispersion as an function of R. -- NOT  working as of now, use other one
def callosvsigprojdentheta(param = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], dR = 0.03):
	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)
        losvsiglist = []
        projdenlist = []
        Rlist = []
        Ri = 0.0000
        while Ri<rlim:	#for each Ri, integrate over r

            	print 'ftheta Ri: ', Ri, rlim
		def losvdispI(r):
			rhosigr2, rhosigt2 = rhosigrt2theta(r, param)
			zz = r*r - Ri*Ri + 1*10**-10
			dz = r / np.sqrt(zz)
			return (zz*rhosigr2 + 0.5*Ri*Ri*rhosigt2) * dz / (r*r)

		def projdenI(r):
			rhor = rhortheta(r, param)
			zz   = r*r - Ri*Ri + 1*10**-10
                        dz   = r / np.sqrt(zz)
			return rhor*dz

		dr = .001 #lms
		ri = Ri + 10**-10
		rlist        = []
		losvIlist    = []
		projdenIlist = []
		while ri <= rlim:
			losvIlist.append(   losvdispI(ri) )
			projdenIlist.append( projdenI(ri) )
			rlist.append(ri)
			ri = ri + dr


		losvsiglist.append( integrate.trapz(losvIlist, rlist) )
		projdenlist.append( integrate.trapz(projdenIlist, rlist) )

		#integration limit for r is from Ri to rlim, 
            	#losvsiglist.append( integrate.romberg(losvdispI, Ri, rlim) )
            	#projdenlist.append( integrate.romberg(projdenI,  Ri, rlim) )
            	Rlist.append(Ri)
            	Ri = Ri + dR

        Rarr = np.asarray(Rlist)
        projdenarr  = np.asarray(projdenlist)
        losvsigarr2 = np.asarray(losvsiglist) / (projdenarr + 10**-10)
        return Rarr, losvsigarr2**0.5, projdenarr/(max(projdenarr) + 10**-10)




#get projected surface density with given data points
def getprojdensity (x, y, opt = 1):
	R = (x*x+y*y)**0.5
	R = R[np.argsort(R)]

        Rlist     = []
        countlist = []
        Ri = R[0]

	if (opt == 1):
	    i  = 0
            binsize = len(x)/100
            while i <= len(R):
                if (i+binsize <= len(R)-1):
                        num  = len(R[i:(i+binsize)])
			area = np.pi * ( (R[i+binsize])**2 - (R[i])**2 )
                        countlist.append(num/area)
                        Rmean = np.mean(R[i:(i+binsize)])
                        Rlist.append(Rmean) 
                elif (i < len(R)):
                        num = len(R[i:len(R)])
			area = np.pi * ( (max(R))**2 - (R[i])**2 )
                        countlist.append(num/area)
                        Rmean = np.mean(R[i:len(R)])
                        Rlist.append(Rmean) 

                i = i + binsize

            return np.asarray(Rlist), np.asarray(countlist)


	elif (opt == 2):
	    dR = .1
            while Ri < max(R):
                if(Ri+dR <= max(R)):
                        R1 = R[(R <= Ri+dR)]
                        if(len(R1) > 0):
                                Rmean = np.mean(R1)  #Ri+dR*.5 #np.mean( R1 )
                                Rarea = ( .5*(Ri**2 + (Ri+dR)**2) )**(1./2)
                                Rlist.append(Rarea)
                                area = (Ri+dR)**2 - (Ri)**2
                                countlist.append(len(R1)/(area*np.pi))
                                R = R[len(R1):len(R)]
                elif(Ri < max(R) and Ri+dR > max(R) ):
                        R1 = R[(R>Ri)]
                        Rmean = np.mean(R1) #Ri+dR*.5 #np.mean( R1 )
                        Rarea = ( .5*(Ri**2 + (Ri+dR)**2) )**(1./2)
                        Rlist.append(Rarea)
                        area = (Ri+dR)**2 - (Ri)**2
                        countlist.append(len(R1)/(area*np.pi))
                Ri = Ri+dR

            return np.array(Rlist), np.array(countlist)

	else:
		print "Error: choose bin type, either 1 or 2"
		return 0


#get los velocity dispersion with given data
def getlosvdisp(x,y,vz):

	R  = (x*x+y*y)**0.5
	Rargsort = np.argsort(R)
	R  = R[Rargsort]
	VZ = vz[Rargsort]

	Rlist     = []
	VZstdlist = []
	Ri = 0.0
	#dR = 0.1
	i  = 0
	binsize = len(VZ)/30
	while i <= len(R):

	        if (i+binsize <= len(VZ)-1):
        	        std = np.std(VZ[i:(i+binsize)], ddof=0)
                	VZstdlist.append(std)
                	Rmean = np.mean(R[i:(i+binsize)])
                	Rlist.append(Rmean) #R[i+binsize/2])
        	elif (i < len(VZ)):
                	std = np.std(VZ[i:len(VZ)], ddof=0)
                	VZstdlist.append(std)
                	Rmean = np.mean(R[i:len(VZ)])
                	Rlist.append(Rmean) # R[(i + (len(VZ)-i)/2)])

        	i   = i + binsize

	return np.asarray(Rlist), np.asarray(VZstdlist)


#calculate los v dispersion: 
def callosvdisp(param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], dR = 0.02):

	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)

        xlim = rlim/rs  #in unit of rs

        dX = dR/rs
        projrholist = []
	losvsiglist = []
        Xlist = []
        Xi = 0.00

        while Xi<xlim:
                def IntegrandNorm(r):   #It's using r, because I define below has limits in r, not x
                        xi = r/rs
                        def rhofunc(vt,vr):
                                return f(xi, vr, vt, param) * vt * 2*np.pi

                        vmax = vesc(xi,[rlim, Vmax, rmax])
                        def bounds_vr():
                                return [0, vmax]

                        def bounds_vt(vr):
                                lim = (vmax*vmax - vr*vr)**0.5
                                return [0, lim]

                        rho = integrate.nquad(rhofunc, [bounds_vt, bounds_vr])[0]
                        zz = (xi*xi - Xi*Xi)
                        dz = rs*xi / np.sqrt(zz)
                        return 2*rho*dz


		def Integrandlosv(r):   #It's using r, because I define below has limits in r, not x
                        xi = r/rs
                        def sigr(vt,vr):
                                return f(xi, vr, vt, param) * vt * (vr*vr) * 2*np.pi 

			def sigt(vt,vr):
                                return f(xi, vr, vt, param) * vt * (vt*vt) * 2*np.pi
			
			def rhofunc2(vt,vr):
                                return f(xi, vr, vt, param) * vt * 2*np.pi		


                        vmax = vesc(xi,[rlim, Vmax, rmax])
                        def bounds_vr():
                                return [0, vmax]  #should be +/-vmax, but the function is symmetric, so.

                        def bounds_vt(vr):
                                lim = (vmax*vmax - vr*vr)**0.5
                                return [0, lim]

                        rhosigr2 = integrate.nquad(sigr, [bounds_vt, bounds_vr])[0]
			rhosigt2 = integrate.nquad(sigt, [bounds_vt, bounds_vr])[0]
			rhor     = integrate.nquad(rhofunc2, [bounds_vt, bounds_vr])[0]
                        zz = (xi*xi - Xi*Xi)
                        dz = rs*xi / np.sqrt(zz)
			result = ( zz*rhosigr2 + 0.5*Xi*Xi*rhosigt2 ) * dz / (xi*xi)  #!!! 0.5 or no???
                        return 2*result


                norm  = integrate.quad(IntegrandNorm, Xi*rs, rlim)[0]
		I     = integrate.quad(Integrandlosv, Xi*rs, rlim)[0]
                losvsiglist.append(I/norm)
		projrholist.append(norm)
                Xlist.append(Xi)
                print 'fvrvt Xi: ', Xi, I
                Xi = Xi + dX
        return np.asarray(Xlist)*rs,(np.asarray(losvsiglist))**0.5, np.asarray(projrholist)/(max(projrholist)+10**-8)



#it returns R, losvdispersion, and projected density
def callosvdisp2(param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], dR=.02, dx=.05, dv=0.1):

        a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)

        xlim = rlim/rs  #in unit of rs

        dX = dR/rs
	losvsiglist = []
	projdenlist = []
	Xlist = []
        Xi = 0.00
	
	xi, rhor, rhosigr2, rhosigt2 = rhosigrt2(param, dx, dv)
        while Xi<xlim:  #for each Xi, do x, vr, vt, integration
		
	    #integration limit is from R to inf, 
	    ulim = len(xi)
            llim = sum(xi<Xi)
	    if ulim>llim:
		xi 	 = xi[llim:ulim]
		rhor 	 = rhor[llim:ulim]
		rhosigr2 = rhosigr2[llim:ulim]
		rhosigt2 = rhosigt2[llim:ulim]
		zz   = xi*xi - Xi*Xi
		dz   = rs*xi / np.sqrt(zz)

		#for each R, integrate over r
		losvIntegrant = rhor * ( zz*rhosigr2 + Xi*Xi*rhosigt2 ) * dz / (xi*xi)  
		projdenIntegrant = rhor*dz
		losvsiglist.append( integrate.simps(losvIntegrant, xi) )
		projdenlist.append( integrate.simps(projdenIntegrant, xi) )		

		Xlist.append(Xi)
		print 'Xi: ', Xi
	    Xi = Xi + dX

	Rarr = rs*np.asarray(Xlist)
	projdenarr = np.asarray(projdenlist)
	losvsigarr2 = np.asarray(losvsiglist) / (projdenarr+10**-8)
	return Rarr, losvsigarr2**0.5, projdenarr/(max(projdenarr) + 10**-8)


#it returns 3 arries: x, rho*vrdispersion^2, rho*vtdispersion^2
def rhosigrt2(param = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], dx = 0.05, dv = 0.1):
	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
        rs = rmax/2.16           # rmax=2.16*rs
        Ps = (Vmax/0.465)**2.0   # Vmax=0.465*sqrt(Ps)

        xlim = rlim/rs  #in unit of rs

        dvr= dv
        dvt= dvr
	rhorlist = []
        sigrlist = []
	sigtlist = []
        xlist    = []
        xi       = 0.0000001
        while xi<xlim: #for each xi, do vr, vt, integration

            vmax = vesc(xi,[rlim, Vmax, rmax])
	    sigrlistr = []
            rhorlistr = []
            sigtlistr = []
            vrlist    = []
            vr        = 0.0
            while vr < vmax: #for each vr, do vt integration

	            vtmax = (vmax * vmax - vr*vr)**0.5
                    sigrlistt = []
		    rhorlistt = []
		    sigtlistt = []
                    vtlist    = []
                    vt        = 0.0
                    while vt < vtmax:
			f0   = f(xi, vr, vt, param) 
          	        rhor = f0 * vt
			sigr = rhor * vr * vr
			sigt = rhor * vt * vt
			vt = vt + dvt  #vt increment
			rhorlistt.append(rhor)
			sigtlistt.append(sigt)
			sigrlistt.append(sigr)
			vtlist.append(vt)

		    #for each vr, integrate over vt	
		    rhorlistr.append( integrate.simps(rhorlistt,vtlist) )
		    sigrlistr.append( integrate.simps(sigrlistt,vtlist) )
		    sigtlistr.append( integrate.simps(sigtlistt,vtlist) )
		    vrlist.append(vr)
		    vr = vr + dvr #vr increment

	    #for each xi, integrate over vr
	    rhorlist.append( integrate.simps(rhorlistr,vrlist) )
            sigrlist.append( integrate.simps(sigrlistr,vrlist) )
            sigtlist.append( integrate.simps(sigtlistr,vrlist) )
            xlist.append(xi)
	    xi = xi + dx

	xarr = np.asarray(xlist)
	rhorarr = np.asarray(rhorlist)
	sigrarr = np.asarray(sigrlist)
	sigtarr = np.asarray(sigtlist)
	return xarr, rhorarr, sigrarr, sigtarr #note sigr/sigt is rho*sig^2


#integrate directly to get velocity dispersion (WRONG: can't fuse rho(r) integrand with sigrt integrand) 
def losvprojden(param = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5], dR = 0.05):
	a,d,e, Ec, rlim, b, q, Jb, Vmax, rmax = param
	rs = rmax/2.16           # rmax=2.16*rs

	losvsiglist = []
        projdenlist = []
        Rlist = []
        Ri = 0.0000
        while Ri<rlim:  #for each Ri, integrate over r
		def losvI(theta, v,r):
			zz = r*r - Ri*Ri + 10**-10  #add 10**-10 to prevent divide by zero
                	dz = r / np.sqrt(zz)
			f  = ftheta( (r/rs), v, theta, param )
			Ja = f * zz*pow(v,4)*np.cos(theta)*np.cos(theta)*np.sin(theta)/(r*r)
			Jb = f * Ri*Ri*pow(v,4)*pow(np.sin(theta),3)/(r*r)
			print 'using Ri', Ri, theta, v, r
			#missing rhor again! WRONG
			rhor = ftheta( (r/rs), v, theta, param ) * v*v * np.sin(theta) 
			return (Ja+Jb)*dz

		def projI(theta,v,r):
			zz = r*r - Ri*Ri + 10**-10
                        dz = r / np.sqrt(zz)
			I  = ftheta( (r/rs), v, theta, param ) * v*v * np.sin(theta) * dz
			return I
	
		def bounds_theta(v,r):
			return [0, np.pi]
		def bounds_v(r):
			vmax = vesc( (r/rs), [rlim, Vmax, rmax] )
	 		return [0, vmax]
		def bounds_r():
			return [Ri, rlim]
		
		def vmax(r):
                        return vesc( (r/rs), [rlim, Vmax, rmax] )

		def vmin(r):
                        return 0.0
		def tmax(v,r):
			return np.pi
		def tmin(v,r):
			return 0.0

		#Ilos  = integrate.nquad( losvI, [bounds_theta, bounds_v, bounds_r], epsrel = 1e-6 )[0] 
		#Iproj = integrate.nquad( projI, [bounds_theta, bounds_v, bounds_r], epsrel = 1e-6 )[0]
		
		Ilos = integrate.tplquad(losvI, Ri, rlim, vmin,vmax, tmin,tmax, epsrel=1e-4, epsabs=0)[0]
		Iproj= integrate.tplquad(projI, Ri, rlim, vmin,vmax, tmin,tmax, epsrel=1e-4, epsabs=0)[0]

		losvsiglist.append( Ilos  )
		projdenlist.append( Iproj )
		Rlist.append(Ri)
		print "3quad ftheta: ", Ri, rlim
		Ri = Ri + dR

	Rarr	    = np.asarray(Rlist)
	projdenarr  = np.asarray(projdenlist)
        losvsigarr2 = np.asarray(losvsiglist) / (projdenarr+10**-8)
        return Rarr, losvsigarr2**0.5, projdenarr/(max(projdenarr) + 10**-8)
	



def romberg( f, a, b, n = 100 ):
    """Estimate the integral of f(x) from a to b using Romberg Integration.

    USAGE:
        r = romberg( f, a, b, n )

    INPUT:
        f       - function to integrate,
        [a, b]  - the interval of integration,
        n       - number of levels of recursion

    OUTPUT:
        numpy float array - Romberg integration array; most accurate
                            answer should be at bottom of right-most column,

    NOTES:
        Based on an algorithm in "Numerical Mathematics and Computing"
        4th Edition, by Cheney and Kincaid, Brooks-Cole, 1999.

    AUTHOR:
        Jonathan R. Senning <jonathan.senning@gordon.edu>
        Gordon College
        February 17, 1999
        Converted to Python August 2008
    """

    r = np.array( [[0] * (n+1)] * (n+1), float )
    h = b - a
    r[0,0] = 0.5 * h * ( f( a ) + f( b ) )

    powerOf2 = 1
    for i in xrange( 1, n + 1 ):

        # Compute the halved stepsize and use this to sum the function at
        # all the new points (in between the points already computed)

        h = 0.5 * h

        sum = 0.0
        powerOf2 = 2 * powerOf2
        for k in xrange( 1, powerOf2, 2 ):
            sum = sum + f( a + k * h )

        # Compute the composite trapezoid rule for the next level of
        # subdivision.  Use Richardson extrapolation to refine these values
        # into a more accurate form.

        r[i,0] = 0.5 * r[i-1,0] + sum * h

        powerOf4 = 1
        for j in xrange( 1, i + 1 ):
            powerOf4 = 4 * powerOf4
            r[i,j] = r[i,j-1] + ( r[i,j-1] - r[i-1,j-1] ) / ( powerOf4 - 1 )
    print size(r)
    return r
