from matplotlib import rc
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.special import gamma
import os

filepath = os.path.dirname(os.path.abspath(__file__))

f=np.genfromtxt(filepath + "\\CovidData.txt", int, "#", "," , skip_header=0)

newcases=f[:,4]
cumucases=f[:,5]
newdeaths=f[:,6]
cumudeaths=f[:,7]


print(len(newcases))
g=0
while g <len(newcases):
    if newcases[g]<1:
        newcases=np.delete(newcases,g)
    g+=1
dd=np.linspace(0,len(newcases),len(newcases))

plt.step(dd,newcases)

plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside

#xmin,xmax,ymin,ymax=plt.axis([0,475,0,32000])

plt.xlabel('Days') # Labels x-axis

plt.ylabel('Number of New Cases') # Labels y-axis
plt.show()

mus=[]

N=len(newcases)
days=np.linspace(1,N,N)
maxs=[14000,21000,26485,37875]
ends=[0,217,218,418,419,634,635,765]
sigmas=[]
h=0
c=0
while h<3:
    udata=[]
    ydata=[]
    xdata=[]
    tdata=[]
    
    g=0
    while g<N:
        if days[g]>ends[c] and days[g]<ends[c+1]:
            udata.append(np.sqrt(newcases[g]))
            ydata.append(newcases[g])
            tdata.append(days[g])
        if newcases[g]==ends[c+1]:
                break
        g+=1
        
    # def f(x, A, mu, sigma):  
    #     return A*np.exp(-((x-mu)**2)/(2*sigma**2))
    
    # A0 = maxs[h]
    # mu0 = (ends[c+1]-ends[c])/2
    # sigma0 = (ends[c+1]-ends[c])/(2*np.sqrt(2))
    
    
    if h==3:
        def f(x, A, mu, sigma,D):
            return A*np.exp(-((x-mu)**2)/(2*sigma**2))+D
        A0 = maxs[h]
        mu0 = (ends[c+1]-ends[c])/2
        sigma0 = (ends[c+1]-ends[c])/(2*np.sqrt(2))
        D0=100
        p0 = [A0, mu0, sigma0,D0]
        name = ['A', 'mu', 'sigma','Vertical Shift']
    else:
        def f(x, A, mu, sigma):  
            return A*np.exp(-((x-mu)**2)/(2*sigma**2))
        A0 = maxs[h]
        mu0 = (ends[c+1]-ends[c])/2
        sigma0 = (ends[c+1]-ends[c])/(2*np.sqrt(2))
        p0 = [A0, mu0, sigma0]
        name = ['A', 'mu', 'sigma']
    
    tmodel = np.linspace(ends[c], ends[c+1], 1000)
    
    ystart=f(tmodel,*p0)
    
    
    ## Curve Fit Function---------------------------------------------------------------------------------------------------------------------------------
    
    popt,pcov=curve_fit(f,tdata,ydata,p0,sigma=udata,absolute_sigma=True)
    dymin = (ydata-f(tdata,*popt))/udata
    min_chisq = sum(dymin*dymin)
    dof=len(tdata) - len(popt)
    
    print('Chi square: ',min_chisq)
    print('Number of degrees of freedom: ',dof)
    print('Chi square per degree of freedom: ',min_chisq/dof)
    print()
    
    mn=min_chisq / dof
    
    print('Fitted parameters with 68% C.I.: ')
    
    for i, pmin in enumerate(popt):
        print('%2i %-10s %12f +/- %10f'%(i, name[i], pmin, np.sqrt(pcov[i,i])*np.sqrt(mn)))
        
    print()
    
    perr=np.sqrt(np.diag(pcov))
    print(perr)
    
    print("Correlation matrix")
    print("               ")
    
    for i in range(len(popt)): print('%-10s'%(name[i],)),
    print()
    
    for i in range(len(popt)):
        print('%10s'%(name[i])),
        for j in range(i+1):
            print('%10f'%(pcov[i,j]/np.sqrt(pcov[i,i]*pcov[j,j]),)),
        print()
    
    mus.append(popt[1])
    sigmas.append(popt[2])
    yfit=f(tmodel,*popt)
    
    
    ## Plotting-------------------------------------------------------------------------------------------------------------------------------------------
    plt.scatter(tdata,ydata)

    plt.xlabel('Days') # Labels x-axis

    plt.ylabel('New Cases') # Labels y-axis    
    
    plt.plot(tmodel,yfit,'-r') # Plots 
    
    
    plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside
    
    plt.show()    
    
    
    h+=1
    c+=2
    
p0=[]

print(mus)
print(sigmas)

mus.append(days[-1])
mus.insert(0,0)
sigmas.insert(0,0)

sigmas.append(0)
u_lambdas=[]
lambdas=[]
endpoints1=[]
endpoints2=[]
d=0
print(mus,sigmas)
#while d<5:
    
    #if d==0:   #d==0:#     
        #endpoints1.append(mus[0])
    
    #elif d%2==0:
        #endpoints1.append(mus[d-1]+3*sigmas[d-1])
    
    ##elif d%2==0:
        ##endpoints1.append(mus[d-1]+3*sigmas[d-1])
    #else:
        #endpoints2.append(mus[d])
    #d+=1
endpoints1.append(mus[0])
endpoints1.append(mus[1]+3*sigmas[1])
endpoints1.append(mus[2]+3*sigmas[2])

endpoints2.append(mus[1])
endpoints2.append(mus[2])
endpoints2.append(mus[3])



print(endpoints1,endpoints2)
for h in range(3):

    i=0
  
    udata=[]
    ydata=[]
    xdata=[]
    tdata=[]
    
    while i<N:
        if days[i]>endpoints1[h] and days[i]<endpoints2[h]:
            udata.append(np.sqrt(newcases[i]))
            ydata.append(newcases[i])
            tdata.append(days[i])
        if days[i]==endpoints2[h]:
            break
        i=i+1    
    
    
    
    def fx(x, A, lambdaa, D):  
        return A*np.exp(np.asarray(lambdaa)*x)+D
    
    A0 = 0
    lambda0 = 0
    D0=0
    p0 = [A0, lambda0, D0]
    name = ['A', 'lambda','D']
    
    
    
    tmodel = np.linspace(endpoints1[h], endpoints2[h], 1000)
    
    
    ystart=fx(tmodel,*p0)
    
    
    ## Curve Fit Function---------------------------------------------------------------------------------------------------------------------------------
    
    popt,pcov=curve_fit(fx,tdata,ydata,p0,sigma=udata,absolute_sigma=True)
    dymin = (ydata-fx(tdata,*popt))/udata
    min_chisq = sum(dymin*dymin)
    dof=len(tdata) - len(popt)
    
    print(len(popt))
    print('Chi square: ',min_chisq)
    print('Number of degrees of freedom: ',dof)
    print('Chi square per degree of freedom: ',min_chisq/dof)
    print()
    
    mn=min_chisq / dof
    
    print('Fitted parameters with 68% C.I.: ')
    
    for i, pmin in enumerate(popt):
        print('%2i %-10s %12f +/- %10f'%(i, name[i], pmin, np.sqrt(pcov[i,i])*np.sqrt(mn)))
        
    print()
    
    perr=np.sqrt(np.diag(pcov))
    print(perr)
    
    print("Correlation matrix")
    print("               ")
    
    for i in range(len(popt)): print('%-10s'%(name[i],)),
    print()
    
    for i in range(len(popt)):
        print('%10s'%(name[i])),
        for j in range(i+1):
            print('%10f'%(pcov[i,j]/np.sqrt(pcov[i,i]*pcov[j,j]),)),
        print()
    
    yfit=fx(tmodel,*popt)
    
    lambdas.append(popt[1])
    
    u_lambdas.append(np.sqrt(pcov[1,1])*np.sqrt(mn))
    
    ## Plotting-------------------------------------------------------------------------------------------------------------------------------------------

    plt.scatter(tdata,ydata,label="New Cases")
    
    plt.plot(tmodel,yfit,'-r',label="Exponential Curvefit") # Plots 
    
    
    plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside
    
    plt.xlabel('Days') # Labels x-axis

    plt.ylabel('New Cases') # Labels y-axis
    
    plt.legend()
    plt.show()
print(lambdas,'+/-',u_lambdas)

rates=[endpoints2[0],endpoints2[1],endpoints2[2]]


udata=[]
ydata=[]
tdata=[]
i=0
while i<len(rates):
    
    udata.append(u_lambdas[i])
    ydata.append(lambdas[i])
    tdata.append(rates[i])
    
    i=i+1    

def fxx(x, m, c):  
    return np.asarray(m)*x+c

m0 = 1
c0 = 0
p0 = [m0, c0]
name = ['m', 'c']



tmodel = np.linspace(0, 500, 1000)


ystart=fxx(tmodel,*p0)


## Curve Fit Function---------------------------------------------------------------------------------------------------------------------------------

popt,pcov=curve_fit(fxx,tdata,ydata,p0,sigma=udata,absolute_sigma=True)
dymin = (ydata-fxx(tdata,*popt))/udata
min_chisq = sum(dymin*dymin)
dof=len(tdata) - len(popt)

print(len(popt))
print('Chi square: ',min_chisq)
print('Number of degrees of freedom: ',dof)
print('Chi square per degree of freedom: ',min_chisq/dof)
print()

mn=min_chisq / dof

print('Fitted parameters with 68% C.I.: ')

for i, pmin in enumerate(popt):
    print('%2i %-10s %12f +/- %10f'%(i, name[i], pmin, np.sqrt(pcov[i,i])*np.sqrt(mn)))
    
print()

perr=np.sqrt(np.diag(pcov))
print(perr)

print("Correlation matrix")
print("               ")

for i in range(len(popt)): print('%-10s'%(name[i],)),
print()

for i in range(len(popt)):
    print('%10s'%(name[i])),
    for j in range(i+1):
        print('%10f'%(pcov[i,j]/np.sqrt(pcov[i,i]*pcov[j,j]),)),
    print()

yfit=fxx(tmodel,*popt)

grad=popt[0]
cval=popt[1]
u_grad=np.sqrt(pcov[0,0])*np.sqrt(mn)
u_cval=np.sqrt(pcov[1,1])*np.sqrt(mn)

## Plotting-------------------------------------------------------------------------------------------------------------------------------------------


plt.plot(tmodel,yfit,'-r',label="Linear Curvefit") # Plots 


plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside

plt.xlabel('Wave Number') # Labels x-axis

plt.ylabel('Infection Rate') # Labels y-axis

xmin,xmax,ymin,ymax=plt.axis([0,500,0.0416,0.0657])

#plt.errorbar(4.0,fxx(4,grad,cval),xerr = None, yerr = np.sqrt(u_grad**2+u_cval**2), fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Proposed Growth Rate')

plt.errorbar(rates, lambdas, xerr = None, yerr = u_lambdas, fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Infection Rate Data') # Plots errorbar

plt.legend()

plt.show()

#print(fxx(4,0.006953,0.041106))

udata=[]
ydata=[]
tdata=[]
i=0
while i<len(rates):
    
    udata.append(u_lambdas[i])
    ydata.append(lambdas[i])
    tdata.append(rates[i])
    
    i=i+1    

def fxxx(x, A, c):  
    return A*np.log(x)+c

A0 = 0.01
c0 = 0.04
#b0 = 0
p0 = [A0, c0]
name = ['A', 'c']



tmodel = np.linspace(96, 500, 1000)


ystart=fxxx(tmodel,*p0)


## Curve Fit Function---------------------------------------------------------------------------------------------------------------------------------

popt,pcov=curve_fit(fxxx,tdata,ydata,p0,sigma=udata,absolute_sigma=True)
dymin = (ydata-fxxx(tdata,*popt))/udata
min_chisq = sum(dymin*dymin)
dof=len(tdata) - len(popt)

print(len(popt))
print('Chi square: ',min_chisq)
print('Number of degrees of freedom: ',dof)
print('Chi square per degree of freedom: ',min_chisq/dof)
print()

mn=min_chisq / dof

print('Fitted parameters with 68% C.I.: ')

for i, pmin in enumerate(popt):
    print('%2i %-10s %12f +/- %10f'%(i, name[i], pmin, np.sqrt(pcov[i,i])*np.sqrt(mn)))
    
print()

perr=np.sqrt(np.diag(pcov))
print(perr)

print("Correlation matrix")
print("               ")

for i in range(len(popt)): print('%-10s'%(name[i],)),
print()

for i in range(len(popt)):
    print('%10s'%(name[i])),
    for j in range(i+1):
        print('%10f'%(pcov[i,j]/np.sqrt(pcov[i,i]*pcov[j,j]),)),
    print()

yfit=fxxx(tmodel,*popt)

grad=popt[0]
cval=popt[1]
#bval=popt[2]
u_grad=np.sqrt(pcov[0,0])*np.sqrt(mn)
u_cval=np.sqrt(pcov[1,1])*np.sqrt(mn)
#u_bval=np.sqrt(pcov[2,2])*np.sqrt(mn)

## Plotting-------------------------------------------------------------------------------------------------------------------------------------------


plt.plot(tmodel,yfit,'-r',label="Linear Curvefit") # Plots 

xmin,xmax,ymin,ymax = plt.axis([96,500,0.0447,0.066])

plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside

plt.xlabel('Wave Number') # Labels x-axis

plt.ylabel('Infection Rate') # Labels y-axis

#plt.errorbar(4.0,fxxx(4,grad,cval),xerr = None, yerr = np.sqrt(u_grad**2+u_cval**2), fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Proposed Growth Rate')

plt.errorbar(rates, lambdas, xerr = None, yerr = u_lambdas, fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Infection Rate Data') # Plots errorbar

plt.legend()

plt.show()

#print(fxxx(4,0.006953,0.041106))

#udata=[]
#ydata=[]
#tdata=[]
#i=0
#while i<len(rates):
    
    #udata.append(u_lambdas[i])
    #ydata.append(lambdas[i])
    #tdata.append(rates[i])
    
    #i=i+1    

#def fy(x, A, c):  
    #return A*np.sqrt(x)+c

#A0 = 1
#c0 = 0
##b0 = 1
#p0 = [A0, c0]
#name = ['A', 'c' ]



#tmodel = np.linspace(0.5, 500, 1000)


#ystart=fy(tmodel,*p0)


### Curve Fit Function---------------------------------------------------------------------------------------------------------------------------------

#popt,pcov=curve_fit(fy,tdata,ydata,p0,sigma=udata,absolute_sigma=True)
#dymin = (ydata-fy(tdata,*popt))/udata
#min_chisq = sum(dymin*dymin)
#dof=len(tdata) - len(popt)

#print(len(popt))
#print('Chi square: ',min_chisq)
#print('Number of degrees of freedom: ',dof)
#print('Chi square per degree of freedom: ',min_chisq/dof)
#print()

#mn=min_chisq / dof

#print('Fitted parameters with 68% C.I.: ')

#for i, pmin in enumerate(popt):
    #print('%2i %-10s %12f +/- %10f'%(i, name[i], pmin, np.sqrt(pcov[i,i])*np.sqrt(mn)))
    
#print()

#perr=np.sqrt(np.diag(pcov))
#print(perr)

#print("Correlation matrix")
#print("               ")

#for i in range(len(popt)): print('%-10s'%(name[i],)),
#print()

#for i in range(len(popt)):
    #print('%10s'%(name[i])),
    #for j in range(i+1):
        #print('%10f'%(pcov[i,j]/np.sqrt(pcov[i,i]*pcov[j,j]),)),
    #print()

#yfit=fy(tmodel,*popt)

#Aval=popt[0]
#cval=popt[1]
##bval=popt[2]
#u_A=np.sqrt(pcov[0,0])*np.sqrt(mn)
#u_cval=np.sqrt(pcov[1,1])*np.sqrt(mn)
##u_bval=np.sqrt(pcov[2,2])*np.sqrt(mn)

### Plotting-------------------------------------------------------------------------------------------------------------------------------------------


#plt.plot(tmodel,yfit,'-r',label="Linear Curvefit") # Plots 


#plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside

#plt.xlabel('Wave Number') # Labels x-axis

#plt.ylabel('Infection Rate') # Labels y-axis

#plt.errorbar(4.0,fy(4,Aval,cval),xerr = None, yerr = np.sqrt(u_A**2+u_cval**2), fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Proposed Growth Rate')

#plt.errorbar(rates, lambdas, xerr = None, yerr = u_lambdas, fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Infection Rate Data') # Plots errorbar

#plt.legend()

#plt.show()

#print(fy(4,Aval,cval))

udata=[]
ydata=[]
tdata=[]
i=0
while i<len(rates):
    
    udata.append(u_lambdas[i])
    ydata.append(lambdas[i])
    tdata.append(rates[i])
    
    i=i+1    

def fn(x, A, c):  
    return -A/(np.log(x+np.asarray(1)))+c

A0 = 0.01
c0 = 0.04
#b0 = 0
p0 = [A0, c0]
name = ['A', 'c']



tmodel = np.linspace(96, 500, 1000)


ystart=fn(tmodel,*p0)


## Curve Fit Function---------------------------------------------------------------------------------------------------------------------------------

popt,pcov=curve_fit(fn,tdata,ydata,p0,sigma=udata,absolute_sigma=True)
dymin = (ydata-fn(tdata,*popt))/udata
min_chisq = sum(dymin*dymin)
dof=len(tdata) - len(popt)

print(len(popt))
print('Chi square: ',min_chisq)
print('Number of degrees of freedom: ',dof)
print('Chi square per degree of freedom: ',min_chisq/dof)
print()

mn=min_chisq / dof

print('Fitted parameters with 68% C.I.: ')

for i, pmin in enumerate(popt):
    print('%2i %-10s %12f +/- %10f'%(i, name[i], pmin, np.sqrt(pcov[i,i])*np.sqrt(mn)))
    
print()

perr=np.sqrt(np.diag(pcov))
print(perr)

print("Correlation matrix")
print("               ")

for i in range(len(popt)): print('%-10s'%(name[i],)),
print()

for i in range(len(popt)):
    print('%10s'%(name[i])),
    for j in range(i+1):
        print('%10f'%(pcov[i,j]/np.sqrt(pcov[i,i]*pcov[j,j]),)),
    print()

yfit=fn(tmodel,*popt)

grad=popt[0]
cval=popt[1]
#bval=popt[2]
u_grad=np.sqrt(pcov[0,0])*np.sqrt(mn)
u_cval=np.sqrt(pcov[1,1])*np.sqrt(mn)
#u_bval=np.sqrt(pcov[2,2])*np.sqrt(mn)

## Plotting-------------------------------------------------------------------------------------------------------------------------------------------


plt.plot(tmodel,yfit,'-r',label="Linear Curvefit") # Plots 


plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside

plt.xlabel('Wave Number') # Labels x-axis

plt.ylabel('Infection Rate') # Labels y-axis

xmin,xmax,ymin,ymax = plt.axis([96,500,0.0447,0.066])

#plt.errorbar(4.0,fn(4,grad,cval),xerr = None, yerr = np.sqrt(((1/np.log(5))*u_grad)**2+u_cval**2), fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Proposed Growth Rate')

plt.errorbar(rates, lambdas, xerr = None, yerr = u_lambdas, fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Infection Rate Data') # Plots errorbar

plt.legend()

plt.show()





#udata=[]
#ydata=[]
#tdata=[]
#i=0
#while i<len(rates):
    
    #udata.append(u_lambdas[i])
    #ydata.append(lambdas[i])
    #tdata.append(rates[i])
    
    #i=i+1    

#def fn1(x, c):  
    #return -np.exp(np.asarray(-1)*x)+c

##A0 = 0.01
#c0 = 0.04
##b0 = 0
#p0 = [c0 ]
#name = ['c']



#tmodel = np.linspace(17, 500, 1000)


#ystart=fn1(tmodel,*p0)


### Curve Fit Function---------------------------------------------------------------------------------------------------------------------------------

#popt,pcov=curve_fit(fn1,tdata,ydata,p0,sigma=udata,absolute_sigma=True)
#dymin = (ydata-fn1(tdata,*popt))/udata
#min_chisq = sum(dymin*dymin)
#dof=len(tdata) - len(popt)

#print(len(popt))
#print('Chi square: ',min_chisq)
#print('Number of degrees of freedom: ',dof)
#print('Chi square per degree of freedom: ',min_chisq/dof)
#print()

#mn=min_chisq / dof

#print('Fitted parameters with 68% C.I.: ')

#for i, pmin in enumerate(popt):
    #print('%2i %-10s %12f +/- %10f'%(i, name[i], pmin, np.sqrt(pcov[i,i])*np.sqrt(mn)))
    
#print()

#perr=np.sqrt(np.diag(pcov))
#print(perr)

#print("Correlation matrix")
#print("               ")

#for i in range(len(popt)): print('%-10s'%(name[i],)),
#print()

#for i in range(len(popt)):
    #print('%10s'%(name[i])),
    #for j in range(i+1):
        #print('%10f'%(pcov[i,j]/np.sqrt(pcov[i,i]*pcov[j,j]),)),
    #print()

#yfit=fn1(tmodel,*popt)

##grad=popt[0]
##cval=popt[1]
###bval=popt[2]
##u_grad=np.sqrt(pcov[0,0])*np.sqrt(mn)
##u_cval=np.sqrt(pcov[1,1])*np.sqrt(mn)
###u_bval=np.sqrt(pcov[2,2])*np.sqrt(mn)

### Plotting-------------------------------------------------------------------------------------------------------------------------------------------


#plt.plot(tmodel,yfit,'-r',label="Linear Curvefit") # Plots 


#plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside

#plt.xlabel('Wave Number') # Labels x-axis

#plt.ylabel('Infection Rate') # Labels y-axis

##plt.errorbar(4.0,fn1(4,grad,cval),xerr = None, yerr = np.sqrt(((1/np.log(5))*u_grad)**2+u_cval**2), fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Proposed Growth Rate')

#plt.errorbar(rates, lambdas, xerr = None, yerr = u_lambdas, fmt = '', marker='_', ls = 'None',capsize=2.3, ecolor = 'b',label='Infection Rate Data') # Plots errorbar

#plt.legend()

#plt.show()