import pylab as pp
import numpy as np
from scipy import integrate, interpolate
from scipy import optimize

ar = open("SocialMediaUsers.csv").read().split("\n")[1:]
ar = ar[:len(ar)-1]

twitter = []
facebook = []
snapchat = []
xSnapchat = []
c = 0

for i,a in enumerate(ar):
  arI = a.split(",")[1:]
  twitter.append(float(arI[0])/7401)
  facebook.append(float(arI[1])/7401)

  if arI[2] != "0":
      xSnapchat.append(c)
      snapchat.append(float(arI[2])/7401)
      c = c + 1

x = np.arange(0, len(ar))

x_data = x
y_data = np.array(twitter)

def f(y, t, k):
  
    s = y[0]
    i = y[1]
    r = y[2]
    d = y[3]

    alpha= k[0]
    beta = k[1]
    gamma = k[2]
    theta = k[3]
    epsilon = k[4]
    
    
    dsdt = -alpha * s + epsilon * i
    didt = -beta * i + alpha * s + theta * r - epsilon * i
    drdt = beta * i - theta * r - gamma * r
    dddt = gamma * r

##    beta = k[0]
##    delta = k[1]
##    mu = k[2]
##
##    dsdt = -(beta * i * s) / (s + i + r + d)
##    didt = (beta * i * s) / (s + i + r + d) - delta * i - mu * i
##    drdt = delta * i
##    dddt = mu * i
    
    return [dsdt, didt, drdt, dddt]

def SolvedFunc(x, paramAr):
    funcLambda = lambda y, t: f(y, t, paramAr)
    solved = integrate.odeint(funcLambda ,y0 ,x)
    return solved[:,1]

def SolvedFunc0(x, paramAr):
    funcLambda = lambda y, t: f(y, t, paramAr)
    solved = integrate.odeint(funcLambda ,y0 ,x)
    return solved[:,0]

def SolvedFunc2(x, paramAr):
    funcLambda = lambda y, t: f(y, t, paramAr)
    solved = integrate.odeint(funcLambda ,y0 ,x)
    return solved[:,2]

def SolvedFunc3(x, paramAr):
    funcLambda = lambda y, t: f(y, t, paramAr)
    solved = integrate.odeint(funcLambda ,y0 ,x)
    return solved[:,3]

def Loss(p):
    return y_data - SolvedFunc(x_data,p)
  

guess = [0.01,0.01,0.01, 0.01, 0.01]
y0 = [1-y_data[0],y_data[0],0,0]


print(y0[0] + y0[1])

(params,kvg) = optimize.leastsq(Loss, guess, args=(),
                           Dfun=None, full_output=0, col_deriv=0,
                           ftol=1.49012e-08, xtol=1.49012e-08, gtol=0.0,
                           maxfev=0, epsfcn=None, factor=100, diag=None)

print(params)

tAxis = np.linspace(min(x_data), max(x_data),30) 
iAxis = interpolate.UnivariateSpline(tAxis, SolvedFunc(tAxis ,params), k=3, s=0)
sAxis = interpolate.UnivariateSpline(tAxis, SolvedFunc0(tAxis ,params), k=3, s=0)
rAxis = interpolate.UnivariateSpline(tAxis, SolvedFunc2(tAxis ,params), k=3, s=0)
dAxis = interpolate.UnivariateSpline(tAxis, SolvedFunc3(tAxis ,params), k=3, s=0)



tAxis = np.linspace(min(x_data), max(x_data),200)

pp.plot(x_data, y_data,'.r',tAxis,iAxis(tAxis),"-g",tAxis,sAxis(tAxis),'-b',tAxis,rAxis(tAxis),'-y',tAxis,dAxis(tAxis),'-r')

pp.xlabel('Time (Quarter of a year)',{"fontsize":16})
pp.ylabel("Users (Millions)",{"fontsize":16})
pp.legend(('data','fit',"s","r","d"),loc=0)
pp.suptitle("Twitter SIRD")
cache = "Parameters: " + str(params)
pp.title(cache, size = 6)
pp.savefig("SIRDtwitter")
pp.show()
