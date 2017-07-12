import matplotlib.pyplot as plt
import numpy as np
import math
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
def mydependent(lamb, gstart, rstart, iterations, saveplace):
    rs = []
    gs = []
    rs2 = []
    gs2 = []
    mys = np.linspace(0, 3.7, 10000)
    for my in mys:
        g = gstart
        r = rstart
        gnew = rnew = 0
        for i in range(0, iterations):
            gnew = lamb*g*(1-g) - my*g*r
            rnew = (lamb + 0.2)*r*(1-r) - my*g*r
            g = gnew
            r = rnew
            if(i == iterations - 2):
                rs2.append(r)
                gs2.append(g)
        rs.append(r)
        gs.append(g)

    plt.gca().set_color_cycle(['red', 'grey'])
    plt.title("My")
    plt.plot(mys, rs)
    plt.plot(mys, gs)
    plt.plot(mys, rs2)
    plt.plot(mys, gs2)
    plt.savefig(saveplace + '.pdf')

    plt.legend(['Red Squirrels', 'Grey Squirrels'])
    plt.show()
def plotKillGrey(lambda_g, lambda_r, mu_g, mu_r, saveplace):
    d = np.linspace(0, 1, 10000)
    klambda_g = []
    for i in range(0, 10000):
        klambda_g.append(killsneeded(lambda_r + d[i], lambda_r, mu_g, mu_r))
    plt.grid()
    plt.plot(d, klambda_g)
    plt.xlabel('$\lambda_g - \lambda_r$')
    plt.ylabel('Fraction of grey squirrels that needs to be killed')
    plt.savefig(saveplace + '.pdf')
    plt.show()
def killsneeded(lambda_g, lambda_r, mu_g, mu_r):
    return (lambda_g - 1 - lambda_g*(lambda_r - 1)/mu_r)
def plotFixedPoints(lambda_g, lambda_r, mu_g, mu_r, saveplace):
    r = []
    g = []
    r.append(0)
    g.append(0)
    r.append(0)
    g.append(1-(1/lambda_g))
    r.append(1-(1/lambda_r))
    g.append(0)
    r.append(rzero(lambda_r, lambda_g, mu_g, mu_r))
    g.append(gzero(lambda_r, lambda_g, mu_g, mu_r))
    plt.plot(g, r, 'bo')
    plt.grid()
    plt.xlabel('$g$')
    plt.ylabel('$r$') 
    plt.savefig(saveplace + '.pdf')
    plt.show()
def lambdadependent(myr, myg, gstart, rstart, iterations, saveplace, theoretical):
    rs = []
    gs = []
    for i in range(0, 20):
        rs.append([])
        gs.append([])
    lambs = np.linspace(0.9, 3.0, 10000)
    for lamb in lambs:
        g = gstart
        r = rstart
        gnew = rnew = 0
        for i in range(0, iterations):
            gnew = (lamb + 0.25)*g*(1-g) - myg*g*r
            rnew = (lamb)*r*(1-r) - myr*g*r
            g = gnew
            r = rnew
            if(g < 0):
                g = 0.5
            if(r < 0):
                r = 0.5
        for i in range(0, 20):
            gnew = (lamb + 0.25)*g*(1-g) - myg*g*r
            rnew = (lamb)*r*(1-r) - myr*g*r
            g = gnew
            r = rnew
            if(g < 0):
                g = 0.5
            if(r < 0):
                r = 0.5
            gs[i].append(g)
            rs[i].append(r)

#    plt.gca().set_color_cycle(['red', 'grey'])
    for i in range(0, 20):
        red_line, = plt.plot(lambs, rs[i], 'r-', label='Red Squirrel')
        grey_line, = plt.plot(lambs, gs[i], ls='-.', c='grey', label='Grey Squirrel')
        
    if(theoretical):
        
        fixedr = []
        fixedg = []
        for lamb in lambs:
            if(rzero(lamb, lamb + 0.25, myg, myr) > 1 or rzero(lamb, lamb + 0.25, myg, myr) < -0.5):
                if(len(fixedr) == 0):
                    fixedr.append(-0.4)
                else:
                    fixedr.append(fixedr[len(fixedr) - 1])
            else:
                fixedr.append(rzero(lamb, lamb + 0.25, myg, myr))
        for lamb in lambs:
            if(gzero(lamb, lamb + 0.25, myg, myr) > 1 or gzero(lamb, lamb + 0.25, myg, myr) < -0.5):
                if(len(fixedg) == 0):
                    fixedg.append(-0.4)
                else:
                    fixedg.append(fixedg[len(fixedg) - 1])
            else:
                fixedg.append(gzero(lamb, lamb + 0.25, myg, myr))

        lilac_line, = plt.plot(lambs, fixedr, ls=':', c='magenta', label='Equation (5) Red Squirrels')
        black_line, = plt.plot(lambs, fixedg, ls='--', c='black', label='Equation (5) Grey Squirrel')

#    plt.title("Squirrel populations")
    if(theoretical):
        plt.legend(handles=[red_line, grey_line, lilac_line, black_line])
    else:
        plt.legend(handles=[red_line, grey_line])
    plt.xlabel("$\lambda_r$")
    plt.ylabel("Population")
    plt.grid()

    plt.savefig(saveplace + '.pdf')
    plt.show()

def rzero(lambdar, lambdag, myg, myr):
    return ((lambdar - 1)*lambdag - myr*(lambdag - 1))/(lambdar*lambdag - myr*myg)

def gzero(lambdar, lambdag, myg, myr):
    return ((lambdag - 1)*lambdar - myg*(lambdar - 1))/(lambdar*lambdag - myr*myg)

def diskrim(p, q):
    return ((p/2)**2 - q)

def alfa1(lambdar, r, myr, g):
    return lambdar*(1 - 2*r) - myr*g

def beta1(lambdag, g, myg, r):
    return lambdag*(1 - 2*g) - myg*r

def gamma1(r, g, myr, myg):
    return r*g*myr*myg

def eigenValues(lambdar, lambdag, myr, myg):
    r = rzero(lambdar, lambdag, myg, myr)
    g = gzero(lambdar, lambdag, myg, myr)
    alfa = alfa1(lambdar, r, myr, g)
    beta = beta1(lambdag, g, myg, r)
    gamma = gamma1(r, g, myr, myg)
    dis = diskrim(-(alfa + beta), alfa*beta - gamma)
    return ((alfa + beta)/2 + dis**0.5, (alfa + beta)/2 - dis**0.5) 

def plotEigenvalues(myg, myr, name, points, interval):
    a = np.zeros((points, points))
    lambdars2 = np.linspace(interval[0], interval[1], points)
    lambdags = np.linspace(interval[0], interval[1], points)
    lambdars = lambdars2[::-1]
    for i, lambdar in enumerate(lambdars):
        for j, lambdag in enumerate(lambdags):
            r = rzero(lambdar, lambdag, myg, myr)
            g = gzero(lambdar, lambdag, myg, myr)
            alfa = alfa1(lambdar, r, myr, g)
            beta = beta1(lambdag, g, myg, r)
            gamma = gamma1(r, g, myr, myg)
            dis = diskrim(-(alfa + beta), alfa*beta - gamma)
            
            if(dis < 0):
                a1 = abs(complex((alfa + beta)/2, (-dis)**0.5))
                a2 = abs(complex((alfa + beta)/2,-(-dis)**0.5))
                a[i][j] = 0
            else:
                a1 = abs((alfa + beta)/2 + (dis)**0.5)
                a2 = abs((alfa + beta)/2 - (dis)**0.5)
            if(a1 < 1 and a2 < 1):
                a[i][j] = 1
            elif((a1 < 1 and a2 > 1) or (a1 > 1 and a2 < 1)):
                a[i][j] = 2
            elif(a1 > 1 and a2 > 1):
                a[i][j] = 3
            else:
                a[i][j] = 0
                       
    plt.imshow(a, extent = [interval[0], interval[1], interval[0], interval[1]])

#    plt.title('Type of Point')
    plt.xlabel('$\lambda_g$')
    plt.ylabel('$\lambda_r$')
    plt.hot()
    plt.grid()
    plt.savefig(saveplace + '.pdf')
    plt.savefig(name)
    plt.show()

#plotEigenvalues(0.5, 0.5, 'sinks0,5_0,5', 1000, (1.0, 3.0))
#plotEigenvalues(1.5, 1.5, 'sink1,5_1,5', 1000, (1.0, 3.0))
#plotEigenvalues(0.5, 0.5, 'points0.5_0.5', 2000, (1.0, 3.0))
#print (eigenValues(1.5, 1.1, 1.5, 1.5))
#lambdadependent(1.5, 1.5, gzero(1.1, 1.35, 1.5, 1.5), rzero(1.1, 1.35, 1.5, 1.5), 100, 'test')
lambdadependent(0.5, 0.5, 0.5, 0.5, 100, "lambdamyT0,5_0,5", True)
lambdadependent(0.5, 0.5, 0.5, 0.5, 100, "lambdamyF0,5_0,5", False)
#lambdadependent(1.5, 1.5, 0.01, 0.9, 100, "lambdamyT0,9_0,9", True)
#lambdadependent(0.6, 0.8, 0.1, 0.5, 100, "lambdamy0.6_0.8")
#lambdadependent(0.7, 0.8, 0.1, 0.5, 100, "lambdamy0.7_0.8")
#lambdadependent(0.9, 1.0, 0.1, 0.5, 100, "lambdamy0.9_1.0")
#lambdadependent(0.8, 1.0, 0.1, 0.5, 100, "lambdamy0.8_1.0")
#lambdadependent(1.5, 1.5, 0.9, 0.1, 100, "lambdamy1.5_1.5")

#plotFixedPoints(1.75, 1.5, 0.5, 0.5, 'fp1,75_1,5')
#plotKillGrey(1.3, 1.05, 0.5, 0.5, 'killGrey1,05')   





