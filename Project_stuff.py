#Project Stuff

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

x = sym.symbols('x')

# %% 
####### Old Data
'''data = np.array([2282,2299,2275,2267,2279,2300,2265,2286,2245,
                 2217,2251,2201,2215,2225,2189,2210,2204,2177,
                 2193,2184,2210,2206,2200,2186,2182])

dates = np.array(['May 22-23','May 24-25','May 26-27','May 28-29','May 30-31','June 01-02',
                  'June 03-04','June 05-06','June 07-08','June 09-10',
                  'June 11-12','June 13-14','June 15-16','June 17-18','June 19-20',
                  'June 21-22','June 23-24','June 25-26','June 27-28','June 29-30',
                  'July 01-02','July 03-04','July 05-06','July 07-08','July 09-10'])'''


####### New Data
data = np.array([2302,2274,2272,2296,2265,2299,2267,2270,2234,2215,
                 2225,2195,2205,2186,2196,2194,2210,2185,2188])

dates = np.array(['May 16-18','May 19-21','May 22-24','May 25-27','May 28-30',
                  '31-June 02','June 03-05','June 06-08','June 09-11',
                  'June 12-14','June 15-17','June 18-20','June 21-23',
                  'June 24-26','June 27-29','30-July 02','July 03-05',
                  'July 06-08','July 09-11'])

m = len(data)

# %%
plt.grid(True, ls = '--', color = 'gray', axis = 'x')
plt.scatter(np.arange(0,len(data),1),data,c = 'navy')
plt.title('League-Wide Spin Rate Around Enforcement of Foreign Substance Rules')
plt.annotate('June 03: Rumors Around MLB of Enforcement Imminent',xy = (6,2265),
             xytext=(6.5,2290),fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.annotate('June 15: Official Announcement from MLB', xy = (10,2225),
             xytext=(10.5,2260),fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.annotate('June 21: Enforcement Begins',xy=(12,2205),xytext=(12.5,2240),
             fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.xticks(np.arange(0,len(data)),dates)
plt.xticks(fontsize=8, rotation=60)
plt.xlim(-1,m+1)
plt.ylim(2160,2320)
plt.ylabel('Spin Rate (rotations per minute)')
plt.xlabel('Date (2021)')
plt.savefig("scatter.svg", format = 'svg', dpi=300)
plt.show()

#%% Divided Differences
def myDD(xs,ys):
    x = sym.symbols('x')
    n = len(xs)
    
    C = np.zeros(shape=(len(xs),len(xs)))
                 
    C[:,0]=ys
    
    for j in np.arange(1,np.size(xs)):
        for i in np.arange(0,np.size(xs)-j):
            C[i,j]=(C[i+1,j-1]-C[i,j-1])/(xs[i+j]-xs[i])
            
    
    p = 0*x
    for i in np.arange(0,len(xs)):
        prod = x**0
        for j in np.arange(0,i):
            prod = prod*(x-xs[j])
        
        p = p + C[0,i]*prod
        
        
    for i in np.arange(0,len(xs)):
        for j in np.arange(0,len(xs)):
            if i+j>=len(xs):
                C[i,j]=np.nan
                
    return C, p

poly = myDD(np.arange(0,len(data),1),data)[1]
X = np.linspace(0,m,500)
Y = [poly.subs(x,val) for val in X]

print('Min:',min(Y),'Max:',max(Y))

plt.plot(X,Y)
plt.scatter(np.arange(0,m,1),data,c = 'navy')
plt.grid(True, ls = '--', color = 'gray', axis = 'x')
plt.xlim(-1,m+1)
plt.ylim(2160,2320)
plt.title('League-Wide Spin Rate Around Enforcement of Foreign Substance Rules')
plt.xticks(np.arange(0,m),dates)
plt.xticks(fontsize=8, rotation=60)
plt.ylabel('Spin Rate (rotations per minute)')
plt.xlabel('Date')
plt.savefig("poly.svg", format = 'svg', dpi=300)
plt.show()

# %% Romberg Integration

def myRomberg(f,I,M):
    x = sym.symbols('x')
    
    a = I[0]
    b = I[1]
    
    R0 = np.array([])
    
    for i in np.arange(0,M+1,1):
        xs = np.linspace(a,b,(2**i-1)+2)
        
        h = (b-a)/(2**i)
    
        ys = np.array([])
        for j in np.arange(0,len(xs)):
            if f.subs(x,xs[j])==sym.nan:
                ys = np.append(ys,sym.limit(f,x,xs[j]))
            if f.subs(x,xs[j])!=sym.nan:
                ys = np.append(ys,f.subs(x,xs[j]))
            
        foo = 0
        goo = 0
    
        for k in np.arange(1,len(xs)-1):
            goo = goo+ys[k]
            
        foo = (ys[0]+ys[-1]+2*goo)*h/2
        
        R0 = np.append(R0,foo)
        
    R = np.zeros(shape = (M+1,M+1))
    R[:,0]=R0
    
    for j in np.arange(1,M+1,1):
        for i in np.arange(j,M+1,1):
            R[i,j]=(4**j)/(4**j-1)*R[i,j-1]-(R[i-1,j-1])/(4**j-1)
            
    for i in np.arange(0,M+1,1):
        for j in np.arange(0,M+1,1):
            if j > i:
                R[i,j]=np.nan
            
    return R

area_romberg = myRomberg(poly,np.array([0,m]),7)[-1,-1]
print(area_romberg)
print(area_romberg/len(data))

area_1r = myRomberg(poly,np.array([3,(m-1)/2]),7)[-1,-1]
area_2r = myRomberg(poly,np.array([(m-1)/2,(m-1)-3]),7)[-1,-1]

print(area_1r/((m-1)/2-3))
print(area_2r/((m-1)/2-3))

perc_diff_rom = 1-(area_2r/((m-1)/2-3))/(area_1r/((m-1)/2-3))
print(perc_diff_rom)

plt.plot(X,Y)
plt.scatter(np.arange(0,m,1),data,c = 'navy')
plt.grid(True, ls = '--', color = 'gray', axis = 'x')
plt.xlim(-1,m+1)
plt.ylim(2000,2400)
plt.title('League-Wide Spin Rate Around Enforcement of Foreign Substance Rules')
plt.xticks(np.arange(0,m),dates)
plt.xticks(fontsize=8, rotation=60)
plt.ylabel('Spin Rate (rotations per minute)')
plt.xlabel('Date')

plt.fill_between(
        x= X.astype(np.float64), 
        y1= np.asarray(Y).astype(np.float64), 
        where= (0 < X)&(X < 18),
        color= "b",
        alpha= 0.2)

plt.annotate('Area = 491,415.078',xy=(5,2100),xytext=(5,2100),fontsize=12)
plt.annotate('Area per unit = 25,863.951', xy=(4,2070), xytext = (4,2070), fontsize = 12)

plt.savefig("poly_area.svg", format = 'svg', dpi=300)
plt.show()

# %% ############# Cubic Spline ########################

#input myCubicSpline(np.arange(0,len(xs),1),data)
def myCubicSpline(xs,ys):
    n = len(xs)
    
    x,a,b,c,d = sym.symbols('x a b c d')
    a = sym.symbols('a:%d'%(n-1)) #creates a0,a1,...,an-1
    b = sym.symbols('b:%d'%(n-1))
    c = sym.symbols('c:%d'%(n-1))
    d = sym.symbols('d:%d'%(n-1))
    
    equations = np.array([])
    sn_arr = np.array([])
    for i in np.arange(0,n-1,1):
        sn_arr = np.append(sn_arr,a[i]+b[i]*x+c[i]*x**2+d[i]*x**3)
    
    for i in np.arange(0,len(xs)-1):
        initl = (sn_arr[i]-ys[i]).subs(x,xs[i])
        conns = (sn_arr[i]-ys[i+1]).subs(x,xs[i+1])
        
        if i<(len(xs)-2):
            smooth = sym.diff(sn_arr[i],x,1).subs(x,xs[i+1])-sym.diff(sn_arr[i+1],x,1).subs(x,xs[i+1])
            doub_smooth = sym.diff(sn_arr[i],x,2).subs(x,xs[i+1])-sym.diff(sn_arr[i+1],x,2).subs(x,xs[i+1])
            
        equations = np.append(equations,[initl,conns,smooth,doub_smooth])    
        
    #Normal equations
    equations = np.append(equations, sym.diff(sn_arr[0],x,2).subs(x,xs[0]))
    equations = np.append(equations, sym.diff(sn_arr[-1],x,2).subs(x,xs[-1]))
    
    Sols = sym.solve(equations,dict = True)
    
    sol_arr = np.array([])
    for i in np.arange(0,n-1,1):
        sol_arr = np.append(sol_arr,Sols[0][a[i]]+Sols[0][b[i]]*x+Sols[0][c[i]]*x**2+Sols[0][d[i]]*x**3)
    
    return sol_arr

spline_funcs = myCubicSpline(np.arange(0,m,1),data)

# %% Plotting spline functions
x = sym.symbols('x')
xs = np.arange(0,len(data))

func = 0*x

func = sym.Piecewise((func,False),(spline_funcs[0],x<=xs[1]))

for i in np.arange(1,len(xs)-2,1):
    func = sym.Piecewise((func,x<=xs[i]),(spline_funcs[i],x<=xs[i]+1))
    
func = sym.Piecewise((func,x<=xs[-1]-1),(spline_funcs[-1],x>xs[-1]-1))

X = np.linspace(0,m-1,500)
Y = [func.subs(x,val) for val in X]

plt.plot(X,Y)
plt.xlim(-1,m+1)
plt.ylim(2160,2320)
plt.grid(True, ls = '--', color = 'gray', axis = 'x')
plt.title('League-Wide Spin Rate Around Enforcement of Foreign Substance Rules')
'''plt.annotate('June 03: Rumors Around MLB of Enforcement Imminent',xy = (1,2265),
             xytext=(1.5,2290),fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.annotate('June 15: Official Announcement from MLB', xy = (7,2215),
             xytext=(7.5,2260),fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.annotate('June 21: Enforcement Begins',xy=(10,2210),xytext=(10.5,2240),
             fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))'''

plt.xticks(np.arange(0,m),dates)
plt.xticks(fontsize=8, rotation=60)
plt.ylabel('Spin Rate (rotations per minute)')
plt.xlabel('Date')
plt.scatter(np.arange(0,len(data),1),data,c = 'navy')
plt.savefig("spline.svg", format = 'svg', dpi=300)
plt.show()

# %% Close-up
plt.plot(X,Y)
plt.xlim(9.9,10.1)
plt.ylim(2210,2235)
plt.grid(True, ls = '--', color = 'gray', axis = 'x')
plt.title('League-Wide Spin Rate Around Enforcement of Foreign Substance Rules')

plt.xticks(np.arange(10,11),dates[10:11])
#plt.xticks(fontsize=8, rotation=60)
plt.ylabel('Spin Rate (rotations per minute)')
plt.xlabel('Date')
plt.scatter(np.arange(0,len(data),1),data,c = 'navy')
plt.savefig("smoothness.svg", format = 'svg', dpi=300)
plt.show()
# %% Composite Trapezoidal Rule

def CTR(f,I,n):
    x = sym.symbols('x')
    
    a = I[0]
    b = I[1]
    
    h = (b-a)/n
    
    xs = np.linspace(a,b,n+1)
    
    ys = np.array([])
    for j in np.arange(0,len(xs)):
        ys = np.append(ys,f.subs(x,xs[j]))
    
    foo = 0
    goo = 0
    
    for k in np.arange(1,len(xs)-1):
        goo = goo+ys[k]
    
    foo = (ys[0]+ys[-1]+2*goo)*h/2
    
    return foo

n_s = 100

area_spline = CTR(func,np.array([0,m]),n_s)
print(area_spline)
print(area_spline/len(data))

area_1_spline = CTR(func,np.array([0,m/2]),n_s)
#print(area_1_spline/10)
area_2_spline = CTR(func,np.array([m/2,m]),n_s)
#print(area_2_spline/9)

area_1_spline_norm = area_1_spline/(m/2)
area_2_spline_norm = area_2_spline/(m/2)

perc_dec = 1-(area_2_spline_norm/area_1_spline_norm)
print(perc_dec)

plt.plot(X,Y)
plt.scatter(np.arange(0,m,1),data,c = 'navy')
plt.grid(True, ls = '--', color = 'gray', axis = 'x')
plt.xlim(-1,m+1)
plt.ylim(2050,2350)
plt.title('League-Wide Spin Rate Around Enforcement of Foreign Substance Rules')
plt.xticks(np.arange(0,m),dates)
plt.xticks(fontsize=8, rotation=60)
plt.ylabel('Spin Rate (rotations per minute)')
plt.xlabel('Date')

plt.fill_between(
        x= X.astype(np.float64), 
        y1= np.asarray(Y).astype(np.float64), 
        where= (0 < X)&(X < 18),
        color= "b",
        alpha= 0.2)

plt.annotate('Area = 42,421.473',xy=(5,2150),xytext=(5,2150),fontsize=12)
plt.annotate('Area per unit = 2,232.709', xy=(4,2100), xytext = (4,2100), fontsize = 12)

plt.savefig("ctr_area.svg", format = 'svg', dpi=300)
plt.show()

#%% Analysis Plots

plt.plot(X,Y)
plt.xlim(-1,m+1)
plt.ylim(2100,2310)
plt.title('League-Wide Spin Rate Around Enforcement of Foreign Substance Rules')
'''plt.annotate('June 03: Rumors Around MLB of Enforcement Imminent',xy = (1,2265),
             xytext=(1.5,2290),fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.annotate('June 15: Official Announcement from MLB', xy = (7,2215),
             xytext=(7.5,2260),fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.annotate('June 21: Enforcement Begins',xy=(10,2210),xytext=(10.5,2240),
             fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))'''

plt.xticks(np.arange(0,len(data)),dates)
plt.xticks(fontsize=8, rotation=60)
plt.ylabel('Spin Rate (rotations per minute)')
plt.xlabel('Date')
plt.scatter(np.arange(0,len(data),1),data,c = 'navy')

plt.fill_between(
        x= X.astype(np.float64), 
        y1= np.asarray(Y).astype(np.float64), 
        where= (0 < X)&(X < 9),
        color= "b",
        alpha= 0.2)

plt.fill_between(
        x= X.astype(np.float64), 
        y1= np.asarray(Y).astype(np.float64), 
        where= (9 < X)&(X < 18),
        color= "r",
        alpha= 0.2)

plt.savefig("area_pre.svg", format = 'svg', dpi=300)

plt.annotate('2267.529',xy=(2.25,2200),xytext=(2.25,2200),fontsize=12)
plt.annotate('2197.886',xy=(11.5,2150),xytext=(11.5,2150),fontsize=12)

plt.savefig("area_post.svg", format = 'svg', dpi=300)

plt.show()

#%% Differential Equation

def myRKMO4(x_prime,t0,x0,I,h):
    a, b = I[0], I[1]
    N = int(abs((b-a)/h))
    #print(N)
    
    #need to make sure x0 corresponds with correct t
    if t0 == a:
        ts = np.linspace(a,b,N+1)
        
    if t0 == b:
        ts = np.linspace(b,a,N+1)
    
    #print(ts)
    ws = np.array([])
    ws = np.append(ws,x0)

    for i in np.arange(1,N+1,1):
        k1 = h*x_prime(ts[i-1],ws[i-1])
        k2 = h*x_prime(ts[i-1]+h/2 , ws[i-1]+k1/2)
        k3 = h*x_prime(ts[i-1]+h/2 , ws[i-1]+k2/2)
        k4 = h*x_prime(ts[i-1]+h , ws[i-1]+k3)
        new_w = ws[i-1]+(1/6)*(k1+2*k2+2*k3+k4)
        ws = np.append(ws,new_w)
    
    if t0 == a:
        return ws
        
    if t0 == b:
        return np.flip(ws) #reverse order

RMSE = np.array([])
intervs = np.arange(10,-10,-0.1)

for r in intervs:
    
    x = sym.symbols('x')
    t = sym.symbols('t')
    
    pp = lambda t,p: r*p*(1-(p/(max(data)+10)))*(1-(p/(min(data)-10)))
    Ivl = np.array([0,len(data)])
    h = 0.1
    
    ans = myRKMO4(pp,0,data[0],Ivl,h)
    check_ans = ans[0::10][0:-1]
    
    RMSE = np.append(RMSE,np.sqrt(np.mean((check_ans-data)**2)))
    
opt_r = np.argmin(RMSE)
pp = lambda t,p: intervs[opt_r]*p*(1-(p/(max(data)+10)))*(1-(p/(min(data)-10)))
ans = myRKMO4(pp,0,data[0],Ivl,h)
    
    

plt.grid(True, ls = '--', color = 'gray', axis = 'x')
plt.scatter(np.arange(0,len(data),1),data,c = 'navy')
plt.title('League-Wide Spin Rate Around Enforcement of Foreign Substance Rules')
plt.annotate('June 03: Rumors Around MLB of Enforcement Imminent',xy = (6,2265),
             xytext=(6.5,2290),fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.annotate('June 15: Official Announcement from MLB', xy = (10,2225),
             xytext=(10.5,2260),fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.annotate('June 21: Enforcement Begins',xy=(12,2205),xytext=(12.5,2240),
             fontsize = 8, arrowprops=dict(facecolor='black', width=0.1, 
                                                              headwidth=1,shrink=0.05))
plt.xticks(np.arange(0,len(data)),dates)
plt.xticks(fontsize=8, rotation=60)
plt.xlim(-1,m+1)
plt.ylim(2160,2320)
plt.ylabel('Spin Rate (rotations per minute)')
plt.xlabel('Date (2021)')
plt.plot(np.linspace(0,len(data),int(len(data)/h)+1),ans)
plt.savefig("logistic.svg", format = 'svg', dpi=300)
plt.show()