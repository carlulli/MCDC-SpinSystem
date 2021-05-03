import numpy as np

#input array, data set
#return value, mean value 
def mean(x):
    summ = 0
    for i in x:
        summ += i
    return summ/len(x)

#input array, data set
#return value, variance
def var(x):
    summ = 0
    aver_x = mean(x)
    for i in x:
        summ += (i - aver_x)**2

    return 1/((len(x) - 1))*summ

#input array --> data set
#return value --> variance of the mean
def mean_var(x):
    
    return var(x)/len(x)

#input array, data set
#return value, chisquared test
def chisquared(x):
    aver = mean(x)
    N = len(x)
    chisq = np.array([ (x[k] - aver)*(x[k] -aver)/aver for k in range(N)])
    return np.sum(chisq)


#input array,value --> data set,maximum correlation time
#return array --> autocorrelation function
def autocorrelation(x,t):
    N = int(len(x))
    t = int(t)
    aver = mean(x)
    autocorr = np.array([ (x[t+k] - aver)*(x[k] - aver) for k in range(N-t)])
    var = np.sum(autocorr)
    return var/(N-t)  

#input array,value --> data set, autocorrelation length
#return value --> integer integration time
def integ_time(func,tmax):
    N  = len(func)
    return 0.5 + np.array([autocorrelation(func,t)/autocorrelation(func,0) for t in range(tmax)]).sum()

#input array,value --> data set, autocorrelation time
#return value --> mean variance when accounting for autocorrelation
def mean_var_autocorr(func,tmax):
    if tmax<0:
        print("Error, Choose a positive tau\n")
        return -1
    N  = len(func)
    ac_t = integ_time(func,tmax)
    mvar_ac = mean_var(func)*2*ac_t
    return mvar_ac

#input array, value --> data set, size of bin
#return array --> binned data
def binning (x,m):
    if len(x)%m != 0:
        print("Error, choose different bin size")
        return -1
    num_bins = int(len(x)/m)
    bins = np.zeros(num_bins)
    for k in range(num_bins):
        bins[k] = 1/m*np.array([x[n]for n in range(m*k, m*(k+1))]).sum()
    return bins

#input array, value --> data, nu = len(data)/m where m is the size of the bin
#return array --> rebinned array with JK method
def JK_binning(x,nu):
    if isinstance(nu, int) == False:
        print("wrong nu choose appropriate integer\n")
        return -1
    jk_bin = np.zeros(nu)
    for k in range(nu):
        jk_bin[k] = 1/(nu-1)*(np.array([i  for i in x]).sum() - x[k])
    return jk_bin

#input array --> data set
#return value --> variance for the JK binned data
def Jk_var(x):
    nu = len(x)
    return (nu-1)/nu*var(x)

#input array,value,value -->  dataset, length of bin array, length of bin
def SBS_binning(x,k,mu):
    nu = len(x)
    x_bin = np.zeros(range(k))
    x_bin_temp = np.zeros(range(mu))
    for l in range(k):
        
        for m in range(mu):
            rand = np.random.randrange(0,1)
            index = np.ceil(rand*nu)
            x_bin_temp[m] = np.array([x[index]])
        
        x_bin[l] = 1/mu*np.sum(x_bin_temp)
      
    return x_bin



