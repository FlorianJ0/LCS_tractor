import numpy as np
import tricubic

n = 3
f = np.zeros((n,n*2,n), dtype='float')

for i in range(n):
  for j in range(n*2):
    for k in range(n):
      f[i][j][k] = i+j+k #some function f(x,y,z) is given on a cubic grid indexed by i,j,k
      
ip = tricubic.tricubic(list(f), [n,2*n,n]) #initialize interpolator with input data on cubic grid

print f.shape, f.dtype




# print f

res=ip.ip(list(np.array([2,5,4])))
# for i in range(100000):
  # res = ip.ip(list(np.random.rand(3)*(n-1))) #interpolate the function f at a random point in space
print res

for i in xrange(10*n):
  ii=0.1*i
  iii=int(ii)
  print ip.ip(list(np.array([2,5,ii])))
  print f[2,5,iii]
  print '^^^'