#eq.py
import time
start=time.time()
e=50
for i in range(1+e):
    a='*'*i
    b='.'*(50-i)
    c=2*i
    print("\r{}%{}-->{}\r".format(c,a,b))
    time.sleep(0.1)
end=time.time()
m=end-start
print("耗时{}s".format(m))
