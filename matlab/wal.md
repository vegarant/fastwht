# WAL
---
The Walsh function
## Syntax
---
```
y = wal(N,n,k)
y = wal(N,n,t)
```

## Description 
---
y = wal(N,n,k). Here we assume N = 2^r for some positive integer r, and that 0<= n,k < N are integers.
   The function evaluates the [WAL function](walsh_function.pdf) at 
    WAL(n,k/N). If n and k are vector input, it is evaluated at all the corresponding values. 

y = wal(N,n,t). Here we assume N = 2^r for some positive integer r, that n is an integer in the interval `[0,N)` and that lies i the interval `[0,1)`.
   The function evaluates the [WAL function](walsh_function.pdf) at 
    WAL(n,t). If n and k are vector input, it is evaluated at all the corresponding values. 

## Examples
---
```
>> N = 4;
>> n = 0:N-1;
>> k = n;
>> wal(N,n,k)
ans =
     1     1     1     1
     1     1    -1    -1
     1    -1    -1     1
     1    -1     1    -1
>> wal(N,n,3)
ans =
     1
    -1
     1
    -1
>> 
>> eps = 1e-14;
>> t = linspace(0,1-eps,N);
>> wal(N,3,t)
ans =
     1    -1     1    -1
>> wal(N,2,0.2) 
ans =
     1
```

