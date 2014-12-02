def root_bisec(f,a,b,error=10**-7){
    '''Root of a function using the bisection method
    
    The aruments expected are:
    f: a function with the format double f(double x)
    a,b: bounds with a separated root.
    error: precision for the root determination default value 10^-7
    '''
    ext=f(a)*f(b)
    if ext>0:
        return nul
    elif ext==0:
        for t in [a,b]:
            if f(t)==0:
                return t
    else:
        while abs(a-b)>error:
            x=(a+b)/2
            ext=f(a)*f(x)
            if ext==0:
                return x
            elif ext>0:
                a=x
            else
                b=x
        return x
    }
    
def root_bisec(f,df,a,error=10**-7){
    '''Root of a function using the Newton-Raphson method
    
    The aruments expected are:
    f: a function with the format double f(double x)
    df: derivative of the function with the same format
    a: initial approximation of the root.
    error: precision for the root determination default value 10^-7
    '''
    ext=f(a)
    lx=a
    if ext==0:
        return a
    cond=True
    while cond:
        x=lx-f(lx)/df(lx)
        cond=abs(x-lx)
        lx=x
    return x
    }
