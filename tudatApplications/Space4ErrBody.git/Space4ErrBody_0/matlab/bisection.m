function root = bisection(f,a,b)

if f(a)*f(b)>0
    disp('Wrong choice bro')
else
    root = (a + b)/2;
    err = abs(f(root));
    while err > 1e-7
        if f(a)*f(root)<0
            b = root;
        else
            a = root;
        end
        root = (a + b)/2;
        err = abs(f(root));
    end
end