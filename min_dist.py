from sage.modular.dirichlet import DirichletCharacter
from sage.modular.dims import dimension_new_cusp_forms

set_verbose(-2)
gp.default("realprecision",20)

H = DirichletGroup(1, base_ring=CyclotomicField(1))
chi = DirichletCharacter(H, H._module([]))
pi = numerical_approx(pi, digits=30)

MM = 47
result = []
arr_upper = []
min_arr = []
for M in range(7, MM + 1):
    var('z')
    N = Newforms(chi, 2 * M + 2, names="a")
    dim = dimension_new_cusp_forms(Gamma0(1), 2 * M + 2)
    f = N[0]

    # Iterating over every cusp form
    for j in range(0, dim):
        # L is the L function associated with each new form
        L = f.lseries(j)

        # "r" is the period polynomial and is defined as in (2.6) from [Conrey]
        r = 0
        for m in range(M):
            r = r + (-1)**m * ((2 * pi * z)**(2 * m + 1) * L(2 * M - 2 * m) / factorial(2 * m + 1))
        
        if r != 0:
            arr = r.roots(ring=CC)
            arr_arg = []
            for x in arr:
                # rotate
                arr_arg.append(arg(x[0]))
            # print(M, arr_arg)
            arr_upper = []

            for x in arr_arg:
                if 10**(-6) <= x <= 3.14 - 10**(-6):
                    arr_upper.append(x)
            arr_upper.sort()
            print(M-5, len(arr_upper), arr_upper)

            result.append([M, arr_upper])

            min = pi
            for k in range(0,len(arr_upper)-1):
                if ((arr_upper[k+1]-arr_upper[k])<min):
                    min=arr_upper[k+1]-arr_upper[k];
            print(str(2*M+2) + ": "+ str(min));
            min_arr.append(min)

tmp = pi
for k in range(0,len(min_arr)-1):
    if (min_arr[k]<tmp):
        tmp = min_arr[k]
print(tmp)
