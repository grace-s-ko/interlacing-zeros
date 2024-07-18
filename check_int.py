from sage.modular.dirichlet import DirichletCharacter
from sage.modular.dims import dimension_new_cusp_forms

set_verbose(-2)
gp.default("realprecision", 20)

H = DirichletGroup(1, base_ring=CyclotomicField(1))
chi = DirichletCharacter(H, H._module([]))
pi = numerical_approx(pi, digits=30)

def zeros(NN, MM):
    result = []
    for M in range(NN, MM + 1):
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
                print(M)

                
                result.append([M, arr_upper])
    return result

k_zeros = zeros(7,39)
kp_zeros = zeros(7,44)

def check_int(k_arr, kp_arr):
  interlacing = True
  if(k_arr[0] < kp_arr[0]) :
    for i in range(1, len(k_arr[1])-2):
        for j in range(1, len(kp_arr[1])-2):
            if (k_arr[1][i] > kp_arr[1][j]) & (k_arr[1][i+1] < kp_arr[1][j+1]):
                interlacing = False
                print('False at ', i, j, k_arr[1][i], kp_arr[1][j], k_arr[1][i+1], kp_arr[1][j+1])
    return interlacing
  
def check_str_int(k_arr, kp_arr):
  interlacing = True
  if(k_arr[0] < kp_arr[0]) :
    print('At level', k_arr[0], kp_arr[0])
    if (k_arr[1][0] <= kp_arr[1][0]) & (k_arr[1][len(k_arr[1])-1] >= kp_arr[1][len(kp_arr[1])-1]):
        interlacing = False
    for i in range(1, len(k_arr[1])-2):
        for j in range(1, len(kp_arr[1])-2):
            if (k_arr[1][i] > kp_arr[1][j]) & (k_arr[1][i+1] < kp_arr[1][j+1]):
                interlacing = False
                print('False at ', i, j, k_arr[1][i], kp_arr[1][j], k_arr[1][i+1], kp_arr[1][j+1])
    return interlacing

for i in range(len(k_zeros)):
  for j in range(len(kp_zeros)):
    if check_str_int(k_zeros[i], kp_zeros[j]) == False:
        print('False at ', i, j)
    
