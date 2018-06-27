from math import floor

remove_spaces = lambda inlst: [i for i in inlst if i != ' ']

def fix_signs(inlst):
    i = 0
    while i < len(inlst):
        if inlst[i] in '+-':    # first sign is detected
            sign = -1 if inlst[i] == '-' else 1    # sign variable assigned
            while i+1 < len(inlst) and inlst[i+1] in '+-':    # while more signs are present
                if inlst[i+1] == '-':    # invert the sign if a minus is detected
                    sign *= -1
                del inlst[i+1]    # delete each excessive sign
            inlst[i] = '-' if sign == -1 else '+'    # change the only sign left's value accordingly
        i += 1    # keep checking for other signs
    return inlst

def fix_dict(indict):
    if type(indict) == dict:
        return frozenset(indict.items())
    else:
        return indict

def get_coefficient(inlst, i):
    coeff = ''
    k = i - 1
    while k >= 0 and inlst[k] in '1234567890.':    # keep going backwards to get the full coefficient
        coeff = inlst[k] + coeff
        k -= 1
    coeff = '1' if not coeff else coeff    # if no coefficient is specified, 1 is assigned
    if k >= 0 and inlst[k] == '-':    # check for a minus sign
        coeff = '-' + coeff
    k = 0 if k < 0 else k    # value correction for convert()
    coeff = float(coeff)
    return (coeff, k)

def get_exponent(inlst, i):
    exp = ''
    if i+1 < len(inlst) and inlst[i+1] == '^':
        k = i + 2
        while k < len(inlst) and inlst[k] in '1234567890':    # keep going forward to get the full exponent
            exp += inlst[k]
            k += 1
    else:
        k = i + 1    # value correction for convert()
    exp = 1 if not exp else exp    # if no exponent is specified, 1 is assigned
    exp = int(exp)    # exponents are assumed to be positive integers
    return (exp, k)

def convert(inlst):
    exps = {}
    i = 0
    while i < len(inlst):
        if inlst[i] == 'x':    # if an x-term is detected
            (coeff, x_start) = get_coefficient(inlst, i)    # get its coefficient
            (exp, x_end) = get_exponent(inlst, i)    # get its exponent
            if exp not in exps:
                exps[exp] = coeff
            else:
                exps[exp] += coeff
            del inlst[x_start:x_end]
            i = x_start
        i += 1
    return exps

def solve_x0_terms(inlst):
    out = 0
    current_term = ''
    while inlst:
        item = inlst.pop(0)
        if item in '#+-':
            out += float(current_term) if current_term else 0
            current_term = '-' if item == '-' else ''
        elif item in '1234567890.':
            current_term += item
    out += float(current_term) if current_term else 0
    return out

def divide_func(exps, div):    # uses polynomial long division
    newexps = {}
    for current_exp in range(max(exps)-max(div), -1, -1):
        if max(exps) - max(div) != current_exp:    # bugfix: FOR Loop coud be changed to something more efficient (needs testing with high exponents)
            continue
        newexps[current_exp] = exps[max(exps)] / div[max(div)]
        for exp, coeff in div.items():
            m_coeff = exp + current_exp
            if m_coeff not in exps:
                exps[m_coeff] = 0
            exps[m_coeff] -= (newexps[current_exp] * coeff)
            if exps[m_coeff] == 0:
                del exps[m_coeff]    # deletion required because of max() in the main loop that could return a coeff with value 0
    return newexps if not exps or not exps[0] else {}    # if there is a reminder, return an empty dict; could be changed to return reminder

def n_factors(n):
    if type(n) == float and not n.is_integer():
        yield (n, 1)
        raise StopIteration
    else:
        n = int(n)
        yield (n, 1)
    if n % 2 == 0:
        for i in range(floor(abs(n/2)), 0, -1):
            if n % i == 0:
                yield (i, int(n/i))
    else:
        tn = floor(abs(n/2))
        for i in range( (tn - 1 if tn % 2 == 0 else tn), 0, -2 ):
            if n % i == 0:
                yield (i, int(n/i))

def x2terms(exps):
    a = exps[2] if 2 in exps else 0
    b = exps[1] if 1 in exps else 0
    c = exps[0] if 0 in exps else 0
    return a,b,c

def delta_calc(a,b,c):
    return b**2 - 4*a*c

def pow_diff(poly):
    out = ()
    root1 = (poly[max(poly)]) ** (1.0 / max(poly))
    root2 = abs(poly[0]) ** (1.0 / max(poly))
    if root1.is_integer() and root2.is_integer():
        root1, root2 = int(root1), int(root2)
        if max(poly) % 2 == 0:
            if poly[0] < 0:
                xm, x0 = poly[max(poly)]**(1.0/2), abs(poly[0])**(1.0/2)
                out = (( { int(max(poly)/2):xm, 0:x0 }, 1 ), ( { int(max(poly)/2):xm, 0:-x0 }, 1 ))
        else:
            if poly[0] < 0:
                out = [( { 1:root1, 0:-root2}, 1 )]
            else:
                out = [( { 1:root1, 0:root2 }, 1 )]
    return out

def binomial_mult(poly, expsort):
    out = ()
    for x0t1, x0t2 in n_factors(poly[0]):
        for xmt1, xmt2 in n_factors(poly[expsort[0]]):
            if (xmt1*x0t2)+(xmt2*x0t1) == poly[expsort[1]]:
                p_div1 = { expsort[1]:xmt1, 0:x0t1 }
                p_div2 = { expsort[1]:xmt2, 0:x0t2 }
                out = (( p_div1, 1 ), ( p_div2, 1 ))
    return out

def factorize(poly_list, func):
    poly = poly_list.pop(0)
    tmexp = max(poly)
    div_polys = []
    if tmexp and poly[0] == 0:
        func.add({ 1:1, 0:0 }, min([e if e > 0 else tmexp for e in poly]))
        poly = divide_func(poly, { min([e if e > 0 else tmexp for e in poly]):1 })
        poly_list.append(poly)
    elif len(poly) == 2 and poly[0]:    # x^2 - 1 -> (x + 1)(x - 1), x^3 - 1, x^3 + 1, etc.
        div_polys = pow_diff(poly)
    elif len(poly) == 3 and poly[0]:    # x^2 + 2x + 1 -> (x + 1)^2, 3x^2 + 7x + 2 -> (3x + 1)(x + 2), etc. max exp can be > 2
        expsort = sorted(poly)[::-1]
        if expsort[0] % 2 == 0 and expsort[0]-expsort[1] == expsort[1]-expsort[2]:
            div_polys = binomial_mult(poly, expsort)
    for p, e in div_polys:
        for div_i in range(e):
            poly = divide_func(poly, p)
            if max(p) > 2 or (max(p) == 2 and p[0] and delta_calc(*x2terms(p)) >= 0):
                poly_list.append(p)
            else:
                func.add(p, 1)
    #elif: other polynomials
    if poly_list:
        factorize(poly_list, func)
    elif poly != {0:1}:
        func.add(poly, 1)

class frozendict_exps:
    def __init__(self):
        self.data = {}
    def __repr__(self):
        return repr(self.data)
    def add(self, indict, exp):
        self.indict = fix_dict(indict)
        if self.indict in self.data:
            self.data[self.indict] += exp
        else:
            self.data[self.indict] = exp
    def remove(self, indict, exp):
        self.indict = fix_dict(indict)
        del self.data[self.indict]
    def edit(self, indict, newdict, exp):
        self.remove(indict)
        self.add(newdict, exp)

class Function(frozendict_exps):
    def __init__(self, data):
        self.data = {}
        if type(data) == dict:
            self.exps = data
        else:
            self.eqt = remove_spaces(data)
            self.eqt = fix_signs(self.eqt)
            self.exps = convert(self.eqt)    # self.eqt is referenced and edited directly by convert()
            self.exps[0] = solve_x0_terms(self.eqt)    # x-terms have already been removed from self.eqt
        self.out = ""
    def __repr__(self):
        return repr(self.data)
