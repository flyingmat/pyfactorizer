from math import floor

remove_spaces = lambda inlst: [i for i in inlst if i != ' ']

def sf2i(inp):
    if float(inp).is_integer():
        return str(int(inp))
    else:
        return str(inp)

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
    if 0 not in newexps:
        newexps[0] = 0
    return newexps if not exps or not exps[0] else {}    # if there is a reminder, return an empty dict; could be changed to return reminder

def n_factors(n):
    if type(n) == float and not n.is_integer():
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
    if max(poly) % 2 == 0:
        root_exp = (1.0 / 2)
    else:
        root_exp = (1.0 / max(poly))
    root1 = (abs(poly[max(poly)]) ** root_exp) * (-1 if poly[max(poly)] < 0 else 1)
    root2 = (abs(poly[0]) ** root_exp) * (-1 if poly[0] < 0 else 1)
    if root1.is_integer() and root2.is_integer():
        root1, root2 = int(root1), int(root2)
        if max(poly) % 2 == 0:
            if poly[0]*poly[max(poly)] < 0:
                xm, x0 = root1, root2
                out = (( { int(max(poly)/2):xm, 0:x0 }, 1 ), ( { int(max(poly)/2):(xm if xm > 0 else -xm), 0:(x0 if xm < 0 else -x0) }, 1 ))
        else:
            out = [( { 1:root1, 0:root2}, 1 )]
    return out

def binomial_mult_3(poly, expsort):
    out = ()
    for x0t1, x0t2 in n_factors(poly[0]):
        for xmt1, xmt2 in n_factors(poly[expsort[0]]):
            if (xmt1*x0t2)+(xmt2*x0t1) == poly[expsort[1]]:
                p_div1 = { expsort[1]:xmt1, 0:x0t1 }
                p_div2 = { expsort[1]:xmt2, 0:x0t2 }
                out = (( p_div1, 1 ), ( p_div2, 1 ))
    return out

def binomial_pow3(poly, expsort):
    out = ()
    if expsort[0] % 3 == 0:
        root1 = (abs(poly[expsort[0]]) ** (1.0/3)) * (-1 if poly[expsort[0]] < 0 else 1)
        root2 = (abs(poly[0]) ** (1.0/3)) * (-1 if poly[0] < 0 else 1)
        if root1.is_integer() and root2.is_integer():
            if poly[expsort[1]] == 3*(root1**2)*root2 and poly[expsort[2]] == 3*(root2**2)*root1:
                out = [({ expsort[2]:root1, 0:root2 }, 3)]
    return out

def binomial_mult_4(poly, expsort):
    out = ()
    if poly[expsort[0]] / poly[expsort[2]] == poly[expsort[1]] / poly[expsort[3]]:
        cfs = [poly[e] for e in expsort]
        for (n3, _) in n_factors( max(abs(cfs[0]), abs(cfs[1])) - min(abs(cfs[0]), abs(cfs[1])) ):
            if 0 == cfs[0] % n3 == cfs[1] % n3:
                n1 = int(cfs[0]/n3)
                n2 = int(cfs[1]/n3)
                if cfs[3] % n2 == 0:
                    n4 = int(cfs[3]/n2)
                    out = [({ min(expsort[1],expsort[2]):n1, 0:n2 }, 1), ({ max(expsort[1],expsort[2]):n3, 0:n4 }, 1)]
                    break
    return out

def bf_int_coordinates(exps):
    out_cord = set()
    for i in range(2,101):
        k = 1/i
        if check_fact(exps,k):
            yield k
        if check_fact(exps,-k):
            yield -k
    for i in range(1,1001):
        if check_fact(exps,i):
            yield i
        if check_fact(exps,-i):
            yield -i

def check_fact(exps,fact):
    out = 0
    for exp in exps:
        out += exps[exp] * (fact**exp)
    return round(out,15) == 0

def factorize(poly_stack, func):
    poly = poly_stack.pop()
    tmexp = max(poly)
    div_polys = []
    common_factor = 1
    checknegative = set([c < 0 for c in poly.values()])
    # factorizing checks
    for (i, _) in n_factors(min([abs(v) for v in poly.values() if v != 0])):    # if common factor in poly, divide e.g. 2x^2+4 -> 2(x^2+2)
        checkmult = set()                                 # check performed on every iteration because of coeffs changing with division
        for coeff in poly.values():
            checkmult.add(coeff % i)
        if len(checkmult) == 1 and 0 in checkmult:
            common_factor = i if checknegative != set([True]) else -i
            break
    if common_factor != 1:
        div_polys = [ ({ 0:common_factor }, 1) ]
    elif len(poly) > 2 and tmexp and poly[0] == 0:    # x^5 + x^3 -> x^3(x^2 + 1)
        div_polys = [ ({ 1:1, 0:0 }, min([e if e > 0 else tmexp for e in poly])) ]
    elif len(poly) == 2 and max(poly) > 1 and poly[0]:    # x^2 - 1 -> (x + 1)(x - 1), x^3 - 1, x^3 + 1, etc.
        div_polys = pow_diff(poly)
    elif len(poly) == 3 and poly[0]:    # x^2 + 2x + 1 -> (x + 1)^2, 3x^2 + 7x + 2 -> (3x + 1)(x + 2), etc. max exp can be > 2
        expsort = sorted(poly)[::-1]
        if expsort[0] % 2 == 0 and expsort[0]-expsort[1] == expsort[1]-expsort[2]:
            div_polys = binomial_mult_3(poly, expsort)
    elif len(poly) == 4 and poly[0]:
        expsort = sorted(poly)[::-1]
        if expsort[0]-expsort[1] == expsort[1]-expsort[2] == expsort[2]-expsort[3]:
            div_polys = binomial_pow3(poly, expsort)
        if not div_polys:    # 6x^6 + 4x^4 + 15x^2 + 10 would trigger the first check but not the second when using ELIF (one doesn't exlude the other)
            if expsort[0]-expsort[2] == expsort[1]-expsort[3]:
                div_polys = binomial_mult_4(poly,expsort)
    if not div_polys and tmexp > 2:
        # bruteforce
        div_count = tmexp
        for xv in bf_int_coordinates(poly):
            div_polys.append(({ 1:1, 0:-xv }, 1))
            div_count -= 1
            if div_count == 0:
                break
    for p, e in div_polys:
        for div_i in range(e):
            poly = divide_func(poly, p)
            if (max(p) > 2) or (max(p) == 2 and p[0] and delta_calc(*x2terms(p)) >= 0):
                poly_stack.append(p)
            else:
                func.add(p, 1)
    if div_polys and ((max(poly) > 2) or (max(poly) == 2 and poly[0] and delta_calc(*x2terms(poly)) >= 0)):
        poly_stack.append(poly)
    else:
        if len(poly) == 2 and not poly[0]:    # fix for ax^2 -> x^2 divided by a -> poly = {2:1,0:0}:1, should be {1:1,0:0}:2
            func.add({ 1:1, 0:0 }, max(poly))
        else:
            func.add(poly, 1)
    if poly_stack:
        factorize(poly_stack, func)

def polyformat(polys, x0t):
    out = ['','']
    brackets = False
    if len(polys) > 1 or x0t != 1:
        brackets = True
    out[0] += sf2i(x0t) if x0t not in (1,-1) or len(polys) == 0 else '-' if x0t == -1 else ''
    for poly, exp in polys.items():
        poly = dict(poly)
        if len(poly) == 2 and not poly[0]:
            out[1] = 'x'
            if exp > 1:
                out[1] += '^' + str(exp)
        else:
            current_poly = ''
            if exp > 1:
                brackets = True
            expsort = sorted(poly)[::-1]
            for e in expsort:
                current_poly += '- ' if poly[e] < 0 else '+ ' if poly[e] > 0 else ''
                if e != 0:
                    current_poly += sf2i(abs(poly[e])) if poly[e] not in (1,-1) else ''
                    current_poly += 'x'
                    current_poly += '^' + sf2i(e) + ' ' if e != 1 else ' '
                else:
                    current_poly += sf2i(abs(poly[e])) if poly[e] else ''
            if current_poly[0] == '+':
                current_poly = current_poly[2:]
            elif current_poly[0] == '-' and brackets:
                current_poly = '-' + current_poly[2:]
            current_poly = '(' + current_poly + ')' if brackets else current_poly
            current_poly += '^' + sf2i(exp) if exp != 1 else ''
            out.append(current_poly)
    return ''.join(out)

class Function():
    def __init__(self, data):
        self.data = {}
        self.x0t = 1
        if type(data) == dict:
            self.exps = data
        else:
            self.eqt = remove_spaces(data)
            self.eqt = fix_signs(self.eqt)
            self.exps = convert(self.eqt)    # self.eqt is referenced and edited directly by convert()
            if 0 not in self.exps:    # 0 may already be in exps because of x^0 terms
                self.exps[0] = 0
            self.exps[0] += solve_x0_terms(self.eqt)    # x-terms have already been removed from self.eqt
        self.out = ""
    def __repr__(self):
        return repr(self.data)
    def add(self, indict, exp):
        if len(dict(indict)) == 1:
            self.x0t *= ((dict(indict))[0] ** exp)    # number-only terms (x^0) are managed separately
        else:
            self.indict = fix_dict(indict)
            if self.indict in self.data:
                self.data[self.indict] += exp
            else:
                self.data[self.indict] = exp
    def factorize(self):
        if set(self.exps.values()) != set([0]):
            factorize([self.exps], self)
            #print(self.data)
            return polyformat(self.data, self.x0t)
        else:
            return '0'
