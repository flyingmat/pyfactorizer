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

def get_coefficient(inlst, i):
    coeff = ''
    k = i - 1
    while k >= 0 and inlst[k] in '1234567890.':    # keep going backwards to get the full coefficient
        coeff = inlst[k] + coeff
        k -= 1
    coeff = 1 if not coeff else coeff    # if no coefficient is specified, 1 is assigned
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
            self.eqt = data
        else:
            self.eqt = remove_spaces(data)
            self.eqt = fix_signs(self.eqt)
            self.exps = convert(self.eqt)    # self.eqt is referenced and edited directly by convert()
            self.exps[0] = solve_x0_terms(self.eqt)    # x-terms have already been removed from self.eqt
        self.out = ""
    def __repr__(self):
        return repr(self.data)
