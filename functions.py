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
