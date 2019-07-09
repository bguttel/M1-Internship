import numpy as np
import os
from os import listdir
from os.path import isfile, join
import decimal

class jobfile:
    def __init__(self):
        self.preview = []
        self.content = []
        self._cl = '#'+'-'*79 #commentline
        self._el = '' #emptyline

    def addel(self):
        c = self.content
        c.append(self._el)

    def addcl(self):
        c = self.content
        c.append(self._cl)

    def addline(self,line):
        c = self.content
        c.append(line)

    def generatepreview(self):
        p = []
        c = self.content
        for li in c:
            li = self.linetostring(li)
            p.append(li+'\n')
        self.preview = p
        return p

    def __repr__(self):
        p = self.preview
        p_flat = ''
        for li in p:
            p_flat +=li
        return p_flat

    def getdecimalsamount(self,number):
        if isinstance(number,int):
            n = 0
        else:
            number = str(number)
            num = decimal.Decimal(number)
            n = num.as_tuple().exponent
            n = abs(n)
        return n

    def correctnumberformat(self,number):
        n_decimal = self.getdecimalsamount(number)
        # n_decimal = 1
        format = ',.{}f'.format(n_decimal)
        n_corrected = str('{:{}}'.format(number,format))
        # n_corrected = str('{:,.1f}'.format(number,format))
        n_corrected = n_corrected.replace(',',' ')
        return n_corrected

    def correctalignment(self,line,n_space = 25):
        l = line.split('\t')
        line_output = ''
        for li in l:
            line_output += '{:<{}}'.format(li,n_space)
        return line_output

    def linetostring(self,line):
        line_str = ''
        if type(line) != list:
            line = [line]

        for param in line:
            num_bool = isinstance(param,(tuple,np.ndarray,list))
            if num_bool:
                p_corrected = '('
                for ni in param:
                    ni = self.correctnumberformat(ni)
                    p_corrected += ni + ', '
                p_corrected = p_corrected[0:-2] + ')' #remove last ', '
            else:
                p_corrected = str(param)
            line_str += p_corrected + '\t'

        line_str = self.correctalignment(line_str)
        return line_str

    def generatefile(self,filename,path):
        extension = '.njf'
        if extension not in filename:
            filename += extension

        cwd_files = get_path_files(path)
        filename = find_unused_name(filename,cwd_files,extension)

        full_path = join(path,filename)

        if not self.preview:
            self.generatepreview()

        f = open(full_path,'w')
        for li in self.preview:
            f.write(li)
        f.close()

    # -- Higher level commands
    def addheader(self):
        c = self.content
        line = '# run nbwrite path/YYYYMMDD/filename.njf -1=shintaro:g1  -2=shintaro:g2'
        c.append(line)

    def addfooter(self):
        c = self.content
        c.append('.write')
        c.append('current\tauto')
        c.append('.write')

def get_path_files(path):
    files = [f for f in listdir(path) if isfile(join(path, f))]
    return files

def find_unused_name(name, list_names, extension='.npj'):
    """
    name with a given extension
    """
    l = len(extension)
    if extension not in name:
        name += extension

    only_name = name[:-l]

    i = 0
    new_name = only_name + '_' + str(i) + extension
    while new_name in list_names:
        new_name = only_name + '_' + str(i) + extension
        i+=1

    return new_name

if __name__ == '__main__':
    path = '/Users/Wanxii/Desktop/Test scripts'
    filename = 'test_jobfile'
    f = jobfile()
    f.addheader()
    f.addcl()
    f.addline(['dose',[+20,1000000],'hello'])
    f.addel()
    f.addcl()
    f.addfooter()
    f.generatepreview()
    # f.generatefile(filename,path)
    print(f)
