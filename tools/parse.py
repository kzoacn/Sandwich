with open('result.txt') as f:
    lines = f.readlines()


paras=[]
for i,l in enumerate(lines):
    if l.startswith('FAEST'):
        para={}
        para['name']=l.split()[0]
        para['sig size']=int(l.split()[1])
        para['sk']=int(l.split()[2])
        para['pk']=int(l.split()[3])
        para['keygen']=lines[i+1].strip()
        para['sign']=lines[i+2].strip()
        para['verify']=lines[i+3].strip()
        paras.append(para)

#output markdown table
print('| name | sig size | sk | pk | keygen | sign | verify |')
print('| ---- | -------- | -- | -- | ------ | ---- | ------ |')
for p in paras:
    print('|',p['name'],'|',p['sig size'],'|',p['sk'],'|',p['pk'],'|',p['keygen'],'|',p['sign'],'|',p['verify'],'|')


        

