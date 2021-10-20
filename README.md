# **Circular-DNA Design**

## Introduction

Our Circular DNA Design Software is made from the Python coding language. 
This software is built for GotCha's mechanism based on Rolling Circle Amplification (RCA). 
With user input of target micro-RNA sequence for capture, Circular DNA Design is able to synthesise different permutations of bases to form circular DNA sequences suitable for GotCha with the help of rules we established in circular DNA design. 
With various rules established, we successfully avoid formation of major loops in circular DNA sequences while maximising thermodynamic stability. 


## Code


### Input

In the input, we allow the user to input the name of the Circular DNA, which will then be the name of the txt. file exported after running the code. <br/>
Then, we show the default Immobilisation Probe binding site on the circular DNA we designed. <br/>
After which, according to the needs of the user, they can input the total length of the Circular DNA Sequence of their choice. 


```

name = input('Name of Circular-DNA = ')

probe_sequence = 'CAACCACACTGGCAAGAGGCAAAAAAAAAAAAAAA' 

probe_binding_site = 'GTGGTTGTCTTCT' 
mirna1 = input('Micro-RNA sequence (IN CAPITAL LETTERS) = ')
al = input('Length of all DNA sequence (>'+str(13+len(mirna1))+') (suggest=75) = ')

```

### Output

For the output of the code, we first show the DNA sequence of our immobilisation probe. <br/>
Then, we create a list, of which the different permutations of bases are added into the list and exported as a txt. file. <br/>


```

print('\n ===================================================================================')
print('\n                                  Probe sequence')
print('\n                       CAACCACACTGGCAAGAGGCAAAAAAAAAAAAAAA')
AUGC_table = {'A':'T','U':'A','G':'C','C':'G'}
ATGC_table = {'A':'T','T':'A','G':'C','C':'G'}
base = 'TCGA'
codons = ''
mirna = []
for a in mirna1:
    mirna.append(''.join(AUGC_table[a]))
fn_name = name + '_circular_dna.txt'
txtfn = open(fn_name,'w')
ff = []
q = ''
print('\n ===================================================================================')
print('\n                           Circular DNA sequence('+str(int(al))+')')
print('\n       GTGGTTGTCTTCT      +  '+q.join(mirna)+'   +  _________________')
print('   Probe_binding_site(13) +  MicroRNA_binding_site('+str(len(mirna))+') +  Unctionless_dna_sequence('+str(int(al)-13-int(len(mirna)))+')')
print('\n ===================================================================================')

```

### Parameters
For the parameters used in our code, we first create lists for probe sequence complementary, probe binding site, circular probe, and functionless sequence bases. <br/>

```

probe_sequence_complementary_list = []
probe_binding_site_mirna_complementary_list = []
circular_list = []
nonss = []


probe_count = -1
probe_counts = ''
probe_inverted_count = len(probe_sequence)
probe_inverted_counts = ''
circular_count = 0
circular_counts = ''
probe_binding_site_mirna_inverted_count = len(probe_binding_site)+len(mirna)
probe_binding_site_mirna_inverted_counts = ''
alll = int(al)



probe_complementary = ''
probe_inverted_complementary = ''
probe_binding_site_mirna_inverted_complementary = ''
nons = ''

```



### Probe Design
For the Circular DNA Probe Design, we considered: <br/>
complementary and reverse complementary sequence for the immobilisation probe and <br/>
complementary and reverse complementary sequence for immobilisation probe binding site.  <br/>

```

for a in probe_sequence:
    probe_sequence_complementary_list += ATGC_table[a]
    probe_count += 1
    probe_inverted_count -= 1
    probe_counts += '{} '.format(str(probe_count))
    probe_inverted_counts += '{} '.format(str(probe_inverted_count))
    
probe_complementary_dict = dict(zip(probe_counts.split(' '),probe_sequence_complementary_list))
probe_complementary_inverted_dict = dict(zip(probe_inverted_counts.split(' '),probe_sequence_complementary_list))
    
                          
for a in probe_complementary_dict.keys():
    if int(a)+2 <= probe_count:
        x = probe_complementary_dict[a]
        y = probe_complementary_dict[str(int(a)+1)]
        z = probe_complementary_dict[str(int(a)+2)]
        probe_complementary += '{}{}{}\t'.format(x,y,z)
  
nons += probe_complementary


for a in probe_complementary_inverted_dict.keys():
    if int(a)+2 <= probe_count:
        x = probe_complementary_inverted_dict[a]
        y = probe_complementary_inverted_dict[str(int(a)+1)]
        z = probe_complementary_inverted_dict[str(int(a)+2)]
        probe_inverted_complementary += '{}{}{}\t'.format(x,y,z)

nons += probe_inverted_complementary


for a in probe_binding_site:
    circular_list += a
    probe_binding_site_mirna_complementary_list += ATGC_table[a]
    circular_count += 1
    probe_binding_site_mirna_inverted_count -= 1
    circular_counts += '{} '.format(str(circular_count))
    probe_binding_site_mirna_inverted_counts += '{} '.format(str(probe_binding_site_mirna_inverted_count))

```

### Reverse complementary sequence of miR
Then, we also took into account the reverse complementary sequence of the target miRNA itself. 

```

for a in mirna:
    circular_list += a
    probe_binding_site_mirna_complementary_list += ATGC_table[a]
    circular_count += 1
    probe_binding_site_mirna_inverted_count -= 1
    circular_counts += '{} '.format(str(circular_count))
    probe_binding_site_mirna_inverted_counts += '{} '.format(str(probe_binding_site_mirna_inverted_count))

codons_dict = dict(zip(circular_counts.split(' '),circular_list))       
probe_binding_site_mirna_complementary_inverted_dict = dict(zip(probe_binding_site_mirna_inverted_counts.split(' '),probe_binding_site_mirna_complementary_list))

for a in probe_binding_site_mirna_complementary_inverted_dict.keys():
    if int(a)+2 <= len(probe_binding_site)+len(mirna)-1:
        x = probe_binding_site_mirna_complementary_inverted_dict[a]
        y = probe_binding_site_mirna_complementary_inverted_dict[str(int(a)+1)]
        z = probe_binding_site_mirna_complementary_inverted_dict[str(int(a)+2)]
        probe_binding_site_mirna_inverted_complementary += '{}{}{}\t'.format(x,y,z)
            
nons += probe_binding_site_mirna_inverted_complementary      
AT = 0
w = 0

```

### Reverse Complementary Sequence of Functionless DNA Sequence 
Lastly, we made sure not to repeat reverse complementary sequences of the functionless DNA sequence in the Circular DNA probe. 

```

nonss = nons.split('\t')
circular_countss = []
circular_countss = circular_counts.split(' ')
circular_countss.pop()
def f(n,fn_name,AT,ff,base,nons,circular_counts,codons_dict,nonss,circular_list,ATGC_table,circular_count,circular_countss,alll):
    for codon in base:
        non = 0
        if circular_count >= n:
            circular_count = n-1
            circular_countss = circular_countss[:n-1]
            circular_list = circular_list[:n-1]
            codons_dict = dict(zip(circular_countss,circular_list))
            nonss = nonss[n+1:]     
        if codon == 'A': 
            x = codon
            y = codons_dict[str(circular_count)]
            z = codons_dict[str(circular_count-1)]
            codons = '{}{}{}'.format(z,y,x)
            if codons == 'AAA':
                non += 1
            for b in nonss:
                if codons == b:
                    non += 1
        if codon == 'T': 
            x = codon
            y = codons_dict[str(circular_count)]
            z = codons_dict[str(circular_count-1)]
            codons = '{}{}{}'.format(z,y,x)
            if codons == 'TTT':
                non += 1
            for b in nonss:
                if codons == b:
                    non += 1
        if codon == 'G': 
            x = codon
            y = codons_dict[str(circular_count)]
            z = codons_dict[str(circular_count-1)]
            codons = '{}{}{}'.format(z,y,x)
            if codons == 'GGG':
                non += 1
            for b in nonss:
                if codons == b:
                    non += 1
        if codon == 'C': 
            x = codon
            y = codons_dict[str(circular_count)]
            z = codons_dict[str(circular_count-1)]
            codons = '{}{}{}'.format(z,y,x)
            if codons == 'CCC':
                non += 1
            for b in nonss:
                if codons == b:
                    non += 1             
        if non == 0:
            add_codons = '{}{}{}\t'.format(ATGC_table[x],ATGC_table[y],ATGC_table[z])
            nons += add_codons
            nonss = nons.split('\t') 
            circular_count += 1
            circular_list += codon 
            circular_counts = '{} '.format(str(circular_count)) 
            circular_countss += circular_counts.split(' ')
            circular_countss.pop()
            codons_dict = dict(zip(circular_countss,circular_list))
            if circular_count == alll:
                AT = 0
                for a in circular_list:
                    if a == 'A' or a == 'T':
                        AT += 1
                if AT >= alll/2-2 and AT <= alll/2+2:
                    ff.append(''.join(circular_list))
            else:
                f(n+1,fn_name,AT,ff,base,nons,circular_counts,codons_dict,nonss,circular_list,ATGC_table,circular_count,circular_countss,alll) 
    return (' ')
        
        
print(f(circular_count+1,fn_name,AT,ff,base,nons,circular_counts,codons_dict,nonss,circular_list,ATGC_table,circular_count,circular_countss,alll))
        
if len(ff) == 0:
    print('\nERROR\n\n[Maybe you should invert the micro-RNA and try it again.]')
else:
    print('\nSUCCESS')    
    
print('\nTotal '+str(len(ff))+' sequence ')

print('\nOutput in '+fn_name)


txtfn.write('{}\n'.format(ff))
txtfn.close()
        
```      
