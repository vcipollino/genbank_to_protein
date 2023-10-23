#### import libraries
import sys

#testing new branch 

####USER DEFINED VARAIBLES
#input_file = 'corona_virus_genbank'
input_file = sys.argv[1]
output_DNA = input_file + '.DNA'
output_RNA = input_file + '.RNA'
output_protein = input_file + '.protein'
output_correct_translation = input_file + '.correct_translation'


####USER DEFINED FUNCTIONS##
#Parse GENBANK format file into a dictionary where geneid=key and coordinates=value
def parse_genbank(genbank_data):
    parsing = True
    gene_next_line = False
    gene_coordinate_dictionary ={}
    for line in genbank_data:
        line=line.strip() 
        if line.startswith("FEATURES"):
            parsing = True
        elif line.startswith("ORIGIN"):
            parsing = False
        if parsing:
            line = line.strip() #remove whitespace in line
            #Get the gene ID if the line begins with 'gene'
            if line.startswith("gene"):
               line_elements = line.split()
               both_coordinates = line_elements[1].split('..')
               start = both_coordinates[0]
               stop = both_coordinates[1]
               clean_coordinates = start + "," + stop #merge coordinates delimited with a comma
               gene_next_line = True
            #if CDS as different coordinates extract and use those coordinates
            elif line.startswith("CDS"):
                line_elements_2 = line.split() #split line into CDS title and join(..coordinates..)
                slippage= False 
                if line_elements[1] != line_elements_2[1]:
                    slippage = True   # if line_elements of the CDS and gene are the different run slippage
                if slippage:
                    i_1 = line.find("(") #extract coordinates delimited by a comma in the parenthesis 
                    i_2 = line.find(")")
                    correct_coordinate = line[i_1 +1:i_2]
                    correct_coordinate = correct_coordinate.split(',')
                    position_1 = correct_coordinate[0].split('..')
                    position_2 = correct_coordinate[1].split('..')
                    #store start and stop of each coordinate
                    start1 = position_1[0] 
                    stop1 = position_1[1] 
                    start2 = position_2[0]
                    stop2 = position_2[1]
            #join the  start to stop of coordinate 1 and start to stop of coordinate 2
                    clean_coordinates2 = start1 + "," + stop1 + ";" + start2 + "," + stop2
                    gene_coordinate_dictionary[gene_name] = clean_coordinates2

            #Get the coordiantes after the line  ID if the line begins with 'gene'    
            elif gene_next_line:
                gene_name = line
                gene_name = gene_name.strip()
                gene_name = gene_name.replace("\"","")
                gene_name = gene_name.replace("/gene=","")                              
                gene_next_line = False #stops after line with coordinates
                #Load the gene_id-->coordinates(comma-delimited) into a dictionary
                gene_coordinate_dictionary[gene_name] = clean_coordinates
    return(gene_coordinate_dictionary)
    
    
#MAIN
#input_file = 'corona_virus_genbank'
with open(input_file, 'r') as input_file:
    corona = [line.strip() for line in input_file]#copy input file to a list called corona
    
#Close input files
input_file.close()

#Parse GENBANK file into a dictionary where geneid=key and coordinates=value
gene_coordinates = parse_genbank(corona)


##### PHASE TWO: extract coding DNA from genome to make DNA sequence FASTA format file #########

#create and open corona_DNA_seq file
DNA_file = open(output_DNA, 'w')

# create new list with only DNA sequence called corona_DNA
for index, line in enumerate(corona): #assign index to each line
    if line.startswith('ORIGIN'):
        start = index #assign line number of line starting with 'ORIGIN' as start index
    if line.startswith('//'):
        end = index + 1 #assign line number of line starting with '//' as end index {+1 because last line is not included in range}
corona_DNA = corona[start:end] 
corona_DNA = ''.join(corona_DNA) #make corona_DNA one long list
corona_DNA = corona_DNA[6:-2] #remove 'ORIGIN' and '//'  from list

corona_DNA = corona_DNA.replace(" ","") #remove spaces from list

corona_DNA = ''.join(filter(lambda x: not x.isdigit(), corona_DNA)) #filter through data by removing all digit values

#loop through corona_DNA and print gene ids followed by values with index in range of start and stop coordinates

for value in gene_coordinates:
    coordinate = gene_coordinates[value]#assign coordinates only to varaiable
    coordinate = coordinate.split(";")
    sequences= []
    sequence2 = ''
    if len(coordinate) >= 2: 
        for coord in coordinate:
            coord= coord.split(",")
            start = int(coord[0]) #assign start and stop coordinates
            stop = int(coord[1])
            sequence = corona_DNA[start-1:stop] 
            sequence2 = sequence2 + sequence #combine the first sequences into one sequence
            #print(sequence2) #currently the last value in sequence2 is the correct one for output
            #sequences.append(sequence2)
    else: 
        for coord in coordinate:
            coord= coord.split(",")
            start = int(coord[0]) #assign start and stop coordinates
            stop = int(coord[1])
            sequence = corona_DNA[start-1:stop] 
            sequences.append(sequence)
    sequences.append(sequence2)
    DNA_file.write(f">{value}\n")
    for sequence in sequences:
        DNA_file.write(f"{sequence}")
    DNA_file.write(f"\n")

#### PHASE THREE: transcribe each DNA sequence for each gene and write to RNA sequence file########

#load DNA sequence file into a list called DNA                                                                        ^
DNA_file = output_DNA
#open DNA sequence file and load file into DNA
with open(DNA_file, 'r') as input_file: 
    DNA = [line.strip() for line in input_file]
#close file
input_file.close()

#for each line in DNA convert Ts to Us and save into a new list called RNA
#new list called RNA, which is a copy of DNA
RNA = list(DNA)

for i in range(len(RNA)):
    if not RNA[i].startswith('>'):
        RNA[i] = RNA[i].upper()
# RNA = [line.upper()for line in RNA]#if list is lower case, make upper case

for base, sequence in enumerate(RNA): #loop through each sequence 
    RNA[base] = sequence.replace('T', 'U')#replace Ts with Us


#open output file 
RNA_file = open(output_RNA, 'w')

#write RNA to output file in FASTA format
#print first line and replace 'dna' to 'rna'
# label_change = {'cdna' : 'rna', 'DNA' : 'rna'}
for line in RNA:
    line = line.replace('cdna','crna')
    line = line.replace('DNA', 'crna')
    RNA_file.write(f"{line}\n")

#close output (RNA) file
RNA_file.close()


#### PHASE FOUR: translate each RNA sequence for each gene to protein ######

#load RNA sequence file into a list called RNA 
RNA_file = output_RNA
#open DNA sequence file and load file into DNA
with open(RNA_file, 'r') as input_file: 
    RNA = [line.strip() for line in input_file]
#close file
input_file.close()

#download genetic code dictionary 
rna2protein = {'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L',
'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',
'UAU':'Y', 'UAC':'Y', 'UAA':'', 'UAG':'',
'UGU':'C', 'UGC':'C', 'UGA':'', 'UGG':'W',
'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
'AUU':'I', 'AUC':'I', 'AUA':'I', 'AUG':'M',
'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',
'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}


#open output file 
protein_file = open(output_protein,'w')

#match codon in RNA_string to codon in rna2protein dictionary
for line in RNA: 
    if line.startswith('>'):
        line = line.replace('rna','protein')
        protein_file.write(f"{line}\n") #write FASTA header to output file
    if line.startswith(('A','G','C','U')):
        for i in range(0,len(line),3):#start after >RNA and read every 3 letters
            codon = line[i:i+3]#assign every 3 nucleotides as a codon
            if len(codon)==3: #only match amino_acid to three letter codon
                if codon in rna2protein:
                    amino_acid = rna2protein[codon]#search for codon in rna2protein dictionary
                protein_file.write(f"{amino_acid}")#write amino acids to output file
        protein_file.write(f"\n")
            
#close protein_sequence file
protein_file.close()
#Purpose: parse genbank file to extrapolate correct translations

############################ MAIN ############################

# open writeable output file for correct translation 
output_file = open(output_correct_translation, 'w')


#create empty dictionary
genbank_correct_translation = {}

############################ save gene names to variable gene_identifier ############################
line_after = False
parse = False
for line in corona: # loop through each line in genbank file using list called corona
    if line.startswith("gene"): # if line starts with gene then go to next line
        line_after = True
        parse = False
        correct_translation = '' #reset correct_translation 
    elif line_after: #at next line strip all characters beside gene identifier and save to variable called gene_identifier
        gene_identifier = line
        gene_identifier = gene_identifier.strip()
        gene_identifier = gene_identifier.replace("\"","")
        gene_identifier = gene_identifier.replace("/gene=","") 
        line_after = False
        
############################ save to translation of each gene to variable correct_translated ############################

    if line.startswith("/translation"):#if line starts '/translation' then parse is True
        parse = True
    # elif '"' in line:
    #     parse = True
    elif line.startswith("3'UTR"):
        parse = False
    if parse: #parse out translation by stripping '/translation=' and remove quotation
        translation = line.replace("\"","")
        translation = translation.replace("/translation=", "")
        translation = translation.upper()
        correct_translation = correct_translation + translation #for each translation each line is added to correct_translation
        genbank_correct_translation[gene_identifier]= correct_translation #write to dictionary with gene identifer as the key and the translation as the value

############################ write gene name and correct translations in FASTA format ############################

#for key in genbank_correct_translation dictionary print '>'gene_indentifier followed by translation to outputfile follo
for value in genbank_correct_translation:
    translate = genbank_correct_translation[value]
    output_file.write(f">{value}\n{translate}\n")


