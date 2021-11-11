import Bio

def check_key(key,dict):
    if key in dict.keys():
        return True
    else:
        return False

def get_id(feature):
    if check_key('locus_tag',feature.qualifiers):
        return feature.qualifiers['locus_tag'][0]
    elif check_key('db_xref',feature.qualifiers):
        return feature.qualifiers['db_xref'][0]
    elif check_key('protein_id',feature.qualifiers):
        return feature.qualifiers['protein_id'][0]
    elif check_key('old_locus_tag',feature.qualifiers):
        return feature.qualifiers['old_locus_tag'][0]
    else:
        return "N/A"

def get_name(feature):
    if check_key('gene',feature.qualifiers):
        return feature.qualifiers['gene'][0]
    elif check_key('gene_synonym',feature.qualifiers):
        return feature.qualifiers['gene_synonym'][0]
    elif check_key('product',feature.qualifiers):
        return feature.qualifiers['product'][0]
    else:
        return "N/A"

def get_description(feature):
    if check_key('product',feature.qualifiers):
        return feature.qualifiers['product'][0]
    elif check_key('note',feature.qualifiers):
        return feature.qualifiers['note'][0]
    # elif check_key('gene',feature.qualifiers):
    #     return feature.qualifiers['gene'][0]
    else:
        return "N/A" 

def get_strand(feature):
    if feature.location.strand == 1:
        return "+"
    elif feature.location.strand == -1:
        return "-"
    else:
        return