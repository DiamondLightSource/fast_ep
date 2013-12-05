
def pdb_file_nres(file_name):
    from iotbx import pdb
    pdb_obj = pdb.hierarchy.input(file_name=file_name)
    nres = sum([len(chain.as_sequence())
                for chain in pdb_obj.hierarchy.chains()])
    return nres

if __name__ == '__main__':
    import sys
    print pdb_file_nres(sys.argv[1])
                        
