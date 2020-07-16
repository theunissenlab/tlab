from numpy import corrcoef

def similarity_index(strf1,strf2):
    s1 = strf1.data.flatten()
    s2 = strf2.data.flatten()
    return corrcoef(s1,s2)[0,1]