import os

_ROOT = os.path.abspath(os.path.dirname(__file__))

mlst_urls = {'Acinetobacter baumannii Oxf': [
    'http://pubmlst.org/data/profiles/abaumannii.txt',
    'http://pubmlst.org/data/alleles/abaumannii/Oxf_gltA.tfa',
    'http://pubmlst.org/data/alleles/abaumannii/Oxf_gyrB.tfa',
    'http://pubmlst.org/data/alleles/abaumannii/Oxf_gdhB.tfa',
    'http://pubmlst.org/data/alleles/abaumannii/Oxf_recA.tfa',
    'http://pubmlst.org/data/alleles/abaumannii/Oxf_cpn60.tfa',
    'http://pubmlst.org/data/alleles/abaumannii/Oxf_gpi.tfa',
    'http://pubmlst.org/data/alleles/abaumannii/Oxf_rpoD.tfa',
    ],
    'Acinetobacter baumannii Pas': [
        'http://pubmlst.org/data/profiles/abaumannii_2.txt',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_cpn60.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_fusA.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_gltA.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_pyrG.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_recA.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_rplB.tfa',
        'http://pubmlst.org/data/alleles/abaumannii_2/Pas_rpoB.tfa',
    ],
    'Enterococcus faecalis': [
        'http://pubmlst.org/data/profiles/efaecalis.txt',
        'http://pubmlst.org/data/alleles/efaecalis/gdh.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/gyd.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/pstS.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/gki.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/aroE.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/xpt.tfa',
        'http://pubmlst.org/data/alleles/efaecalis/yqiL.tfa',
    ],
    'Enterococcus faecium': [
        'http://pubmlst.org/data/profiles/efaecium.txt',
        'http://pubmlst.org/data/alleles/efaecium/atpA.tfa',
        'http://pubmlst.org/data/alleles/efaecium/ddl.tfa',
        'http://pubmlst.org/data/alleles/efaecium/gdh.tfa',
        'http://pubmlst.org/data/alleles/efaecium/purK.tfa',
        'http://pubmlst.org/data/alleles/efaecium/gyd.tfa',
        'http://pubmlst.org/data/alleles/efaecium/pstS.tfa',
        'http://pubmlst.org/data/alleles/efaecium/adk.tfa',
    ],
    'Staphylococcus aureus': [
        'http://pubmlst.org/data/profiles/saureus.txt',
        'http://pubmlst.org/data/alleles/saureus/arcc.tfa',
        'http://pubmlst.org/data/alleles/saureus/aroe.tfa',
        'http://pubmlst.org/data/alleles/saureus/glpf.tfa',
        'http://pubmlst.org/data/alleles/saureus/gmk_.tfa',
        'http://pubmlst.org/data/alleles/saureus/pta_.tfa',
        'http://pubmlst.org/data/alleles/saureus/tpi_.tfa',
        'http://pubmlst.org/data/alleles/saureus/yqil.tfa',
    ], }

