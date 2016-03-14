#import urllib.request
import sys
import urllib
# import path
import pkg_resources
import cPickle
import jellyfish
from Bio import SeqIO
import os
from straintypemer import _ROOT
from straintypemer import mlst_urls


def update_mlst_resources():
    resource_path = os.path.join(_ROOT, "data/mlst_resources/")
    if os.path.exists(resource_path) is False:
        os.mkdir(resource_path)
    file_lists = {}
    for k in mlst_urls.keys():
        strain_path = os.path.join(resource_path, k)
        if os.path.exists(strain_path) is False:
            os.mkdir(strain_path)

        sys.stderr.write("Updating '{0}' resources\n".format(k))
        sys.stderr.write("...saving_to '{0}' resources\n".format(resource_path))
        file_lists.update({k : []})
        for url in mlst_urls[k]:
            file_lists[k].append(url.split("/")[-1].lower().replace("_.","."))
            sys.stderr.write('\tretrieving: {0}\n'.format(url))
            sys.stderr.write('\tsaving: {0}\n'.format(url.split("/")[-1].lower().replace("_.",".")))
            urllib.urlretrieve(url, os.path.join(strain_path, url.split("/")[-1].lower().replace("_.",".") ))
    pickle_profiles(file_lists, resource_path)
    return


def pickle_profiles(file_lists, resource_path, kmer_size=31):
    jellyfish.MerDNA_k(kmer_size)
    #instantiate the pickle obj
    mlst_profiles_dict = {}
    for species, file_list in file_lists.iteritems():

        if species not in mlst_profiles_dict:
            mlst_profiles_dict.update({species : {"ST" : {}, "GENES" : {}, "GENE_ORDER" : None}})

        number_of_genes = len(file_list) - 1
        for i, l in enumerate(open([os.path.join(resource_path, species, f)
                                    for f in file_list if '.txt' in f][0])):
            line = l.strip().split("\t")
            if i == 0:
                gene_list = line[1:number_of_genes + 1]
            else:
                profile = ":".join(line[1:number_of_genes + 1])
                st = line[0]
                mlst_profiles_dict[species]["ST"].update({profile : st})
                mlst_profiles_dict[species]["GENE_ORDER"] = gene_list

        for _file in [f for f in file_list if f[-4:] == '.tfa']:
            for seq_record in SeqIO.parse(os.path.join(resource_path, species,_file), 'fasta'):
                seq_num = seq_record.name.split("_")[-1]
                gene_name = "_".join(seq_record.name.replace("__","_").replace("-","_").split("_")[:-1])
                if gene_name not in mlst_profiles_dict[species]["GENES"]:
                   mlst_profiles_dict[species]["GENES"].update({gene_name : {seq_num : set([])}})
                else:
                   mlst_profiles_dict[species]["GENES"][gene_name].update({seq_num : set([])})
                for j in range(0, len(seq_record.seq) - kmer_size + 1):
                    kmer = seq_record.seq[j : j + kmer_size]
                    mer = jellyfish.MerDNA(str(kmer))
                    mer.canonicalize()
                    mlst_profiles_dict[species]["GENES"][gene_name][seq_num].add(str(mer))
            sys.stderr.write("\tparsing: {0} : {1}\n".format(species, gene_name))
    cPickle.dump(mlst_profiles_dict, open(os.path.join(resource_path, "mlst_profiles.pkl"), "wb") )
    for f in file_lists.itervalues():
        os.remove(f)
    return

