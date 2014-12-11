import numpy as np
import pandas
import os

def write_plink(G, F, only_snps_matching=None, intermediate_file_name='dataset', final_file_name='dataset_recode',make_bed=False, maf = 0.05):
    if only_snps_matching is not None:
        G_filtered = G.filter(regex=only_snps_matching)
    else:
        G_filtered = G

    print G_filtered.shape[1], "SNPs selected"
    fn_map = intermediate_file_name+'.map'
    fn_ped = intermediate_file_name+'.ped'
    write_MAP(G_filtered, filename=fn_map)
    write_PED(G_filtered, F, filename=fn_ped)
    if make_bed:
        import external.plink as pn
        plinkstr = '%s --noweb --file %s --compound-genotypes --make-bed --out %s' % (pn.executable,intermediate_file_name, final_file_name)
        if maf>0.0:
            plinkstr = '%s --maf %.4f'%(plinkstr,maf)
        os.system(plinkstr)
        os.remove(fn_map)
        os.remove(fn_ped)
    
def write_MAP(G, filename='dataset.map', default_chr=1):
    snp_names = np.array(G.columns[:,None], dtype=str)
    chromosomes = np.ones((snp_names.shape[0],1), dtype=int)
    genetic_distances = np.arange(snp_names.shape[0])[:, None]
    bp_positions = np.arange(snp_names.shape[0])[:, None]
    map_data = np.concatenate((chromosomes, snp_names, genetic_distances, bp_positions), axis=1)
    np.savetxt(filename, map_data, delimiter='\t', fmt='%s')
    print 'MAP file written'

def write_PED(G, F, filename='dataset.ped'):
    pheno = np.array(np.ones((F.shape[0], 1)) * -9, dtype=int)
    ped_data = np.concatenate((F['Family'].values[:, None], F.index[:, None],
                               F['Paternal ID'].values[:, None], F['Maternal ID'].values[:, None],
                               F['Sex'].values[:, None], np.array(pheno, dtype=str)), axis=1)

    genotype = np.empty(G.shape, dtype='|S2')
    # From now on let's assume that the minor allele is always G and the major A
    genotype[G.values==2] = 'GG'
    genotype[G.values==1] = 'AG' # AG = GA
    genotype[G.values==0] = 'AA'
    ped_data = np.concatenate((ped_data, genotype), axis=1)
    np.savetxt(filename, ped_data, delimiter='\t', fmt='%s')
    print 'PED file written'

if __name__ == '__main__':
    geno = np.random.binomial(2, 0.5, size=(100, 1000))
    sample_names = ['sample_%d' % i for i in range(100)]
    snp_names = ['snp_%d' % i for i in range(1000)]
    family_names = np.array(['family_%d' % i for i in range(100)], dtype=str)[:, None]
    paternal_ids = np.zeros_like(family_names)
    maternal_ids = np.zeros_like(family_names)
    sex = np.zeros_like(family_names)
    fam_pid_mid_sex = np.concatenate((family_names, paternal_ids, maternal_ids, sex), axis=1)
    G = pandas.DataFrame(data=geno, index=sample_names, columns=snp_names)
    F = pandas.DataFrame(data=fam_pid_mid_sex, index=sample_names, columns=['Family', 'Paternal ID', 'Maternal ID', 'Sex'])

    write_plink(G, F)
