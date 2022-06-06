#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :gem2adata.py
@说明        : trans gem to anndata
@时间        :2022/03/18 23:15:06
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''

import math
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse as sp
import logging

logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG)


def getArgs():
    """recceive the argument

    Returns:
        args: args
    """
    parser = argparse.ArgumentParser(
        description=__doc__, prog='gem2adata.py')
    parser.add_argument('-i', action='store', dest='in_gem', type=str, required=True,
                        help='your gem file')
    parser.add_argument('-b', action='store', dest='bin', type=int, default=50,
                        help='length of bin, default = 50')
    parser.add_argument('-l', action='store', dest='clu_list', type=str,
                        help='the list file of genes to cluster, default is all')
    return parser.parse_args()


def gemZoomer(gem_file, clu_list):
    """zoom the gem to max(full) size

    Args:
        gem_file (str): gem file path
        clu_list (str): cluster list

    Returns:
        gem: zommed gem
    """

    try:
        read_type_dict = {"gene": str, "x": int,
                          "y": int, "MIDcounts": int, "status": str}
        gem = pd.read_csv(gem_file, sep='\t', dtype=read_type_dict)
    except pd.errors.ParserError:
        logging.error(
            '\033[31m Reading gem file Error, Please check ! \033[0m')
    if clu_list != None:
        indices = gem['geneID'].isin(clu_list)
        gem = gem[indices]
    if 'status' not in gem.columns:
        gem['status'] = 'S'
    min_x = gem['x'].min()
    max_x = gem['x'].max()
    min_y = gem['y'].min()
    max_y = gem['y'].max()
    x_scal = max_x-min_x+1
    y_scal = max_y-min_y+1
    logging.info("\033[32m max x:{} max y: {} \033[0m".format(max_x, max_y))
    # pattern = re.compile(r'\d+,\d+~\d+,\d+')
    # if re.match(pattern, region_arg):
    #     info = re.split(',|~', region_arg)
    #     if int(info[0]) >= 0 and int(info[1]) <= x_scal \
    #             and int(info[0]) < int(info[1]) and int(info[2]) >= 0 \
    #             and int(info[3]) <= y_scal and int(info[2]) < int(info[3]):
    #         region = (int(info[0]), int(info[1]), int(
    #             info[2]), int(info[3]), 'PART')
    #     else:
    #         print("\033[1;31mInvalid region info, select complete region.\033[0m")
    #         region = (0, x_scal-1, 0, y_scal-1, 'ALL')
    # elif region_arg == 'ALL':
    # logging.info("\033[32m Split pos into x and y \033[0m")("Select complete region.")
    region = (0, x_scal-1, 0, y_scal-1, 'ALL')

    x_scal = region[1]-region[0]+1
    y_scal = region[3]-region[2]+1
    gem['x'] = gem['x']-min_x-region[0]
    gem['y'] = gem['y']-min_y-region[2]
    gem = gem[(gem.x >= 0) & (gem.x < x_scal) &
              (gem.y >= 0) & (gem.y < y_scal)]

    return gem


# def gemArranger(dir, skiprows):

#     print('step1. capture boundaries')
#     files = os.listdir(dir)
#     if len(files) == 0:
#         exit('\033[1;31mERROR: No gem files in "{}".\033[0m'.format(dir))
#     pattern = re.compile(r'.+\.gem')
#     max_x = 0
#     max_y = 0
#     gem_num = 0
#     for file in files:
#         if re.match(pattern, file):
#             gem_num += 1
#             file_path = dir.rstrip('/\\') + '/' + file
#             gem = pd.read_csv(file_path, delimiter='\t', skiprows=skiprows)
#             x_scal = gem['x'].max() - gem['x'].min()
#             y_scal = gem['y'].max() - gem['y'].min()
#             max_x = max(x_scal, max_x)
#             max_y = max(y_scal, max_y)
#     print('boundaries captured')
#     return max_x, max_y, gem_num


def transGemToAnnData(gem, bin, sample):
    """trans gem data to adata

    Args:
        gem (gem): gem
        bin (int): bin size
        sample (str): sample name

    Returns:
        adata: adata
    """

    logging.info("\033[32m Trans gem to AnnData format. \033[0m")
    # count info
    bin = int(bin)
    x_pixel = math.ceil((gem['x'].max()-gem['x'].min()+1)/bin)
    y_pixel = math.ceil((gem['y'].max()-gem['y'].min()+1)/bin)
    logging.info('\033[32m The compressed matrix shape is x={},y={} \033[0m'.format(
        x_pixel, y_pixel))

    # compress the matrix by dataframe
    gem['cx'] = np.floor(gem['x']/bin).astype(int)
    gem['cy'] = np.floor(gem['y']/bin).astype(int)
    gem['bin'] = gem['cx'].astype(str) + '_' + gem['cy'].astype(str)

    # extract the bin and gene categorical in dataframe
    bin_cat = pd.Categorical(gem['bin'])
    '''
    May be bugs here !!!
    '''
    # bin_list = np.sort(np.array(bin_cat.unique()))
    bin_list = np.sort(gem['bin'].unique())  # WWJ
    gem['bin_index'] = bin_cat.codes

    gene_cat = pd.Categorical(gem['gene'])
    # gene_list = np.sort(np.array(gene_cat.unique()))
    gene_list = np.sort(gem['gene'].unique())  # WWJ
    gem['gene_index'] = gene_cat.codes

    # use coo-matrix to compress the gem quickly
    mat = sp.coo_matrix(
        (gem['MIDcounts'], (gem['bin_index'], gem['gene_index']))).tocsr()  # 构建一个稀疏矩阵

    # extract the unspliced mRNA matrix
    un_gem = gem[gem['status'] == 'Unspliced']
    un_counts = np.append(un_gem['MIDcounts'], 0)
    un_bin = np.append(un_gem['bin_index'], gem['bin_index'].max())
    un_gene = np.append(un_gem['gene_index'], gem['gene_index'].max())
    un_mat = sp.coo_matrix((un_counts, (un_bin, un_gene))).tocsr()
    del un_gem

    # extract the spliced mRNA matrix
    sp_gem = gem[gem['status'].isin(['Spliced', 'Ambiguous'])]
    sp_counts = np.append(sp_gem['MIDcounts'], 0)
    sp_bin = np.append(sp_gem['bin_index'], gem['bin_index'].max())
    sp_gene = np.append(sp_gem['gene_index'], gem['gene_index'].max())
    sp_mat = sp.coo_matrix((sp_counts, (sp_bin, sp_gene))).tocsr()
    del sp_gem

    # set the obs list, which contain the sample infos
    obs = pd.DataFrame(index=bin_list)
    obs['bin_loc'] = bin_list
    loc = obs['bin_loc'].str.split('_', 1, expand=True).astype(int)
    loc.columns = ['cx', 'cy']
    obs = pd.merge(obs, loc, how='left', left_index=True, right_index=True)
    obs['sample'] = sample
    # set the var list, which contain the gene info
    var = pd.DataFrame(index=gene_list)
    var['gene_ID'] = gene_list
    # trans the compressed data to AnnData and save as h5ad file
    adata = ad.AnnData(mat, obs=obs, var=var, layers={
                       'spliced': sp_mat, 'unspliced': un_mat})
    # adata = ad.AnnData(mat,obs=obs,var=var,) ##wwj
    return adata


# def adataMerger(dir, max_x, max_y, gem_num, clu_list, bin, skiprows):

#     print('step2. concat gems')
#     files = os.listdir(dir)
#     files.sort()
#     pattern = re.compile(r'.+\.gem')
#     y_num = math.ceil(gem_num**0.5)
#     gem_count = 0
#     gem_list = []

#     for file in files:
#         if re.match(pattern, file):
#             gem_list.append(file)

#     gene_list = np.array([])
#     for gem_file in gem_list:
#         file_path = dir.rstrip('/\\') + '/' + gem_file
#         try:
#             gem = pd.read_csv(file_path, delimiter='\t', skiprows=skiprows)
#         except pd.errors.ParserError:
#             exit(
#                 '\033[1;31mERROR:You may need to set -s to skip some lines.\033[0m')
#         gene_cat = pd.Categorical(gem['geneID'])
#         gene_list = np.append(gene_list, gene_cat.unique())
#     gene_list = np.unique(gene_list)
#     if clu_list != None:
#         gene_list = np.intersect1d(gene_list, clu_list)

#     start_obs = pd.DataFrame(index=['starter'])
#     start_obs['bin_loc'] = np.nan
#     start_obs['sample'] = np.nan
#     start_obs['cx'] = int(0)
#     start_obs['cy'] = int(0)
#     start_var = pd.DataFrame(index=gene_list)
#     start_var['gene_ID'] = gene_list
#     start_mat = sp.csr_matrix(np.zeros((1, len(gene_list))).astype(int))

#     adata = ad.AnnData(start_mat, obs=start_obs, var=start_var)

#     max_x = max_x + bin if max_x % bin == 0 else math.ceil(max_x/bin)*bin
#     max_y = max_y + bin if max_y % bin == 0 else math.ceil(max_y/bin)*bin

#     for gem_file in gem_list:
#         print('dealing with {}'.format(gem_file))
#         file_path = dir.rstrip('/\\') + '/' + gem_file
#         try:
#             gem = pd.read_csv(file_path, delimiter='\t', skiprows=skiprows)
#         except pd.errors.ParserError:
#             exit(
#                 '\033[1;31mERROR:You may need to set -s to skip some lines.\033[0m')
#         gem['x'] = gem['x'] - gem['x'].min() + max_x * \
#             (math.floor(gem_count/y_num))
#         gem['y'] = gem['y'] - gem['y'].min() + max_y*((gem_count) % y_num)
#         if clu_list != None:
#             indices = gem['geneID'].isin(clu_list)
#             gem = gem[indices]
#         if 'Type' not in gem.columns:
#             gem['Type'] = 'S'
#         new_adata = transGemToAnnData(gem, bin, gem_file.replace('.gem', ''))
#         gem_count += 1
#         adata = ad.concat([adata, new_adata], join='outer', merge='first')

#     return adata


def workflow(args):
    """main workflow

    Args:
        args (args): args

    Returns:
        adata: adata
        sample: sample name
    """
    if args.clu_list != None:
        with open(args.clu_list) as in_file:
            clu_list = in_file.read().splitlines()
    else:
        clu_list = None

    sample = args.in_gem.split('/')[-1].replace('.gem', '')
    # sample = 'wwjtest'
    gem = gemZoomer(args.in_gem, clu_list)
    adata = transGemToAnnData(gem, args.bin, sample)

    sc.pp.filter_cells(adata, min_counts=1)
    adata.uns['gem_type'] = 'bin'
    adata.obsm['X_spatial'] = np.concatenate((adata.obs.cx, adata.obs.cy), axis=0).\
        reshape(2, -1).T.astype(np.float32)
    return adata, sample


if __name__ == '__main__':
    """main
    """
    args = getArgs()
    adata, sample = workflow(args)
    out_data = '{}_bin_{}.h5ad'.format(sample, args.bin)
    adata.write(out_data, as_dense='X')
    logging.info(
        "\033[32m All Done! Please start analysis by {} \033[0m".format(out_data))
