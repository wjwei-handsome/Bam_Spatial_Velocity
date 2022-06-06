#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :gem_format.py
@说明        : format the raw gem file
@时间        :2022/01/06 16:56:00
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''


import numpy as np
import pandas as pd

import argparse
import logging
from tqdm import tqdm

tqdm.pandas(desc="Progressing!")  # set the precess bar description

logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG)


def getArgs():
    # receive the arguments
    parser = argparse.ArgumentParser(description=__doc__, prog='gem_format.py')
    parser.add_argument('-i', action='store', dest='in_gem', type=str, required=True,
                        help='your raw gem file')
    parser.add_argument('-b', action='store', dest='binWidth', type=int, default=50,
                        help='width of bin, default = 50')
    parser.add_argument('-r', action='store', dest='rangefile', type=str,
                        required=True, help='postion range file')
    parser.add_argument('-o', action='store', dest='out_gem', type=str,
                        required=True, help='ouput gem file name')
    return parser.parse_args()

# 0. define map&apply functions


def pos_bin_mapper(x, binWidth):
    """trans pos('18183_12175') to bin('18150_12150'), default binWidth 50  

    Args:
        x ([str in pandas.colnums]): [Any]
        binWidth ([int]): [the width of bin, default 50]

    Returns:
        [str]: [x and y concat with '_' after bin]
    """
    x_num_str_before_bin = x.split('_')[0]
    y_num_str_before_bin = x.split('_')[1]
    x_num_int_after_bin = np.floor(
        int(x_num_str_before_bin)/binWidth).astype(int)*binWidth
    y_num_int_after_bin = np.floor(
        int(y_num_str_before_bin)/binWidth).astype(int)*binWidth
    new_pos = str(x_num_int_after_bin)+'_'+str(y_num_int_after_bin)
    return new_pos


def select_in_range(x):
    """judge if in range

    Args:
        x (str): [x and y concat with '_' after bin]

    Returns:
        [bool]: [boolean]
    """
    flag = x in pos_range
    return flag


def split_pos(x, idx):
    """split pos to x and y

    Args:
        x (str): [x and y concat with '_' after bin]
        idx ([int]): [0 for x, 1 for y]

    Returns:
        [int]: [x or y]
    """
    num_str = x.split('_')[idx]
    num_int = int(num_str)
    return num_int


# 1. read the raw_gem file as DataFrame

def format_gem(args):
    """mian workflow of formatting gem file

    Args:
        args (args): args object

    Returns:
        out_df: output gem DataFrame 
    """
    try:
        logging.info(
            "\033[32m Reading the raw gem file {} \033[0m".format(args.in_gem))
        gem_df = pd.read_csv(args.in_gem, sep='\t', header=None, usecols=[0, 1, 2, 3, 4],
                             names=['gene', 'pos', 'count', 'status', 'chrtype'])[:-1]
        records_num = gem_df.shape[0]
        logging.info("\033[32m Finish reading raw gem file {} for {} records \033[0m".format(
            args.in_gem, records_num))
    except:
        logging.error(
            '\033[31m Reading raw gem file Error, Please check ! \033[0m')
        exit()

    # 2. trans pos to bin

    #binWidth = 50
    logging.info("\033[32m Apply pos to bin{} \033[0m".format(args.binWidth))
    gem_df['pos'] = gem_df['pos'].progress_apply(
        pos_bin_mapper, args=(args.binWidth,))

    # 3. select records which in valid ranges

    # add the 'flag' columns which represtent if in valid rangess
    logging.info("\033[32m Select records which in valid ranges... \033[0m")
    gem_df['flag'] = gem_df['pos'].progress_map(select_in_range)

    # select DataFrame that in ranges(flag==True)
    gem_df_inrange = gem_df[gem_df['flag'] == True][['gene', 'pos', 'status']]
    select_num = gem_df_inrange.shape[0]
    select_pct = select_num/records_num
    logging.info("\033[32m Selected {} records, [{}]\033[0m".format(
        select_num, select_pct))

    # 4. get the count of special gene in special cell has special status
    logging.info(
        "\033[32m Get the count of special gene in special cell has special status \033[0m")
    out_df = gem_df_inrange.groupby(['gene', 'pos', 'status']).size(
    ).reset_index().rename(columns={0: 'MIDcounts'})

    # add 'x' and 'y' columns split by pos
    logging.info("\033[32m Split pos into x and y \033[0m")
    out_df['x'] = out_df['pos'].progress_apply(split_pos, args=(0,))
    out_df['y'] = out_df['pos'].progress_apply(split_pos, args=(1,))

    # 5. output the tsv
    out_order = ['gene', 'x', 'y', 'MIDcounts', 'status']
    #out_df[out_ouder].to_csv('./chr9.out',index=False,sep= '\t')
    return out_df[out_order]


if __name__ == '__main__':
    args = getArgs()
    logging.info("\033[32m Start to format the raw gem file {} with bin{} \033[0m".format(
        args.in_gem, args.binWidth))
    # read the ranges as array in one row
    try:
        logging.info(
            "\033[32m Reading the range file {} \033[0m".format(args.rangefile))
        pos_range = pd.read_csv(args.rangefile, header=None, names=[
                                'pos', ]).values.flatten()
    except:
        logging.error(
            '\033[31m Reading range file Error, Please check the filepath! \033[0m')
        exit()
    out_df = format_gem(args)
    logging.info(
        "\033[32m Progcess done! Start output to {}.... \033[0m".format(args.out_gem))
    out_df.to_csv(args.out_gem, index=False, sep='\t')
    logging.info(
        "\033[32m All Done! Please check gem file {} ,and trans it to adata \033[0m".format(args.out_gem))
